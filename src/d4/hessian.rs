use super::*;

/// Analytical Cartesian Hessian at fixed lattice, stored column-major (3N by 3N).
/// Hard-cutoff and image-selection boundaries must not be crossed.
pub fn hessian(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    cell: Option<(&[f64; 9], [bool; 3])>,
) -> Result<Vec<f64>, &'static str> {
    if model.ewald.is_some() {
        return Err("D4 Fourier Hessians are not supported");
    }
    validate(numbers, positions)?;
    model.validate_selection(numbers.len())?;
    validate_cutoff(model.cutoff)?;
    positions
        .len()
        .checked_mul(positions.len())
        .filter(|&count| count <= isize::MAX as usize / std::mem::size_of::<f64>())
        .ok_or("D4 Hessian dimensions overflow")?;
    let cell = cell.filter(|(_, periodic)| periodic.iter().any(|&active| active));
    let wrapped;
    let positions = if let Some((lattice, periodic)) = cell {
        wrapped = wrap_positions(positions, lattice, periodic)?;
        wrapped.as_slice()
    } else {
        positions
    };
    let (lattice, periodic) = cell.unwrap_or((&[0.0; 9], [false; 3]));
    TAPE.with(|tape| unsafe { *tape.get() = Tape::default() });
    let coordinates: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(index, &value)| Dual::variable(value, positions.len(), index))
        .collect();
    let energy = energy(
        numbers,
        &coordinates,
        charge,
        model,
        param,
        lattice,
        periodic,
    )?;
    energy.hessian()
}

fn distance2(
    coordinates: &[Dual],
    lattice: &[f64; 9],
    first: usize,
    second: usize,
    image: [i32; 3],
) -> Dual {
    let translation = lattice_translation(lattice, image);
    (0..3).fold(Dual::constant(0.0, coordinates[0].size), |sum, axis| {
        let delta = coordinates[3 * first + axis].clone() - coordinates[3 * second + axis].clone()
            + -translation[axis];
        sum + delta.clone() * delta
    })
}

pub(super) fn coordination(
    numbers: &[i32],
    coordinates: &[Dual],
    lattice: &[Dual; 9],
    lattice_values: &[f64; 9],
    periodic: [bool; 3],
    charge_model: bool,
    cutoff: f64,
) -> Vec<Dual> {
    let size = coordinates[0].size;
    let constant = |value| Dual::constant(value, size);
    let images = lattice_indices(lattice_values, periodic, cutoff);
    let positions = positions_values(coordinates);
    let mut result = vec![constant(0.0); numbers.len()];
    for first in 0..numbers.len() {
        let left = numbers[first] as usize - 1;
        for second in 0..=first {
            let right = numbers[second] as usize - 1;
            let radius = if charge_model {
                real(CHARGE_RCOV, left) + real(CHARGE_RCOV, right)
            } else {
                element(left, 0) + element(right, 0)
            };
            let factor = if charge_model {
                1.0
            } else {
                let difference = (element(left, 3) - element(right, 3)).abs();
                4.10451 * (-(difference + 19.08857).powi(2) / (2.0 * 11.28174_f64.powi(2))).exp()
            };
            for &image in &images {
                let scalar: f64 =
                    image_displacement(&positions, lattice_values, first, second, image)
                        .iter()
                        .map(|value| value * value)
                        .sum();
                if scalar < 1e-12 || scalar > cutoff * cutoff {
                    continue;
                }
                let squared = (0..3).fold(constant(0.0), |sum, axis| {
                    let delta = (0..3).fold(
                        coordinates[3 * first + axis].clone()
                            - coordinates[3 * second + axis].clone(),
                        |value, column| {
                            value - lattice[axis + 3 * column].clone() * image[column] as f64
                        },
                    );
                    sum + delta.clone() * delta
                });
                let count =
                    (((squared.sqrt() + -radius) * (-7.5 / radius)).erf() + 1.0) * (0.5 * factor);
                result[first] += count.clone();
                if first != second {
                    result[second] += count;
                }
            }
        }
    }
    if charge_model {
        result
            .into_iter()
            .map(|value| {
                constant((1.0 + 8.0_f64.exp()).ln()) - ((constant(8.0) - value).exp() + 1.0).ln()
            })
            .collect()
    } else {
        result
    }
}

fn charges(
    numbers: &[i32],
    coordinates: &[Dual],
    charge: f64,
    model: Model,
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<Vec<Dual>, &'static str> {
    let constant = |value| Dual::constant(value, coordinates[0].size);
    if let Some(values) = model.fixed_charges {
        return Ok(values.iter().map(|&value| constant(value)).collect());
    }
    if periodic.iter().any(|&active| active) {
        return periodic_dual_charges::<true>(
            numbers,
            coordinates,
            charge,
            &lattice.map(constant),
            lattice,
            periodic,
            model,
        );
    }
    if model.eeqbc {
        return eeqbc_charges(numbers, coordinates, charge, None);
    }
    if !charge.is_finite() {
        return Err("D4 total charge must be finite");
    }
    let atoms = numbers.len();
    let dimension = atoms + 1;
    let cn = coordination(
        numbers,
        coordinates,
        &lattice.map(constant),
        lattice,
        periodic,
        true,
        25.0,
    );
    let mut matrix = vec![constant(0.0); dimension * dimension];
    let mut rhs = vec![constant(0.0); dimension];
    for first in 0..atoms {
        let left = numbers[first] as usize - 1;
        rhs[first] = (cn[first].clone() + 1e-14).sqrt() * eeq(left, 2) + -eeq(left, 0);
        matrix[first * dimension + atoms] = constant(1.0);
        matrix[atoms * dimension + first] = constant(1.0);
        matrix[first * dimension + first] =
            constant(eeq(left, 1) + (2.0 / PI).sqrt() / eeq(left, 3));
        for second in 0..first {
            let right = numbers[second] as usize - 1;
            let gamma = 1.0 / (eeq(left, 3).powi(2) + eeq(right, 3).powi(2)).sqrt();
            let distance = distance2(coordinates, lattice, first, second, [0; 3]).sqrt();
            let interaction = (distance.clone() * gamma).erf() / distance;
            matrix[first * dimension + second] = interaction.clone();
            matrix[second * dimension + first] = interaction;
        }
    }
    rhs[atoms] = constant(charge);
    dual_solve(&mut matrix, &mut rhs)?;
    rhs.pop();
    Ok(rhs)
}

pub(super) fn eeq_pair(
    coordinates: &[Dual],
    lattice: &[f64; 9],
    first: usize,
    second: usize,
    gamma: f64,
    alpha: f64,
) -> Dual {
    let constant = |value| Dual::constant(value, coordinates[0].size);
    let images = closest_images(&positions_values(coordinates), lattice, first, second);
    let inverse = inverse(lattice);
    let reciprocal = std::array::from_fn(|index| inverse[(index % 3) * 3 + index / 3] * 2.0 * PI);
    let volume = determinant(lattice);
    let mut total = constant(0.0);
    for image in &images {
        for direct in fixed_indices(2, true) {
            let squared = distance2(
                coordinates,
                lattice,
                first,
                second,
                std::array::from_fn(|axis| image[axis] - direct[axis]),
            );
            if squared.value <= f64::EPSILON {
                continue;
            }
            let distance = squared.sqrt();
            total +=
                ((distance.clone() * gamma).erf() - (distance.clone() * alpha).erf()) / distance;
        }
        let translation = lattice_translation(lattice, *image);
        for wave_index in fixed_indices(2, false) {
            let wave = lattice_translation(&reciprocal, wave_index);
            let squared: f64 = wave.iter().map(|value| value * value).sum();
            let factor = (-squared / (4.0 * alpha * alpha)).exp() / squared / volume * (4.0 * PI);
            let phase = (0..3).fold(constant(0.0), |sum, axis| {
                sum + (coordinates[3 * first + axis].clone()
                    - coordinates[3 * second + axis].clone()
                    + -translation[axis])
                    * wave[axis]
            });
            total += phase.clone().unary_second(
                phase.value.cos(),
                -phase.value.sin(),
                -phase.value.cos(),
            ) * factor;
        }
    }
    total / images.len() as f64
}

fn weights(atom: usize, other: usize, cn: &Dual, charge: &Dual, model: Model) -> Vec<Dual> {
    let constant = |value| Dual::constant(value, cn.size);
    let width = if model.d4s {
        real(PAIR_WEIGHTS, other + ELEMENTS * atom)
    } else {
        model.wf
    };
    let mut values: Vec<_> = (0..int(NREF, atom) as usize)
        .map(|reference_index| {
            let delta = cn.clone() + -reference(REFCN, atom, reference_index);
            (1..=int(NGW, reference_index + REFERENCES * atom)).fold(
                constant(0.0),
                |sum, gaussian| {
                    sum + (delta.clone() * delta.clone() * (-width * gaussian as f64)).exp()
                },
            )
        })
        .collect();
    let norm = values
        .iter()
        .fold(constant(0.0), |sum, value| sum + value.clone());
    let effective = charge.clone() + element(atom, 2);
    for (index, value) in values.iter_mut().enumerate() {
        let reference_charge =
            reference(if model.eeqbc { EEQBC_REFQ } else { REFQ }, atom, index) + element(atom, 2);
        let zeta = if effective.value < 0.0 {
            constant(model.ga.exp())
        } else {
            ((constant(1.0)
                - ((constant(1.0) - constant(reference_charge) / effective.clone())
                    * (model.gc * element(atom, 4)))
                .exp())
                * model.ga)
                .exp()
        };
        *value = value.clone() / norm.clone() * zeta;
    }
    values
}

fn coefficients(numbers: &[i32], cn: &[Dual], charges: &[Dual], model: Model) -> Vec<Dual> {
    let atoms = numbers.len();
    let mut result = vec![Dual::constant(0.0, cn[0].size); atoms * atoms];
    for first in 0..atoms {
        let left = numbers[first] as usize - 1;
        for second in 0..=first {
            let right = numbers[second] as usize - 1;
            let left_weights = weights(left, right, &cn[first], &charges[first], model);
            let right_weights = weights(right, left, &cn[second], &charges[second], model);
            let mut value = Dual::constant(0.0, cn[0].size);
            for (left_ref, left_weight) in left_weights.iter().enumerate() {
                for (right_ref, right_weight) in right_weights.iter().enumerate() {
                    value += left_weight.clone()
                        * right_weight.clone()
                        * reference_c6(left, left_ref, right, right_ref, model);
                }
            }
            result[first * atoms + second] = value.clone();
            result[second * atoms + first] = value;
        }
    }
    result
}

fn switch(squared: Dual, cutoff: f64, width: f64) -> Dual {
    let constant = |value| Dual::constant(value, squared.size);
    if width <= 0.0 || cutoff <= 0.0 || squared.value.sqrt() <= cutoff - width.min(cutoff) {
        return constant(1.0);
    }
    if squared.value >= cutoff * cutoff {
        return constant(0.0);
    }
    let fraction = (constant(cutoff) - squared.sqrt()) / width.min(cutoff);
    fraction.clone().powf(3.0) * (fraction.clone() * (fraction * 6.0 + -15.0) + 10.0)
}

fn energy(
    numbers: &[i32],
    coordinates: &[Dual],
    charge: f64,
    model: Model,
    param: Param,
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<Dual, &'static str> {
    let constant = |value| Dual::constant(value, coordinates[0].size);
    let atoms = numbers.len();
    let cutoff = model.cutoff;
    let positions = positions_values(coordinates);
    let squared_value = |first, second, image| {
        image_displacement(&positions, lattice, first, second, image)
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
    };
    let cn = coordination(
        numbers,
        coordinates,
        &lattice.map(constant),
        lattice,
        periodic,
        false,
        cutoff.cn,
    );
    let charges = charges(numbers, coordinates, charge, model, lattice, periodic)?;
    let c6 = coefficients(numbers, &cn, &charges, model);
    let mut total = constant(0.0);
    let images = lattice_indices(lattice, periodic, cutoff.disp2);
    for first in 0..atoms {
        for second in 0..=first {
            if !model.owns_pair(first, second) {
                continue;
            }
            let rrij = 3.0
                * element(numbers[first] as usize - 1, 1)
                * element(numbers[second] as usize - 1, 1);
            let radius = param.a1 * rrij.sqrt() + param.a2;
            let weight = if first == second { 0.5 } else { 1.0 };
            for &image in &images {
                let scalar = squared_value(first, second, image);
                if scalar < f64::EPSILON || scalar > cutoff.disp2.powi(2) {
                    continue;
                }
                let squared = distance2(coordinates, lattice, first, second, image);
                let potential = constant(param.s6) / (squared.clone().powf(3.0) + radius.powi(6))
                    + constant(param.s8 * rrij) / (squared.clone().powf(4.0) + radius.powi(8));
                total = total
                    - c6[first * atoms + second].clone()
                        * potential
                        * switch(squared, cutoff.disp2, cutoff.width2)
                        * weight;
            }
        }
    }
    if param.s9.abs() < f64::EPSILON {
        return Ok(total);
    }
    let c6 = coefficients(numbers, &cn, &vec![constant(0.0); atoms], model);
    let images = lattice_indices(lattice, periodic, cutoff.disp3);
    for first in 0..atoms {
        let neighbors: Vec<Vec<_>> = (0..=first)
            .map(|second| {
                images
                    .iter()
                    .copied()
                    .filter(|&image| {
                        let squared = squared_value(first, second, image);
                        squared >= f64::EPSILON && squared <= cutoff.disp3.powi(2)
                    })
                    .collect()
            })
            .collect();
        for second in 0..=first {
            if !model.owns_pair(first, second) {
                continue;
            }
            for third in 0..=second {
                if !model.active(third) {
                    continue;
                }
                let pairs = [(first, second), (first, third), (second, third)];
                let radius: f64 = pairs
                    .iter()
                    .map(|&(left, right)| {
                        param.a1
                            * (3.0
                                * element(numbers[left] as usize - 1, 1)
                                * element(numbers[right] as usize - 1, 1))
                            .sqrt()
                            + param.a2
                    })
                    .product();
                let weight = if first == third {
                    1.0 / 6.0
                } else if first == second || second == third {
                    0.5
                } else {
                    1.0
                };
                for &left_image in &neighbors[second] {
                    let left = distance2(coordinates, lattice, first, second, left_image);
                    for &right_image in &neighbors[third] {
                        let opposite_image =
                            std::array::from_fn(|axis| right_image[axis] - left_image[axis]);
                        let scalar = squared_value(second, third, opposite_image);
                        if scalar < f64::EPSILON || scalar > cutoff.disp3.powi(2) {
                            continue;
                        }
                        let right = distance2(coordinates, lattice, first, third, right_image);
                        let opposite =
                            distance2(coordinates, lattice, second, third, opposite_image);
                        let product2 = left.clone() * right.clone() * opposite.clone();
                        let angular = (left.clone() + opposite.clone() - right.clone())
                            * (left.clone() + right.clone() - opposite.clone())
                            * (right.clone() + opposite.clone() - left.clone())
                            * 0.375
                            / product2.clone().powf(2.5)
                            + constant(1.0) / product2.clone().powf(1.5);
                        let damping = constant(1.0)
                            / ((constant(radius) / product2.sqrt()).powf(param.alpha / 3.0) * 6.0
                                + 1.0);
                        let c9 = (c6[first * atoms + second].clone()
                            * c6[first * atoms + third].clone()
                            * c6[second * atoms + third].clone())
                        .sqrt()
                            * (param.s9 * weight);
                        total += c9
                            * angular
                            * damping
                            * switch(left.clone(), cutoff.disp3, cutoff.width3)
                            * switch(right, cutoff.disp3, cutoff.width3)
                            * switch(opposite, cutoff.disp3, cutoff.width3);
                    }
                }
            }
        }
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn molecular_hessians_match_gradient_differences() {
        check_hessians(None);
    }

    #[test]
    fn periodic_hessians_match_gradient_differences() {
        let lattice = [6.1, 0.3, 0.2, 1.1, 7.2, 0.5, 0.4, 0.6, 8.0];
        for periodic in [[true, false, false], [true, false, true], [true; 3]] {
            check_hessians(Some((&lattice, periodic)));
        }
    }

    fn check_hessians(cell: Option<(&[f64; 9], [bool; 3])>) {
        let numbers = [6, 8, 1];
        let positions = [1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.5, 3.0, 1.5];
        let fixed = [0.2, -0.4, 0.1];
        let param = load_param("pbe", true).unwrap();
        for param in [
            param,
            Param { s9: 0.0, ..param },
            Param {
                s6: 0.0,
                s8: 0.0,
                ..param
            },
        ] {
            for mut base in [
                Model::D4,
                Model::D4S,
                Model::custom(false, 2.0, 1.0, 4.0).unwrap(),
                Model::custom(true, 2.0, 1.0, 6.0).unwrap(),
            ] {
                base.set_cutoff(Cutoff {
                    cn: 7.0,
                    disp2: 8.0,
                    disp3: 6.0,
                    width2: 1.0,
                    width3: 1.0,
                })
                .unwrap();
                let mut eeqbc_model = base;
                eeqbc_model.set_charge_model(1).unwrap();
                for model in [
                    base,
                    eeqbc_model,
                    base.with_fixed_charges(Some(&fixed)).unwrap(),
                ] {
                    let actual = hessian(&numbers, &positions, 0.25, model, param, cell).unwrap();
                    let evaluate = |positions: &[f64]| {
                        if let Some((lattice, periodic)) = cell {
                            periodic_dispersion(
                                &numbers, positions, 0.25, model, param, lattice, periodic,
                            )
                            .unwrap()
                        } else {
                            dispersion(&numbers, positions, 0.25, model, param).unwrap()
                        }
                    };
                    for column in 0..9 {
                        for axis in 0..3 {
                            assert!(
                                (0..3)
                                    .map(|atom| actual[column * 9 + 3 * atom + axis])
                                    .sum::<f64>()
                                    .abs()
                                    < 1e-12
                            );
                        }
                        let mut displaced = positions;
                        displaced[column] += 1e-4;
                        let plus = evaluate(&displaced);
                        displaced[column] -= 2e-4;
                        let minus = evaluate(&displaced);
                        for row in 0..9 {
                            let expected = (plus.gradient[row] - minus.gradient[row]) / 2e-4;
                            assert!(
                                (actual[column * 9 + row] - expected).abs() < 1e-9,
                                "row={row} column={column}: {} != {expected}",
                                actual[column * 9 + row]
                            );
                            assert!(
                                (actual[column * 9 + row] - actual[row * 9 + column]).abs() < 1e-12
                            );
                        }
                    }
                }
            }
        }
    }
}
