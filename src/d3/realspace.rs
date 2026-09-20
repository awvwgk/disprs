use super::*;
use crate::dual::Dual2;
use crate::geometry::{lattice_repetitions, lattice_translation, periodic_reciprocal};

#[allow(clippy::type_complexity)]
#[allow(
    clippy::cast_possible_truncation,
    reason = "periodic_reciprocal rejects degenerate cells, so the image index stays finite"
)]
pub(super) fn prepare(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    periodic: [bool; 3],
    cn_cutoff: f64,
    radius: f64,
) -> Result<(Prepared, Vec<Vec<(usize, [f64; 3], f64)>>), &'static str> {
    let reciprocal = periodic_reciprocal(lattice, periodic)?;
    let (elements, species, _, _) = prepare_with_cutoff(numbers, positions, 0.0)?;
    let atoms = numbers.len();
    let repetitions = lattice_repetitions(lattice, periodic, radius);
    let repetitions: [i32; 3] =
        std::array::from_fn(|axis| repetitions[axis] + i32::from(periodic[axis]));
    let mut neighbors = vec![Vec::new(); atoms];
    let mut coordination = vec![0.0; atoms];
    for first in 0..atoms {
        for second in 0..atoms {
            let base = displacement(positions, first, second);
            let center: [i32; 3] = std::array::from_fn(|column| {
                (0..3)
                    .map(|axis| base[axis] * reciprocal[axis + 3 * column])
                    .sum::<f64>()
                    .round() as i32
            });
            for first_image in -repetitions[0]..=repetitions[0] {
                for second_image in -repetitions[1]..=repetitions[1] {
                    for third_image in -repetitions[2]..=repetitions[2] {
                        let image = [
                            center[0] + first_image,
                            center[1] + second_image,
                            center[2] + third_image,
                        ];
                        let shift = lattice_translation(lattice, image);
                        let vector: [f64; 3] = std::array::from_fn(|axis| base[axis] - shift[axis]);
                        let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
                        if distance2 < 1.0e-12 || distance2 > radius * radius {
                            continue;
                        }
                        neighbors[first].push((second, vector, distance2));
                        if distance2 <= cn_cutoff * cn_cutoff {
                            let radius = value(COVALENT_RADII, elements[first])
                                + value(COVALENT_RADII, elements[second]);
                            coordination[first] +=
                                1.0 / (1.0 + (-16.0 * (radius / distance2.sqrt() - 1.0)).exp());
                        }
                    }
                }
            }
        }
    }
    let weights = elements
        .iter()
        .zip(&coordination)
        .map(|(&element, &cn)| weights(element, cn))
        .collect();
    Ok(((elements, species, coordination, weights), neighbors))
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn periodic_realspace(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    periodic: [bool; 3],
    damping: Damping,
    atm: Option<Atm>,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(PeriodicResult, Vec<f64>, Vec<f64>), &'static str> {
    evaluate(
        numbers,
        positions,
        lattice,
        periodic,
        damping,
        atm,
        cutoff,
        ghosts,
        partition,
        None,
        Model::D3,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn periodic_realspace_hessian(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    periodic: [bool; 3],
    damping: Damping,
    atm: Option<Atm>,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>), &'static str> {
    let mut hessian = Vec::new();
    let result = evaluate(
        numbers,
        positions,
        lattice,
        periodic,
        damping,
        atm,
        cutoff,
        ghosts,
        partition,
        Some(&mut hessian),
        Model::D3,
    )?;
    Ok((result.0.energy, hessian))
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn evaluate(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    periodic: [bool; 3],
    damping: Damping,
    atm: Option<Atm>,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
    hessian: Option<&mut Vec<f64>>,
    model: Model,
) -> Result<(PeriodicResult, Vec<f64>, Vec<f64>), &'static str> {
    model.validate(numbers)?;
    if numbers.is_empty() || positions.iter().any(|value| !value.is_finite()) {
        return Err("D3 requires a nonempty finite structure");
    }
    if [
        cutoff.cn,
        cutoff.disp2,
        cutoff.disp3,
        cutoff.width2,
        cutoff.width3,
    ]
    .iter()
    .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err("D3 cutoffs must be finite and nonnegative");
    }
    validate_ghosts(numbers, ghosts)?;
    let atoms = numbers.len();
    let radius = cutoff
        .cn
        .max(cutoff.disp2)
        .max(if atm.is_some() { cutoff.disp3 } else { 0.0 });
    let ((elements, species, coordination, weights), neighbors) =
        prepare(numbers, positions, lattice, periodic, cutoff.cn, radius)?;
    let derivatives: Vec<_> = elements
        .iter()
        .zip(&coordination)
        .zip(&weights)
        .map(|((&element, &cn), weights)| weight_derivatives(element, cn, weights))
        .collect();
    let mut cn_dual = Vec::new();
    if hessian.is_some() {
        for first in 0..atoms {
            let mut total = Dual2::constant(0.0, positions.len());
            for &(second, vector, distance2) in &neighbors[first] {
                if distance2 > cutoff.cn * cutoff.cn {
                    continue;
                }
                let distance =
                    distance_dual(positions.len(), first, second, vector, distance2).sqrt();
                let radius = value(COVALENT_RADII, elements[first])
                    + value(COVALENT_RADII, elements[second]);
                total += 1.0 / (((radius / distance - 1.0) * -16.0).exp() + 1.0);
            }
            cn_dual.push(total);
        }
    }
    let mut coefficients = vec![(0.0, 0.0, 0.0); atoms * atoms];
    let mut coefficient_dual = Vec::new();
    for first in 0..atoms {
        for second in 0..atoms {
            let smooth = (model == Model::D3S).then(|| {
                smooth::pair_coefficients(
                    elements[first],
                    elements[second],
                    coordination[first],
                    coordination[second],
                )
            });
            coefficients[first * atoms + second] = if let Some(partial) = smooth {
                (partial[0], partial[1], partial[2])
            } else {
                atomic_c6_derivatives(
                    elements[first],
                    elements[second],
                    &weights[first],
                    &weights[second],
                    &derivatives[first],
                    &derivatives[second],
                )
            };
            if hessian.is_some() {
                let partial = smooth.unwrap_or_else(|| {
                    atomic_c6_second_derivatives(
                        elements[first],
                        elements[second],
                        &weights[first],
                        &weights[second],
                        &derivatives[first],
                        &derivatives[second],
                        &weight_second_derivatives(
                            elements[first],
                            coordination[first],
                            &weights[first],
                        ),
                        &weight_second_derivatives(
                            elements[second],
                            coordination[second],
                            &weights[second],
                        ),
                    )
                });
                coefficient_dual.push(compose(
                    partial[0],
                    &[partial[1], partial[2]],
                    &[partial[3], partial[4], partial[4], partial[5]],
                    &[&cn_dual[first], &cn_dual[second]],
                ));
            }
        }
    }
    let mut total_dual = Dual2::constant(
        0.0,
        if hessian.is_some() {
            positions.len()
        } else {
            0
        },
    );
    let mut result = PeriodicResult {
        energy: 0.0,
        gradient: vec![0.0; positions.len()],
        virial: [0.0; 9],
    };
    let mut pair2 = vec![0.0; atoms * atoms];
    let mut pair3 = vec![0.0; atoms * atoms];
    let mut dedcn = vec![0.0; atoms];
    let mut energy_correction = 0.0;
    for first in 0..atoms {
        if !partition.owns_index(first) || ghosts.get(first) == Some(&true) {
            continue;
        }
        for &(second, vector, distance2) in &neighbors[first] {
            if ghosts.get(second) == Some(&true) || distance2 > cutoff.disp2 * cutoff.disp2 {
                continue;
            }
            let (c6, dc6_first, dc6_second) = coefficients[first * atoms + second];
            let (energy, radial, coefficient) = kernel(
                damping,
                elements[first],
                elements[second],
                species[first],
                species[second],
                distance2,
                c6,
            );
            let (switch, dswitch) = smooth_cutoff(distance2, cutoff.disp2, cutoff.width2);
            if hessian.is_some() {
                let distance = distance_dual(positions.len(), first, second, vector, distance2);
                let partial = kernel_second_derivatives(
                    damping,
                    elements[first],
                    elements[second],
                    species[first],
                    species[second],
                    distance2,
                    c6,
                );
                let pair = compose(
                    partial[0],
                    &[partial[1], partial[2]],
                    &[partial[3], partial[4], partial[4], partial[5]],
                    &[&distance, &coefficient_dual[first * atoms + second]],
                );
                let switch = smooth_cutoff_second(distance2, cutoff.disp2, cutoff.width2);
                total_dual +=
                    pair * compose(switch[0], &[switch[1]], &[switch[2]], &[&distance]) * 0.5;
            }
            let increment = 0.5 * switch * energy - energy_correction;
            let total = result.energy + increment;
            energy_correction = (total - result.energy) - increment;
            result.energy = total;
            pair2[first * atoms + second] += 0.25 * switch * energy;
            pair2[second * atoms + first] += 0.25 * switch * energy;
            dedcn[first] += 0.5 * switch * coefficient * dc6_first;
            dedcn[second] += 0.5 * switch * coefficient * dc6_second;
            add_derivative(
                &mut result,
                first,
                second,
                vector,
                0.5 * (switch * radial + dswitch * energy),
            );
        }
        if let Some(atm) = atm {
            let near: Vec<_> = neighbors[first]
                .iter()
                .filter(|&&(atom, _, distance2)| {
                    distance2 <= cutoff.disp3 * cutoff.disp3 && ghosts.get(atom) != Some(&true)
                })
                .collect();
            for (index, &&(second, vector_second, distance_second)) in near.iter().enumerate() {
                for &&(third, vector_third, distance_third) in &near[..index] {
                    let vector_pair =
                        std::array::from_fn(|axis| vector_third[axis] - vector_second[axis]);
                    let distance_pair = vector_pair.iter().map(|value| value * value).sum::<f64>();
                    if distance_pair < 1.0e-12 || distance_pair > cutoff.disp3 * cutoff.disp3 {
                        continue;
                    }
                    let pairs = [(first, second), (first, third), (second, third)];
                    let values = pairs.map(|(left, right)| coefficients[left * atoms + right].0);
                    let radii = pairs.map(|(left, right)| {
                        value(VDW_RADII, pair_index(elements[left], elements[right]))
                    });
                    let (energy, radial, coefficient) = atm_triplet_derivatives(
                        [distance_second, distance_third, distance_pair],
                        values,
                        radii,
                        atm,
                        cutoff,
                    );
                    if hessian.is_some() {
                        let distances = [distance_second, distance_third, distance_pair];
                        let vectors = [vector_second, vector_third, vector_pair];
                        let arguments: [Dual2; 3] = std::array::from_fn(|index| {
                            distance_dual(
                                positions.len(),
                                pairs[index].0,
                                pairs[index].1,
                                vectors[index],
                                distances[index],
                            )
                        });
                        let local = local_atm_partials(distances, values, radii, atm, cutoff);
                        total_dual += compose(
                            local.value,
                            &local.gradient,
                            &local.hessian,
                            &[
                                &arguments[0],
                                &arguments[1],
                                &arguments[2],
                                &coefficient_dual[first * atoms + second],
                                &coefficient_dual[first * atoms + third],
                                &coefficient_dual[second * atoms + third],
                            ],
                        ) / 3.0;
                    }
                    let increment = energy / 3.0 - energy_correction;
                    let total = result.energy + increment;
                    energy_correction = (total - result.energy) - increment;
                    result.energy = total;
                    for (pair, ((left, right), vector)) in pairs
                        .iter()
                        .copied()
                        .zip([vector_second, vector_third, vector_pair])
                        .enumerate()
                    {
                        pair3[left * atoms + right] += energy / 18.0;
                        pair3[right * atoms + left] += energy / 18.0;
                        let (_, dc6_left, dc6_right) = coefficients[left * atoms + right];
                        dedcn[left] += coefficient[pair] * dc6_left / 3.0;
                        dedcn[right] += coefficient[pair] * dc6_right / 3.0;
                        add_derivative(&mut result, left, right, vector, radial[pair] / 3.0);
                    }
                }
            }
        }
    }
    for first in 0..atoms {
        for &(second, vector, distance2) in &neighbors[first] {
            if distance2 > cutoff.cn * cutoff.cn {
                continue;
            }
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (1.0 + (-16.0 * (radius / distance2.sqrt() - 1.0)).exp());
            let radial = -8.0 * radius * count * (1.0 - count) / distance2.powf(1.5);
            add_derivative(&mut result, first, second, vector, dedcn[first] * radial);
        }
    }
    if let Some(hessian) = hessian {
        *hessian = total_dual.hessian;
    }
    Ok((result, pair2, pair3))
}

fn distance_dual(
    size: usize,
    first: usize,
    second: usize,
    vector: [f64; 3],
    distance2: f64,
) -> Dual2 {
    let mut result = Dual2::constant(distance2, size);
    for (axis, &component) in vector.iter().enumerate() {
        result.gradient[3 * first + axis] += 2.0 * component;
        result.gradient[3 * second + axis] -= 2.0 * component;
        add_hessian_block(&mut result.hessian, size, first, second, axis, axis, 2.0);
    }
    result
}

fn compose(value: f64, first: &[f64], second: &[f64], arguments: &[&Dual2]) -> Dual2 {
    let size = arguments[0].gradient.len();
    let mut result = Dual2::constant(value, size);
    for (index, argument) in arguments.iter().enumerate() {
        for row in 0..size {
            result.gradient[row] += first[index] * argument.gradient[row];
            for column in 0..size {
                result.hessian[row * size + column] +=
                    first[index] * argument.hessian[row * size + column];
                for (other, other_argument) in arguments.iter().enumerate() {
                    result.hessian[row * size + column] += second[index * arguments.len() + other]
                        * argument.gradient[row]
                        * other_argument.gradient[column];
                }
            }
        }
    }
    result
}

fn add_derivative(
    result: &mut PeriodicResult,
    first: usize,
    second: usize,
    vector: [f64; 3],
    radial: f64,
) {
    for axis in 0..3 {
        let component = 2.0 * radial * vector[axis];
        result.gradient[3 * first + axis] += component;
        result.gradient[3 * second + axis] -= component;
        for (other, &value) in vector.iter().enumerate() {
            result.virial[axis + 3 * other] += component * value;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directional_realspace_hessians() {
        let numbers = [6, 8];
        let positions = [0.2, 0.3, 0.4, 2.0, 0.5, 0.6];
        let lattice = [7.0, 0.1, 0.0, 0.4, 8.0, 0.2, 0.0, 0.3, 9.0];
        let cutoff = RealspaceCutoff {
            cn: 9.0,
            disp2: 10.0,
            disp3: 8.0,
            width2: 2.0,
            width3: 2.0,
        };
        for periodic in [[true, false, false], [true, true, false], [true; 3]] {
            for kind in 0..=6 {
                let damping = if kind == 6 {
                    Damping::Z {
                        s6: 1.0,
                        s8: 1.0,
                        a1: 200770.0,
                    }
                } else {
                    load_named("pbe", kind, true).unwrap()
                };
                let atm = Some(Atm {
                    s9: 1.0,
                    alpha: 16.0,
                });
                let (energy, hessian) = periodic_realspace_hessian(
                    &numbers,
                    &positions,
                    &lattice,
                    periodic,
                    damping,
                    atm,
                    cutoff,
                    &[],
                    WorkPartition::SERIAL,
                )
                .unwrap();
                let evaluate = |positions: &[f64]| {
                    periodic_realspace(
                        &numbers,
                        positions,
                        &lattice,
                        periodic,
                        damping,
                        atm,
                        cutoff,
                        &[],
                        WorkPartition::SERIAL,
                    )
                    .unwrap()
                    .0
                };
                assert!((energy - evaluate(&positions).energy).abs() < 1.0e-14);
                for column in 0..positions.len() {
                    let mut shifted = positions;
                    shifted[column] += 1.0e-5;
                    let plus = evaluate(&shifted).gradient;
                    shifted[column] -= 2.0e-5;
                    let minus = evaluate(&shifted).gradient;
                    for row in 0..positions.len() {
                        assert!(
                            (hessian[row * positions.len() + column]
                                - (plus[row] - minus[row]) / 2.0e-5)
                                .abs()
                                < 1.0e-8,
                            "kind {kind} row {row} column {column}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn directional_realspace_derivatives_and_pair_sums() {
        let numbers = [6, 8, 1];
        let positions = [0.1, 0.3, 0.2, 2.3, 0.6, 0.5, 0.5, 2.1, 0.8];
        let lattice = [7.0, 0.2, 0.0, 1.0, 8.0, 0.3, 0.4, 0.7, 9.0];
        let cutoff = RealspaceCutoff {
            cn: 10.0,
            disp2: 12.0,
            disp3: 9.0,
            width2: 2.0,
            width3: 2.0,
        };
        for periodic in [[true, false, false], [true, false, true], [true; 3]] {
            for kind in 0..=6 {
                let damping = if kind == 6 {
                    Damping::Z {
                        s6: 1.0,
                        s8: 1.0,
                        a1: 200770.0,
                    }
                } else {
                    load_named("pbe", kind, true).unwrap()
                };
                let atm = Some(Atm {
                    s9: 1.0,
                    alpha: 16.0,
                });
                let evaluate = |positions: &[f64], lattice: &[f64; 9]| {
                    periodic_realspace(
                        &numbers,
                        positions,
                        lattice,
                        periodic,
                        damping,
                        atm,
                        cutoff,
                        &[],
                        WorkPartition::SERIAL,
                    )
                    .unwrap()
                };
                let (result, pair2, pair3) = evaluate(&positions, &lattice);
                assert!((result.energy - pair2.iter().chain(&pair3).sum::<f64>()).abs() < 1.0e-13);
                let mut translated = positions;
                for axis in 0..3 {
                    translated[axis] += 13.0 * lattice[axis];
                }
                assert!((evaluate(&translated, &lattice).0.energy - result.energy).abs() < 1.0e-13);
                for ghosts in [&[][..], &[false, true, false][..]] {
                    let calculate = |partition| {
                        periodic_realspace(
                            &numbers, &positions, &lattice, periodic, damping, atm, cutoff, ghosts,
                            partition,
                        )
                        .unwrap()
                    };
                    let (serial, serial2, serial3) = calculate(WorkPartition::SERIAL);
                    let (left, left2, left3) = calculate(WorkPartition::new(0, 2).unwrap());
                    let (right, right2, right3) = calculate(WorkPartition::new(1, 2).unwrap());
                    assert!((serial.energy - left.energy - right.energy).abs() < 1.0e-13);
                    for ((expected, first), second) in serial
                        .gradient
                        .iter()
                        .chain(&serial.virial)
                        .chain(&serial2)
                        .chain(&serial3)
                        .zip(
                            left.gradient
                                .iter()
                                .chain(&left.virial)
                                .chain(&left2)
                                .chain(&left3),
                        )
                        .zip(
                            right
                                .gradient
                                .iter()
                                .chain(&right.virial)
                                .chain(&right2)
                                .chain(&right3),
                        )
                    {
                        assert!((expected - first - second).abs() < 1.0e-13);
                    }
                }
                for coordinate in 0..positions.len() {
                    let mut shifted = positions;
                    shifted[coordinate] += 1.0e-5;
                    let plus = evaluate(&shifted, &lattice).0.energy;
                    shifted[coordinate] -= 2.0e-5;
                    let minus = evaluate(&shifted, &lattice).0.energy;
                    assert!(
                        (result.gradient[coordinate] - (plus - minus) / 2.0e-5).abs() < 1.0e-8,
                        "kind {kind}, coordinate {coordinate}"
                    );
                }
                for component in 0..9 {
                    let strained = |delta| {
                        let mut shifted = positions;
                        let mut cell = lattice;
                        for atom in 0..numbers.len() {
                            shifted[3 * atom + component % 3] +=
                                delta * positions[3 * atom + component / 3];
                        }
                        for column in 0..3 {
                            cell[component % 3 + 3 * column] +=
                                delta * lattice[component / 3 + 3 * column];
                        }
                        evaluate(&shifted, &cell).0.energy
                    };
                    assert!(
                        (result.virial[component]
                            - (strained(1.0e-5) - strained(-1.0e-5)) / 2.0e-5)
                            .abs()
                            < 1.0e-8,
                        "kind {kind}, strain {component}"
                    );
                }
            }
        }
    }
}
