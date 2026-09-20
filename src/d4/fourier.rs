use super::*;
use crate::d3::fourier::{factorize, reciprocal_derivatives, FourierTerm};

#[allow(clippy::too_many_arguments)]
pub(super) fn dispersion(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    lattice: &[f64; 9],
    periodic: [bool; 3],
    config: EwaldConfig,
) -> Result<Dispersion, &'static str> {
    validate(numbers, positions)?;
    model.validate_selection(numbers.len())?;
    if !periodic.iter().all(|&active| active) {
        return Err("D4 Fourier dispersion requires full 3D periodicity");
    }
    if param.s9.abs() >= f64::EPSILON {
        return Err("D4 Fourier dispersion does not support ATM");
    }
    if model.cutoff.width2 != 0.0 {
        return Err("D4 Fourier dispersion does not support real-space pair smoothing");
    }
    if [
        param.s6,
        param.s8,
        param.s9,
        param.a1,
        param.a2,
        param.alpha,
    ]
    .iter()
    .any(|value| !value.is_finite())
    {
        return Err("D4 damping parameters must be finite");
    }
    let positions = wrap_positions(positions, lattice, periodic)?;
    let atoms = numbers.len();
    let size = positions.len() + 9;
    let response = charge_response(
        numbers,
        &positions,
        charge,
        model,
        Some((lattice, periodic)),
        true,
    )?;
    TAPE.with(|tape| unsafe { *tape.get() = Tape::default() });
    let coordinates: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(index, &value)| Dual::variable(value, size, index))
        .collect();
    let dual_lattice =
        std::array::from_fn(|index| Dual::variable(lattice[index], size, positions.len() + index));
    let cn = periodic_dual_coordination(
        numbers,
        &coordinates,
        &dual_lattice,
        lattice,
        periodic,
        false,
        model.cutoff.cn,
    );
    let cn_gradients: Vec<_> = cn.iter().map(Dual::gradient).collect();
    let elements: Vec<_> = numbers.iter().map(|&number| number as usize - 1).collect();
    let mut species = Vec::new();
    for &element in &elements {
        if !species.contains(&element) {
            species.push(element);
        }
    }
    let mut result = Dispersion {
        energy: 0.0,
        gradient: vec![0.0; positions.len()],
        virial: [0.0; 9],
    };
    let mut dedcn = vec![0.0; atoms];
    let mut dedq = vec![0.0; atoms];
    for (left_index, &left) in species.iter().enumerate() {
        for &right in &species[..=left_index] {
            let pair_species = if left == right {
                vec![left]
            } else {
                vec![left, right]
            };
            let selected: Vec<_> = elements
                .iter()
                .enumerate()
                .filter_map(|(atom, &element)| pair_species.contains(&element).then_some(atom))
                .collect();
            let atom_species: Vec<_> = selected
                .iter()
                .map(|&atom| usize::from(elements[atom] != left))
                .collect();
            let local_positions: Vec<_> = selected
                .iter()
                .flat_map(|&atom| positions[3 * atom..3 * atom + 3].iter().copied())
                .collect();
            let counts: Vec<_> = pair_species
                .iter()
                .map(|&element| int(NREF, element) as usize)
                .collect();
            let dimension: usize = counts.iter().sum();
            let mut matrix = vec![0.0; dimension * dimension];
            for first in 0..counts[0] {
                let offset = if left == right { 0 } else { counts[0] };
                for second in 0..*counts.last().unwrap() {
                    let value = reference_c6(left, first, right, second, model);
                    matrix[first * dimension + offset + second] = value;
                    matrix[(offset + second) * dimension + first] = value;
                }
            }
            let (eigenvalues, vectors) =
                factorize(matrix, dimension, config.rank, config.tolerance);
            let rank = eigenvalues.len();
            let mut factors = vec![0.0; rank * selected.len()];
            let mut partials = vec![[0.0; 2]; factors.len()];
            for (local, &atom) in selected.iter().enumerate() {
                if !model.active(atom) {
                    continue;
                }
                let element = elements[atom];
                let other = if element == left { right } else { left };
                let weights = weight_partials(
                    element,
                    other,
                    cn[atom].value,
                    response.charges[atom],
                    model,
                );
                let offset = if element == left { 0 } else { counts[0] };
                for term in 0..rank {
                    for (reference, weight) in
                        weights.iter().enumerate().take(int(NREF, element) as usize)
                    {
                        let vector = vectors[term * dimension + offset + reference];
                        factors[term * selected.len() + local] += weight[0] * vector;
                        partials[term * selected.len() + local][0] += weight[1] * vector;
                        partials[term * selected.len() + local][1] += weight[2] * vector;
                    }
                }
            }
            let rrij = 3.0 * element(left, 1) * element(right, 1);
            let radius = param.a1 * rrij.sqrt() + param.a2;
            if !radius.is_finite() || radius <= 0.0 {
                return Err("D4 Fourier damping radius must be finite and positive");
            }
            let terms = [
                FourierTerm {
                    prefactor: param.s6,
                    numerator: 0,
                    denominator: 6,
                    radius,
                },
                FourierTerm {
                    prefactor: param.s8 * rrij,
                    numerator: 0,
                    denominator: 8,
                    radius,
                },
            ];
            let zero = terms.map(|term| FourierTerm {
                prefactor: 0.0,
                ..term
            });
            let kernels = if left == right {
                vec![terms]
            } else {
                vec![zero, terms, terms, zero]
            };
            let kcut = if config.kcut > 0.0 {
                config.kcut
            } else {
                terms
                    .iter()
                    .map(|term| -1e-8_f64.ln() / (radius * (PI / term.denominator as f64).sin()))
                    .fold(0.0, f64::max)
            };
            let (local_result, adjoints) = reciprocal_derivatives(
                &local_positions,
                lattice,
                &pair_species,
                &atom_species,
                &kernels,
                &eigenvalues,
                &factors,
                EwaldConfig { kcut, ..config },
                model.partition,
            )?;
            result.energy += local_result.energy;
            for (total, value) in result.virial.iter_mut().zip(local_result.virial) {
                *total += value;
            }
            for (local, &atom) in selected.iter().enumerate() {
                for axis in 0..3 {
                    result.gradient[3 * atom + axis] += local_result.gradient[3 * local + axis];
                }
                for term in 0..rank {
                    let index = term * selected.len() + local;
                    dedcn[atom] += adjoints[index] * partials[index][0];
                    dedq[atom] += adjoints[index] * partials[index][1];
                }
            }
        }
    }
    for atom in 0..atoms {
        for (coordinate, value) in result.gradient.iter_mut().enumerate() {
            *value += dedcn[atom] * cn_gradients[atom][coordinate]
                + dedq[atom] * response.cartesian[atom * positions.len() + coordinate];
        }
        let strain = response_strain(&cn_gradients[atom], &positions, Some((lattice, periodic)));
        for (index, value) in result.virial.iter_mut().enumerate() {
            *value += dedcn[atom] * strain[index] + dedq[atom] * response.strain[9 * atom + index];
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_matches_realspace() {
        let numbers = [6, 8, 1];
        let positions = [1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.5, 3.0, 1.5];
        let lattice = [9.0, 0.0, 0.0, 0.6, 10.0, 0.0, 0.3, 0.5, 11.0];
        let param = load_param("pbe", false).unwrap();
        for mut model in [Model::D4, Model::D4S] {
            model
                .set_cutoff(Cutoff {
                    cn: 20.0,
                    disp2: 600.0,
                    ..model.cutoff
                })
                .unwrap();
            let expected = periodic_dispersion(
                &numbers, &positions, 0.25, model, param, &lattice, [true; 3],
            )
            .unwrap();
            model
                .set_ewald(Some(EwaldConfig {
                    rank: 0,
                    tolerance: 1e-12,
                    kcut: 12.0,
                    mesh: -1,
                }))
                .unwrap();
            let actual = periodic_dispersion(
                &numbers, &positions, 0.25, model, param, &lattice, [true; 3],
            )
            .unwrap();
            assert!(
                (actual.energy - expected.energy).abs() < 2e-9,
                "{} != {}",
                actual.energy,
                expected.energy
            );
            for (left, right) in actual.gradient.iter().zip(&expected.gradient) {
                assert!((left - right).abs() < 2e-9, "{left} != {right}");
            }
            for (left, right) in actual.virial.iter().zip(&expected.virial) {
                assert!((left - right).abs() < 5e-9, "{left} != {right}");
            }
        }
    }

    #[test]
    fn mesh_converges_to_direct() {
        let numbers = [6, 8, 1];
        let positions = [1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.5, 3.0, 1.5];
        let lattice = [9.0, 0.0, 0.0, 0.6, 10.0, 0.0, 0.3, 0.5, 11.0];
        let param = load_param("pbe", false).unwrap();
        let mut model = Model::D4S;
        let mut values = Vec::new();
        for mesh in [-1, 16, 32] {
            model
                .set_ewald(Some(EwaldConfig {
                    rank: 0,
                    tolerance: 1e-12,
                    kcut: 12.0,
                    mesh,
                }))
                .unwrap();
            values.push(
                periodic_dispersion(
                    &numbers, &positions, 0.25, model, param, &lattice, [true; 3],
                )
                .unwrap(),
            );
        }
        let error = |index: usize| (values[index].energy - values[0].energy).abs();
        assert!(error(2) < error(1), "{} >= {}", error(2), error(1));
        assert!(error(2) < 5e-8, "{}", error(2));
        for (left, right) in values[2].gradient.iter().zip(&values[0].gradient) {
            assert!((left - right).abs() < 5e-8, "{left} != {right}");
        }
        for (left, right) in values[2].virial.iter().zip(&values[0].virial) {
            assert!((left - right).abs() < 5e-8, "{left} != {right}");
        }
    }

    #[test]
    fn fourier_responses_and_partitions() {
        let numbers = [6, 8, 1];
        let positions = [1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.5, 3.0, 1.5];
        let lattice = [9.0, 0.0, 0.0, 0.6, 10.0, 0.0, 0.3, 0.5, 11.0];
        let param = load_param("pbe", false).unwrap();
        let fixed = [0.2, -0.4, 0.1];
        for d4s in [false, true] {
            for charge_kind in 0..3 {
                let mut model = Model::custom(d4s, 2.0, 1.0, if d4s { 6.0 } else { 4.0 })
                    .unwrap()
                    .with_ghosts(&[false, false, true]);
                model.set_charge_model(i32::from(charge_kind == 1)).unwrap();
                if charge_kind == 2 {
                    model = model.with_fixed_charges(Some(&fixed)).unwrap();
                }
                for mesh in [-1, 16] {
                    model
                        .set_ewald(Some(EwaldConfig {
                            rank: 2,
                            tolerance: 1e-8,
                            kcut: 5.0,
                            mesh,
                        }))
                        .unwrap();
                    let evaluate = |xyz: &[f64], cell: &[f64; 9], selected| {
                        periodic_dispersion(&numbers, xyz, 0.25, selected, param, cell, [true; 3])
                            .unwrap()
                    };
                    let result = evaluate(&positions, &lattice, model);
                    let step = 1e-5;
                    for coordinate in [0, 4, 7] {
                        let mut shifted = positions;
                        shifted[coordinate] += step;
                        let plus = evaluate(&shifted, &lattice, model).energy;
                        shifted[coordinate] -= 2.0 * step;
                        let minus = evaluate(&shifted, &lattice, model).energy;
                        assert!(
                            ((plus - minus) / (2.0 * step) - result.gradient[coordinate]).abs()
                                < 2e-9
                        );
                    }
                    for component in [0, 4, 7] {
                        let row = component % 3;
                        let column = component / 3;
                        let strained = |amount| {
                            let mut xyz = positions;
                            let mut cell = lattice;
                            for atom in 0..3 {
                                xyz[3 * atom + row] += amount * positions[3 * atom + column];
                                cell[3 * atom + row] += amount * lattice[3 * atom + column];
                            }
                            evaluate(&xyz, &cell, model).energy
                        };
                        assert!(
                            ((strained(step) - strained(-step)) / (2.0 * step)
                                - result.virial[component])
                                .abs()
                                < 2e-9
                        );
                    }
                    let mut energy = 0.0;
                    let mut gradient = [0.0; 9];
                    let mut virial = [0.0; 9];
                    for part in 0..3 {
                        let mut local = model;
                        local.set_work_partition(part, 3).unwrap();
                        let output = evaluate(&positions, &lattice, local);
                        energy += output.energy;
                        for (sum, value) in gradient.iter_mut().zip(output.gradient) {
                            *sum += value;
                        }
                        for (sum, value) in virial.iter_mut().zip(output.virial) {
                            *sum += value;
                        }
                    }
                    assert!((energy - result.energy).abs() < 1e-12);
                    for (left, right) in gradient.iter().zip(result.gradient) {
                        assert!((left - right).abs() < 1e-12);
                    }
                    for (left, right) in virial.iter().zip(result.virial) {
                        assert!((left - right).abs() < 1e-12);
                    }
                }
            }
        }
    }
}
