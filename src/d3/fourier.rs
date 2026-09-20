use super::*;
use crate::geometry::{determinant, inverse, lattice_repetitions, lattice_translation};
use std::f64::consts::PI;

#[derive(Clone, Copy, PartialEq)]
pub struct EwaldConfig {
    pub rank: usize,
    pub tolerance: f64,
    pub kcut: f64,
    pub mesh: i32,
}

pub struct PeriodicResult {
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub virial: [f64; 9],
}

struct LowRank {
    eigenvalues: Vec<f64>,
    vectors: Vec<f64>,
    offsets: Vec<usize>,
    dimension: usize,
}

#[derive(Clone, Copy)]
pub(crate) struct FourierTerm {
    pub(crate) prefactor: f64,
    pub(crate) numerator: i32,
    pub(crate) denominator: i32,
    pub(crate) radius: f64,
}

type SplineData = (Vec<[isize; 3]>, Vec<[[f64; 6]; 3]>, Vec<[[f64; 6]; 3]>);

#[cfg(test)]
pub fn periodic_energy(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    config: EwaldConfig,
) -> Result<f64, &'static str> {
    periodic_energy_partitioned(
        numbers,
        positions,
        lattice,
        damping,
        config,
        WorkPartition::SERIAL,
    )
}

pub fn periodic_energy_partitioned(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    config: EwaldConfig,
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    periodic_energy_partitioned_with_ghosts(
        numbers,
        positions,
        lattice,
        damping,
        config,
        &[],
        partition,
    )
}

pub fn periodic_energy_partitioned_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    config: EwaldConfig,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    periodic_energy_partitioned_with_cutoff(
        numbers,
        positions,
        lattice,
        damping,
        config,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn periodic_energy_partitioned_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    config: EwaldConfig,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let elements: Vec<_> = numbers.iter().map(|&number| number as usize - 1).collect();
    let species = unique_species(&elements);
    let lowrank = decompose(&species, config.rank, config.tolerance);
    let coordination = periodic_coordination(&elements, positions, lattice, cutoff.cn);
    let weights: Vec<_> = elements
        .iter()
        .zip(coordination)
        .map(|(&element, coordination)| weights(element, coordination))
        .collect();
    let atom_species: Vec<_> = elements
        .iter()
        .map(|element| species.iter().position(|known| known == element).unwrap())
        .collect();
    let rank = lowrank.eigenvalues.len();
    let mut factors = vec![0.0; rank * numbers.len()];
    for atom in 0..numbers.len() {
        if ghosts.get(atom) == Some(&true) {
            continue;
        }
        let species_index = atom_species[atom];
        for term in 0..rank {
            for (reference, weight) in weights[atom]
                .iter()
                .enumerate()
                .take(reference_count(elements[atom]))
            {
                factors[term * numbers.len() + atom] += weight
                    * lowrank.vectors
                        [term * lowrank.dimension + lowrank.offsets[species_index] + reference];
            }
        }
    }

    let volume = determinant(lattice).abs();
    let reciprocal = transpose(&inverse(lattice).map(|value| 2.0 * PI * value));
    let kcut = if config.kcut > 0.0 {
        config.kcut
    } else {
        automatic_kcut(&species, damping)
    };
    if config.mesh >= 0 {
        return spme_energy(
            positions,
            lattice,
            damping,
            config.mesh,
            kcut,
            &elements,
            &species,
            &atom_species,
            &lowrank,
            &factors,
            partition,
        );
    }
    let mut energy = 0.0;
    if partition.owns_index(0) {
        for atom in 0..numbers.len() {
            let terms = fourier_terms(damping, elements[atom], elements[atom])?;
            let zero: f64 = terms.iter().map(potential_zero).sum();
            for term in 0..rank {
                let factor = factors[term * numbers.len() + atom];
                energy += 0.5 * lowrank.eigenvalues[term] * factor * factor * zero;
            }
        }
    }

    let mut kernels = Vec::with_capacity(species.len() * species.len());
    for &first in &species {
        for &second in &species {
            kernels.push(fourier_terms(damping, first, second)?);
        }
    }
    let vectors = lattice_points(&reciprocal, kcut);
    energy += crate::parallel::map(vectors.len(), 256, |start, stride| {
        let mut energy = 0.0;
        let mut structure = vec![(0.0, 0.0); rank * species.len()];
        for index in (start..vectors.len()).step_by(stride) {
            let kvector = vectors[index];
            if !partition.owns_index(index) {
                continue;
            }
            let norm = norm(&kvector);
            if norm > kcut {
                continue;
            }
            if kvector
                .iter()
                .find(|&&value| value != 0.0)
                .is_some_and(|&value| value < 0.0)
            {
                continue;
            }
            let volume = volume / if norm == 0.0 { 1.0 } else { 2.0 };
            structure.fill((0.0, 0.0));
            for atom in 0..numbers.len() {
                let phase = dot(&kvector, &positions[3 * atom..3 * atom + 3]);
                let (sine, cosine) = phase.sin_cos();
                for term in 0..rank {
                    let factor = factors[term * numbers.len() + atom];
                    let item = &mut structure[term * species.len() + atom_species[atom]];
                    item.0 += factor * cosine;
                    item.1 += factor * sine;
                }
            }
            for first in 0..species.len() {
                for second in 0..=first {
                    let phi: f64 = kernels[first * species.len() + second]
                        .iter()
                        .map(|term| fourier_transform(*term, norm))
                        .sum();
                    for term in 0..rank {
                        let left = structure[term * species.len() + first];
                        let right = structure[term * species.len() + second];
                        energy -= 0.5 / volume
                            * if first == second { 1.0 } else { 2.0 }
                            * lowrank.eigenvalues[term]
                            * phi
                            * (left.0 * right.0 + left.1 * right.1);
                    }
                }
            }
        }
        energy
    })
    .into_iter()
    .sum::<f64>();
    Ok(energy)
}

pub fn periodic_derivatives_partitioned(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    config: EwaldConfig,
    partition: WorkPartition,
) -> Result<PeriodicResult, &'static str> {
    periodic_derivatives_partitioned_with_ghosts(
        numbers,
        positions,
        lattice,
        damping,
        config,
        &[],
        partition,
    )
}

pub fn periodic_derivatives_partitioned_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    config: EwaldConfig,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<PeriodicResult, &'static str> {
    periodic_derivatives_partitioned_with_cutoff(
        numbers,
        positions,
        lattice,
        damping,
        config,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn periodic_derivatives_partitioned_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    config: EwaldConfig,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<PeriodicResult, &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let elements: Vec<_> = numbers.iter().map(|&number| number as usize - 1).collect();
    let species = unique_species(&elements);
    let atom_species: Vec<_> = elements
        .iter()
        .map(|element| species.iter().position(|known| known == element).unwrap())
        .collect();
    let lowrank = decompose(&species, config.rank, config.tolerance);
    let coordination = periodic_coordination(&elements, positions, lattice, cutoff.cn);
    let weights: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .map(|(&element, &coordination)| weights(element, coordination))
        .collect();
    let derivatives: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .zip(weights.iter())
        .map(|((&element, &coordination), weights)| {
            weight_derivatives(element, coordination, weights)
        })
        .collect();
    let rank = lowrank.eigenvalues.len();
    let atoms = numbers.len();
    let mut factors = vec![0.0; rank * atoms];
    let mut factor_derivatives = vec![0.0; rank * atoms];
    for atom in 0..atoms {
        if ghosts.get(atom) == Some(&true) {
            continue;
        }
        for term in 0..rank {
            for reference in 0..reference_count(elements[atom]) {
                let vector = lowrank.vectors
                    [term * lowrank.dimension + lowrank.offsets[atom_species[atom]] + reference];
                factors[term * atoms + atom] += weights[atom][reference] * vector;
                factor_derivatives[term * atoms + atom] += derivatives[atom][reference] * vector;
            }
        }
    }
    let mut kernels = Vec::with_capacity(species.len() * species.len());
    for &first in &species {
        for &second in &species {
            kernels.push(fourier_terms(damping, first, second)?);
        }
    }
    let kcut = if config.kcut > 0.0 {
        config.kcut
    } else {
        automatic_kcut(&species, damping)
    };
    let (mut result, response) = reciprocal_derivatives(
        positions,
        lattice,
        &species,
        &atom_species,
        &kernels,
        &lowrank.eigenvalues,
        &factors,
        EwaldConfig { kcut, ..config },
        partition,
    )?;
    let dedcn: Vec<_> = (0..atoms)
        .map(|atom| {
            (0..rank)
                .map(|term| response[term * atoms + atom] * factor_derivatives[term * atoms + atom])
                .sum()
        })
        .collect();
    contract_periodic_cn(
        &elements,
        positions,
        lattice,
        cutoff.cn,
        &dedcn,
        &mut result.gradient,
        &mut result.virial,
    );
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn reciprocal_derivatives(
    positions: &[f64],
    lattice: &[f64; 9],
    species: &[usize],
    atom_species: &[usize],
    kernels: &[[FourierTerm; 2]],
    eigenvalues: &[f64],
    factors: &[f64],
    config: EwaldConfig,
    partition: WorkPartition,
) -> Result<(PeriodicResult, Vec<f64>), &'static str> {
    if config.mesh >= 0 {
        return spme_derivatives(
            positions,
            lattice,
            species,
            atom_species,
            kernels,
            eigenvalues,
            factors,
            config,
            partition,
        );
    }
    let atoms = atom_species.len();
    let rank = eigenvalues.len();
    let volume = determinant(lattice).abs();
    let reciprocal = transpose(&inverse(lattice).map(|value| 2.0 * PI * value));
    let kcut = config.kcut;
    let mut energy = 0.0;
    let mut self_energy = 0.0;
    let mut dedcn = vec![0.0; rank * atoms];
    let mut gradient = vec![0.0; 3 * atoms];
    let mut virial = [0.0; 9];
    if partition.owns_index(0) {
        for atom in 0..atoms {
            let zero: f64 = kernels[atom_species[atom] * species.len() + atom_species[atom]]
                .iter()
                .map(potential_zero)
                .sum();
            for term in 0..rank {
                let factor = factors[term * atoms + atom];
                let contribution = 0.5 * eigenvalues[term] * factor * factor * zero;
                energy += contribution;
                self_energy += contribution;
                dedcn[term * atoms + atom] += eigenvalues[term] * factor * zero;
            }
        }
    }
    let vectors = lattice_points(&reciprocal, kcut);
    let partials = crate::parallel::map(vectors.len(), 256, |start, stride| {
        let mut energy = 0.0;
        let mut dedcn = vec![0.0; rank * atoms];
        let mut gradient = vec![0.0; 3 * atoms];
        let mut virial = [0.0; 9];
        let mut phases = vec![Complex::default(); atoms];
        let mut structure = vec![Complex::default(); rank * species.len()];
        let mut phi = vec![0.0; species.len() * species.len()];
        let mut dphi = vec![0.0; phi.len()];
        let mut field = vec![Complex::default(); structure.len()];
        for index in (start..vectors.len()).step_by(stride) {
            let kvector = vectors[index];
            if !partition.owns_index(index) {
                continue;
            }
            let wave = norm(&kvector);
            if wave > kcut {
                continue;
            }
            if kvector
                .iter()
                .find(|&&value| value != 0.0)
                .is_some_and(|&value| value < 0.0)
            {
                continue;
            }
            let volume = volume / if wave == 0.0 { 1.0 } else { 2.0 };
            phases.fill(Complex::default());
            structure.fill(Complex::default());
            for atom in 0..atoms {
                let phase = dot(&kvector, &positions[3 * atom..3 * atom + 3]);
                let (sine, cosine) = phase.sin_cos();
                phases[atom] = Complex {
                    real: cosine,
                    imaginary: sine,
                };
                for term in 0..rank {
                    let item = &mut structure[term * species.len() + atom_species[atom]];
                    *item = item.add(phases[atom].scale(factors[term * atoms + atom]));
                }
            }
            phi.fill(0.0);
            dphi.fill(0.0);
            for left in 0..species.len() {
                for right in 0..=left {
                    for &fourier in &kernels[left * species.len() + right] {
                        let (value, derivative) = fourier_transform_derivative(fourier, wave);
                        phi[left * species.len() + right] += value;
                        dphi[left * species.len() + right] += derivative;
                    }
                    phi[right * species.len() + left] = phi[left * species.len() + right];
                    dphi[right * species.len() + left] = dphi[left * species.len() + right];
                }
            }
            field.fill(Complex::default());
            for left in 0..species.len() {
                for right in 0..species.len() {
                    for term in 0..rank {
                        let index = term * species.len() + left;
                        field[index] = field[index].add(
                            structure[term * species.len() + right]
                                .scale(phi[left * species.len() + right]),
                        );
                    }
                }
            }
            for atom in 0..atoms {
                for term in 0..rank {
                    let value = phases[atom]
                        .multiply(field[term * species.len() + atom_species[atom]].conjugate());
                    let scale = eigenvalues[term] / volume;
                    let contribution = -0.5 * scale * factors[term * atoms + atom] * value.real;
                    energy += contribution;
                    dedcn[term * atoms + atom] -= scale * value.real;
                    for axis in 0..3 {
                        gradient[3 * atom + axis] +=
                            scale * factors[term * atoms + atom] * value.imaginary * kvector[axis];
                    }
                }
            }
            if wave > 0.0 {
                let mut contraction = 0.0;
                for left in 0..species.len() {
                    for right in 0..species.len() {
                        for term in 0..rank {
                            contraction += eigenvalues[term]
                                * dphi[left * species.len() + right]
                                * structure[term * species.len() + left]
                                    .multiply(structure[term * species.len() + right].conjugate())
                                    .real;
                        }
                    }
                }
                for row in 0..3 {
                    for column in 0..3 {
                        virial[row + 3 * column] +=
                            0.5 * contraction * kvector[row] * kvector[column] / (wave * volume);
                    }
                }
            }
        }
        (energy, dedcn, gradient, virial)
    });
    for (partial_energy, partial_cn, partial_gradient, partial_virial) in partials {
        energy += partial_energy;
        for (total, value) in dedcn.iter_mut().zip(partial_cn) {
            *total += value;
        }
        for (total, value) in gradient.iter_mut().zip(partial_gradient) {
            *total += value;
        }
        for (total, value) in virial.iter_mut().zip(partial_virial) {
            *total += value;
        }
    }
    let reciprocal_energy = energy - self_energy;
    for axis in 0..3 {
        virial[axis + 3 * axis] -= reciprocal_energy;
    }
    Ok((
        PeriodicResult {
            energy,
            gradient,
            virial,
        },
        dedcn,
    ))
}

fn unique_species(elements: &[usize]) -> Vec<usize> {
    let mut species = Vec::new();
    for &element in elements {
        if !species.contains(&element) {
            species.push(element);
        }
    }
    species
}

fn decompose(species: &[usize], requested_rank: usize, tolerance: f64) -> LowRank {
    let mut offsets = Vec::with_capacity(species.len());
    let mut dimension = 0;
    for &element in species {
        offsets.push(dimension);
        dimension += reference_count(element);
    }
    let mut matrix = vec![0.0; dimension * dimension];
    for (first_species, &first_element) in species.iter().enumerate() {
        for first_reference in 0..reference_count(first_element) {
            let first = offsets[first_species] + first_reference;
            for (second_species, &second_element) in species.iter().enumerate() {
                for second_reference in 0..reference_count(second_element) {
                    let second = offsets[second_species] + second_reference;
                    matrix[first * dimension + second] = reference_c6(
                        first_element,
                        second_element,
                        first_reference,
                        second_reference,
                    );
                }
            }
        }
    }
    let (eigenvalues, vectors) = factorize(matrix, dimension, requested_rank, tolerance);
    LowRank {
        eigenvalues,
        vectors,
        offsets,
        dimension,
    }
}

pub(crate) fn factorize(
    matrix: Vec<f64>,
    dimension: usize,
    requested_rank: usize,
    tolerance: f64,
) -> (Vec<f64>, Vec<f64>) {
    let original = matrix.clone();
    let (eigenvalues, eigenvectors) = jacobi(matrix, dimension);
    let limit = requested_rank.min(dimension);
    let mut rank = if requested_rank > 0 { limit } else { dimension };
    if requested_rank == 0 {
        let floor = original.iter().map(|value| value.abs()).fold(0.0, f64::max) * 1.0e-6;
        let mut residual = original.clone();
        for term in 0..dimension {
            for row in 0..dimension {
                for column in 0..dimension {
                    residual[row * dimension + column] -= eigenvalues[term]
                        * eigenvectors[term * dimension + row]
                        * eigenvectors[term * dimension + column];
                }
            }
            let error = residual
                .iter()
                .zip(original.iter())
                .map(|(left, right)| left.abs() / right.abs().max(floor))
                .fold(0.0, f64::max);
            if error <= tolerance {
                rank = term + 1;
                break;
            }
        }
    }
    (
        eigenvalues[..rank].to_vec(),
        eigenvectors[..rank * dimension].to_vec(),
    )
}

fn jacobi(mut matrix: Vec<f64>, size: usize) -> (Vec<f64>, Vec<f64>) {
    let mut vectors = vec![0.0; size * size];
    for index in 0..size {
        vectors[index * size + index] = 1.0;
    }
    let norm = matrix.iter().map(|value| value * value).sum::<f64>().sqrt();
    let skip = 0.01 * f64::EPSILON * norm;
    for _ in 0..100 {
        let off = (0..size)
            .flat_map(|row| (row + 1..size).map(move |column| (row, column)))
            .map(|(row, column)| matrix[row * size + column].powi(2))
            .sum::<f64>();
        if (2.0 * off).sqrt() <= f64::EPSILON * norm {
            break;
        }
        for first in 0..size - 1 {
            for second in first + 1..size {
                let offdiag = matrix[first * size + second];
                if offdiag.abs() <= skip {
                    continue;
                }
                let theta =
                    0.5 * (matrix[second * size + second] - matrix[first * size + first]) / offdiag;
                let tangent = if theta >= 0.0 { 1.0 } else { -1.0 }
                    / (theta.abs() + (1.0 + theta * theta).sqrt());
                let cosine = 1.0 / (1.0 + tangent * tangent).sqrt();
                let sine = tangent * cosine;
                for row in 0..size {
                    let left = matrix[row * size + first];
                    let right = matrix[row * size + second];
                    matrix[row * size + first] = cosine * left - sine * right;
                    matrix[row * size + second] = sine * left + cosine * right;
                }
                for column in 0..size {
                    let left = matrix[first * size + column];
                    let right = matrix[second * size + column];
                    matrix[first * size + column] = cosine * left - sine * right;
                    matrix[second * size + column] = sine * left + cosine * right;
                }
                matrix[first * size + second] = 0.0;
                matrix[second * size + first] = 0.0;
                for row in 0..size {
                    let left = vectors[first * size + row];
                    let right = vectors[second * size + row];
                    vectors[first * size + row] = cosine * left - sine * right;
                    vectors[second * size + row] = sine * left + cosine * right;
                }
            }
        }
    }
    let mut order: Vec<_> = (0..size).collect();
    order.sort_by(|&left, &right| {
        matrix[right * size + right]
            .abs()
            .total_cmp(&matrix[left * size + left].abs())
    });
    let eigenvalues = order
        .iter()
        .map(|&index| matrix[index * size + index])
        .collect();
    let mut eigenvectors = Vec::with_capacity(size * size);
    for column in order {
        eigenvectors.extend((0..size).map(|row| vectors[column * size + row]));
    }
    (eigenvalues, eigenvectors)
}

fn periodic_coordination(
    elements: &[usize],
    positions: &[f64],
    lattice: &[f64; 9],
    cutoff: f64,
) -> Vec<f64> {
    let translations = lattice_points(lattice, cutoff);
    let mut coordination = vec![0.0; elements.len()];
    let partials = crate::parallel::map(elements.len(), 64, |start, stride| {
        let mut coordination = vec![0.0; elements.len()];
        for first in (start..elements.len()).step_by(stride) {
            for second in 0..=first {
                let radius = value(COVALENT_RADII, elements[first])
                    + value(COVALENT_RADII, elements[second]);
                for translation in &translations {
                    let vector: [f64; 3] = std::array::from_fn(|axis| {
                        positions[3 * first + axis]
                            - positions[3 * second + axis]
                            - translation[axis]
                    });
                    let distance2 = dot(&vector, &vector);
                    if distance2 > cutoff * cutoff || distance2 < 1.0e-12 {
                        continue;
                    }
                    let count = 1.0 / (1.0 + (-16.0 * (radius / distance2.sqrt() - 1.0)).exp());
                    coordination[first] += count;
                    if first != second {
                        coordination[second] += count;
                    }
                }
            }
        }
        coordination
    });
    for partial in partials {
        for (total, value) in coordination.iter_mut().zip(partial) {
            *total += value;
        }
    }
    coordination
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "the match guard restricts alpha to 4..=64"
)]
fn fourier_terms(
    damping: Damping,
    first: usize,
    second: usize,
) -> Result<[FourierTerm; 2], &'static str> {
    let rrij = 3.0 * value(R4R2, first) * value(R4R2, second);
    match damping {
        Damping::Rational { s6, s8, a1, a2 } => {
            let radius = a1 * rrij.sqrt() + a2;
            Ok([
                FourierTerm {
                    prefactor: s6,
                    numerator: 0,
                    denominator: 6,
                    radius,
                },
                FourierTerm {
                    prefactor: s8 * rrij,
                    numerator: 0,
                    denominator: 8,
                    radius,
                },
            ])
        }
        Damping::Zero {
            s6,
            s8,
            rs6,
            rs8,
            alpha,
        } if (alpha - alpha.round()).abs() <= f64::EPSILON.sqrt()
            && (4.0..=64.0).contains(&alpha.round()) =>
        {
            let vdw = value(VDW_RADII, pair_index(first, second));
            // Bounded above so the pole sum stays short and alpha+2 fits i32,
            // and below so the numerator never lands on the sin(0) pole.
            let alpha6 = alpha.round() as i32;
            let alpha8 = alpha6 + 2;
            Ok([
                FourierTerm {
                    prefactor: s6,
                    numerator: alpha6 - 6,
                    denominator: alpha6,
                    radius: 6.0_f64.powf(1.0 / alpha) * rs6 * vdw,
                },
                FourierTerm {
                    prefactor: s8 * rrij,
                    numerator: alpha8 - 8,
                    denominator: alpha8,
                    radius: 6.0_f64.powf(1.0 / (alpha + 2.0)) * rs8 * vdw,
                },
            ])
        }
        _ => Err("damping function has no analytical Fourier transform"),
    }
}

fn fourier_transform(term: FourierTerm, wave: f64) -> f64 {
    fourier_transform_derivative(term, wave).0
}

fn fourier_transform_derivative(term: FourierTerm, wave: f64) -> (f64, f64) {
    if term.prefactor == 0.0 {
        return (0.0, 0.0);
    }
    let exponent = term.numerator - term.denominator + 2;
    let prefactor =
        4.0 * PI * PI / term.denominator as f64 * term.prefactor * term.radius.powi(exponent);
    if wave <= 0.0 {
        return (
            prefactor * term.radius
                / (PI * (term.numerator + 3) as f64 / term.denominator as f64).sin(),
            0.0,
        );
    }
    let scaled = wave * term.radius;
    let (sum, derivative) = (0..term.denominator / 2)
        .map(|pole| {
            let angle = PI * (2 * pole + 1) as f64 / term.denominator as f64;
            let argument = 0.5 * PI + exponent as f64 * angle + scaled * angle.cos();
            let damping = (-scaled * angle.sin()).exp();
            (argument.sin() * damping, (argument + angle).cos() * damping)
        })
        .fold((0.0, 0.0), |sum, value| (sum.0 + value.0, sum.1 + value.1));
    let value = prefactor * sum / wave;
    (value, (prefactor * term.radius * derivative - value) / wave)
}

fn potential_zero(term: &FourierTerm) -> f64 {
    if term.numerator > 0 {
        0.0
    } else {
        term.prefactor / term.radius.powi(term.denominator)
    }
}

fn automatic_kcut(species: &[usize], damping: Damping) -> f64 {
    let mut cutoff: f64 = 0.0;
    for &first in species {
        for &second in species {
            if let Ok(terms) = fourier_terms(damping, first, second) {
                for term in terms {
                    cutoff = cutoff.max(
                        -1.0e-8_f64.ln() / (term.radius * (PI / term.denominator as f64).sin()),
                    );
                }
            }
        }
    }
    cutoff
}

fn lattice_points(lattice: &[f64; 9], cutoff: f64) -> Vec<[f64; 3]> {
    let repetitions = lattice_repetitions(lattice, [true; 3], cutoff);
    let mut result = Vec::new();
    for first in -repetitions[0]..=repetitions[0] {
        for second in -repetitions[1]..=repetitions[1] {
            for third in -repetitions[2]..=repetitions[2] {
                result.push(lattice_translation(lattice, [first, second, third]));
            }
        }
    }
    result
}

fn transpose(matrix: &[f64; 9]) -> [f64; 9] {
    [
        matrix[0], matrix[3], matrix[6], matrix[1], matrix[4], matrix[7], matrix[2], matrix[5],
        matrix[8],
    ]
}

fn dot(first: &[f64], second: &[f64]) -> f64 {
    first
        .iter()
        .zip(second)
        .map(|(left, right)| left * right)
        .sum()
}

fn norm(vector: &[f64]) -> f64 {
    dot(vector, vector).sqrt()
}

#[derive(Clone, Copy, Default)]
struct Complex {
    real: f64,
    imaginary: f64,
}

impl Complex {
    fn add(self, rhs: Self) -> Self {
        Self {
            real: self.real + rhs.real,
            imaginary: self.imaginary + rhs.imaginary,
        }
    }
    fn sub(self, rhs: Self) -> Self {
        Self {
            real: self.real - rhs.real,
            imaginary: self.imaginary - rhs.imaginary,
        }
    }
    fn multiply(self, rhs: Self) -> Self {
        Self {
            real: self.real * rhs.real - self.imaginary * rhs.imaginary,
            imaginary: self.real * rhs.imaginary + self.imaginary * rhs.real,
        }
    }
    fn scale(self, factor: f64) -> Self {
        Self {
            real: self.real * factor,
            imaginary: self.imaginary * factor,
        }
    }
    fn conjugate(self) -> Self {
        Self {
            real: self.real,
            imaginary: -self.imaginary,
        }
    }
}

fn contract_periodic_cn(
    elements: &[usize],
    positions: &[f64],
    lattice: &[f64; 9],
    cutoff: f64,
    dedcn: &[f64],
    gradient: &mut [f64],
    virial: &mut [f64; 9],
) {
    let translations = lattice_points(lattice, cutoff);
    let partials = crate::parallel::map(elements.len(), 64, |start, stride| {
        let mut gradient = vec![0.0; positions.len()];
        let mut virial = [0.0; 9];
        for first in (start..elements.len()).step_by(stride) {
            for second in 0..=first {
                let radius = value(COVALENT_RADII, elements[first])
                    + value(COVALENT_RADII, elements[second]);
                for translation in &translations {
                    let vector: [f64; 3] = std::array::from_fn(|axis| {
                        positions[3 * first + axis]
                            - positions[3 * second + axis]
                            - translation[axis]
                    });
                    let distance2 = dot(&vector, &vector);
                    if distance2 > cutoff * cutoff || distance2 < 1.0e-12 {
                        continue;
                    }
                    let distance = distance2.sqrt();
                    let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
                    let derivative = -count * (1.0 - count) * 16.0 * radius / distance2;
                    let scale = dedcn[first] + dedcn[second];
                    for axis in 0..3 {
                        let component = derivative * vector[axis] / distance * scale;
                        gradient[3 * first + axis] += component;
                        gradient[3 * second + axis] -= component;
                        for other in 0..3 {
                            let virial_scale = if first == second { dedcn[first] } else { scale };
                            virial[axis + 3 * other] +=
                                derivative * vector[axis] / distance * vector[other] * virial_scale;
                        }
                    }
                }
            }
        }
        (gradient, virial)
    });
    for (partial_gradient, partial_virial) in partials {
        for (total, value) in gradient.iter_mut().zip(partial_gradient) {
            *total += value;
        }
        for (total, value) in virial.iter_mut().zip(partial_virial) {
            *total += value;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn spme_derivatives(
    positions: &[f64],
    lattice: &[f64; 9],
    species: &[usize],
    atom_species: &[usize],
    terms: &[[FourierTerm; 2]],
    eigenvalues: &[f64],
    factors: &[f64],
    config: EwaldConfig,
    partition: WorkPartition,
) -> Result<(PeriodicResult, Vec<f64>), &'static str> {
    let atoms = atom_species.len();
    let rank = eigenvalues.len();
    let kcut = config.kcut;
    let (mesh, count) = spme_mesh(lattice, config.mesh, kcut)?;
    let inverse_lattice = inverse(lattice);
    let reciprocal = transpose(&inverse_lattice.map(|value| 2.0 * PI * value));
    let volume = determinant(lattice).abs();
    let (bases, splines, spline_derivatives) = spline_weights(positions, &inverse_lattice, mesh);
    let euler = euler_factors(mesh);
    let grid_data: Vec<_> = (0..count)
        .map(|index| {
            let grid = grid_coordinates(index, mesh);
            let wave_index: [f64; 3] = std::array::from_fn(|axis| {
                (grid[axis] as isize
                    - mesh[axis] as isize * ((2 * grid[axis]) / mesh[axis]) as isize)
                    as f64
            });
            let wave: [f64; 3] = std::array::from_fn(|axis| {
                (0..3)
                    .map(|column| reciprocal[axis + 3 * column] * wave_index[column])
                    .sum()
            });
            (
                grid,
                euler[0][grid[0]] * euler[1][grid[1]] * euler[2][grid[2]],
                wave,
                norm(&wave),
            )
        })
        .collect();
    let mut kernels = vec![vec![0.0; count]; species.len() * species.len()];
    let mut kernel_derivatives = kernels.clone();
    for left in 0..species.len() {
        for right in 0..=left {
            let pair_terms = terms[left * species.len() + right];
            for index in 0..count {
                let wave_norm = grid_data[index].3;
                for &fourier in &pair_terms {
                    let (value, derivative) = fourier_transform_derivative(fourier, wave_norm);
                    kernels[left * species.len() + right][index] += value;
                    kernel_derivatives[left * species.len() + right][index] += derivative;
                }
                kernels[right * species.len() + left][index] =
                    kernels[left * species.len() + right][index];
                kernel_derivatives[right * species.len() + left][index] =
                    kernel_derivatives[left * species.len() + right][index];
            }
        }
    }
    let mut energy = 0.0;
    let mut reciprocal_energy = 0.0;
    let mut dedcn = vec![0.0; rank * atoms];
    let mut gradient = vec![0.0; 3 * atoms];
    let mut virial = [0.0; 9];
    for atom in 0..atoms {
        let zero: f64 = terms[atom_species[atom] * species.len() + atom_species[atom]]
            .iter()
            .map(potential_zero)
            .sum();
        for term in 0..rank {
            if !partition.owns_index(term) {
                continue;
            }
            let factor = factors[term * atoms + atom];
            energy += 0.5 * eigenvalues[term] * factor * factor * zero;
            dedcn[term * atoms + atom] += eigenvalues[term] * factor * zero;
        }
    }
    let derivative_transform: [[f64; 3]; 3] = std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            mesh[column] as f64 * reciprocal[row + 3 * column] / (2.0 * PI)
        })
    });
    let forward_twiddles = mesh.map(|size| fft_twiddles(size, 1));
    let reverse_twiddles = mesh.map(|size| fft_twiddles(size, -1));
    let partials = crate::parallel::map(rank, 2, |start, stride| {
        let mut energy = 0.0;
        let mut dedcn = vec![0.0; rank * atoms];
        let mut gradient = vec![0.0; 3 * atoms];
        let mut virial = [0.0; 9];
        let mut fft_line = vec![Complex::default(); *mesh.iter().max().unwrap()];
        for term in (start..rank).step_by(stride) {
            if !partition.owns_index(term) {
                continue;
            }
            let mut grids = vec![vec![Complex::default(); count]; species.len()];
            for atom in 0..atoms {
                for zpoint in 0..6 {
                    let z = modulo(bases[atom][2] + zpoint as isize + 1, mesh[2]);
                    for ypoint in 0..6 {
                        let y = modulo(bases[atom][1] + ypoint as isize + 1, mesh[1]);
                        for xpoint in 0..6 {
                            let x = modulo(bases[atom][0] + xpoint as isize + 1, mesh[0]);
                            let weight = splines[atom][0][xpoint]
                                * splines[atom][1][ypoint]
                                * splines[atom][2][zpoint];
                            grids[atom_species[atom]][grid_index(x, y, z, mesh)].real +=
                                factors[term * atoms + atom] * weight;
                        }
                    }
                }
            }
            for grid in &mut grids {
                fft_3d(grid, mesh, &forward_twiddles, &mut fft_line);
            }
            for index in 0..count {
                let (_, correction, wave, wave_norm) = grid_data[index];
                if wave_norm > 0.0 {
                    let mut contraction = 0.0;
                    for left in 0..species.len() {
                        for right in 0..species.len() {
                            contraction += eigenvalues[term]
                                * correction
                                * kernel_derivatives[left * species.len() + right][index]
                                * grids[left][index]
                                    .multiply(grids[right][index].conjugate())
                                    .real;
                        }
                    }
                    for row in 0..3 {
                        for column in 0..3 {
                            virial[row + 3 * column] +=
                                0.5 * contraction * wave[row] * wave[column] / (wave_norm * volume);
                        }
                    }
                }
            }
            for left in 0..species.len() {
                let mut potential = vec![Complex::default(); count];
                for index in 0..count {
                    let correction = grid_data[index].1;
                    for right in 0..species.len() {
                        potential[index] = potential[index].add(
                            grids[right][index]
                                .scale(correction * kernels[left * species.len() + right][index]),
                        );
                    }
                }
                fft_3d(&mut potential, mesh, &reverse_twiddles, &mut fft_line);
                for atom in 0..atoms {
                    if atom_species[atom] != left {
                        continue;
                    }
                    let mut gathered = 0.0;
                    let mut gathered_derivative = [0.0; 3];
                    for zpoint in 0..6 {
                        let z = modulo(bases[atom][2] + zpoint as isize + 1, mesh[2]);
                        for ypoint in 0..6 {
                            let y = modulo(bases[atom][1] + ypoint as isize + 1, mesh[1]);
                            for xpoint in 0..6 {
                                let x = modulo(bases[atom][0] + xpoint as isize + 1, mesh[0]);
                                let potential = potential[grid_index(x, y, z, mesh)].real;
                                let spline = [
                                    splines[atom][0][xpoint],
                                    splines[atom][1][ypoint],
                                    splines[atom][2][zpoint],
                                ];
                                gathered += spline[0] * spline[1] * spline[2] * potential;
                                let mesh_derivative = [
                                    spline_derivatives[atom][0][xpoint] * spline[1] * spline[2],
                                    spline[0] * spline_derivatives[atom][1][ypoint] * spline[2],
                                    spline[0] * spline[1] * spline_derivatives[atom][2][zpoint],
                                ];
                                for axis in 0..3 {
                                    gathered_derivative[axis] += (0..3)
                                        .map(|direction| {
                                            derivative_transform[axis][direction]
                                                * mesh_derivative[direction]
                                        })
                                        .sum::<f64>()
                                        * potential;
                                }
                            }
                        }
                    }
                    let scale = -0.5 / volume * eigenvalues[term];
                    let contribution = scale * factors[term * atoms + atom] * gathered;
                    energy += contribution;
                    dedcn[term * atoms + atom] += 2.0 * scale * gathered;
                    for axis in 0..3 {
                        gradient[3 * atom + axis] +=
                            2.0 * scale * factors[term * atoms + atom] * gathered_derivative[axis];
                    }
                }
            }
        }
        (energy, dedcn, gradient, virial)
    });
    for (partial_energy, partial_cn, partial_gradient, partial_virial) in partials {
        energy += partial_energy;
        reciprocal_energy += partial_energy;
        for (total, value) in dedcn.iter_mut().zip(partial_cn) {
            *total += value;
        }
        for (total, value) in gradient.iter_mut().zip(partial_gradient) {
            *total += value;
        }
        for (total, value) in virial.iter_mut().zip(partial_virial) {
            *total += value;
        }
    }
    for axis in 0..3 {
        virial[axis + 3 * axis] -= reciprocal_energy;
    }
    Ok((
        PeriodicResult {
            energy,
            gradient,
            virial,
        },
        dedcn,
    ))
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "a saturated size is rejected by the following checked_next_power_of_two"
)]
fn spme_mesh(
    lattice: &[f64; 9],
    requested_mesh: i32,
    kcut: f64,
) -> Result<([usize; 3], usize), &'static str> {
    let mut mesh = [0; 3];
    for (direction, size) in mesh.iter_mut().enumerate() {
        let requested = if requested_mesh > 0 {
            requested_mesh as usize
        } else {
            let length = norm(&[
                lattice[3 * direction],
                lattice[3 * direction + 1],
                lattice[3 * direction + 2],
            ]);
            let automatic = (kcut * length / PI).ceil();
            if !automatic.is_finite() {
                return Err("Ewald mesh size is not representable");
            }
            (automatic as usize).max(1)
        };
        *size = requested
            .checked_next_power_of_two()
            .ok_or("Ewald mesh size is not representable")?;
    }
    let count = mesh
        .iter()
        .try_fold(1usize, |count, &size| count.checked_mul(size))
        .ok_or("Ewald mesh size is not representable")?;
    std::alloc::Layout::array::<([usize; 3], f64, [f64; 3], f64)>(count)
        .map_err(|_| "Ewald mesh storage is not representable")?;
    Ok((mesh, count))
}

#[allow(clippy::too_many_arguments)]
fn spme_energy(
    positions: &[f64],
    lattice: &[f64; 9],
    damping: Damping,
    requested_mesh: i32,
    kcut: f64,
    elements: &[usize],
    species: &[usize],
    atom_species: &[usize],
    lowrank: &LowRank,
    factors: &[f64],
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    let (mesh, count) = spme_mesh(lattice, requested_mesh, kcut)?;
    let inverse_lattice = inverse(lattice);
    let reciprocal = transpose(&inverse_lattice.map(|value| 2.0 * PI * value));
    let volume = determinant(lattice).abs();
    let (bases, splines, _) = spline_weights(positions, &inverse_lattice, mesh);
    let euler = euler_factors(mesh);
    let mut kernels = vec![vec![0.0; count]; species.len() * species.len()];
    let mut wave_norms = vec![0.0; count];
    for third in 0..mesh[2] {
        for second in 0..mesh[1] {
            for first in 0..mesh[0] {
                let grid = [first, second, third];
                let wave_index: [f64; 3] = std::array::from_fn(|axis| {
                    let index = grid[axis];
                    (index as isize - mesh[axis] as isize * ((2 * index) / mesh[axis]) as isize)
                        as f64
                });
                let wave: [f64; 3] = std::array::from_fn(|axis| {
                    (0..3)
                        .map(|column| reciprocal[axis + 3 * column] * wave_index[column])
                        .sum()
                });
                let wave_norm = norm(&wave);
                let index = grid_index(first, second, third, mesh);
                wave_norms[index] = wave_norm;
            }
        }
    }
    for left in 0..species.len() {
        for right in 0..=left {
            let terms = fourier_terms(damping, species[left], species[right])?;
            for (index, &wave_norm) in wave_norms.iter().enumerate() {
                let value = terms
                    .iter()
                    .map(|term| fourier_transform(*term, wave_norm))
                    .sum();
                kernels[left * species.len() + right][index] = value;
                kernels[right * species.len() + left][index] = value;
            }
        }
    }

    let rank = lowrank.eigenvalues.len();
    let atoms = elements.len();
    let mut energy = 0.0;
    let forward_twiddles = mesh.map(|size| fft_twiddles(size, 1));
    let reverse_twiddles = mesh.map(|size| fft_twiddles(size, -1));
    for atom in 0..atoms {
        let zero: f64 = fourier_terms(damping, elements[atom], elements[atom])?
            .iter()
            .map(potential_zero)
            .sum();
        for term in 0..rank {
            if !partition.owns_index(term) {
                continue;
            }
            let factor = factors[term * atoms + atom];
            energy += 0.5 * lowrank.eigenvalues[term] * factor * factor * zero;
        }
    }

    energy += crate::parallel::map(rank, 2, |start, stride| {
        let mut energy = 0.0;
        let mut fft_line = vec![Complex::default(); *mesh.iter().max().unwrap()];
        for term in (start..rank).step_by(stride) {
            if !partition.owns_index(term) {
                continue;
            }
            let mut grids = vec![vec![Complex::default(); count]; species.len()];
            for atom in 0..atoms {
                for third in 0..6 {
                    let z = modulo(bases[atom][2] + third as isize + 1, mesh[2]);
                    for second in 0..6 {
                        let y = modulo(bases[atom][1] + second as isize + 1, mesh[1]);
                        for first in 0..6 {
                            let x = modulo(bases[atom][0] + first as isize + 1, mesh[0]);
                            let weight = splines[atom][0][first]
                                * splines[atom][1][second]
                                * splines[atom][2][third];
                            grids[atom_species[atom]][grid_index(x, y, z, mesh)].real +=
                                factors[term * atoms + atom] * weight;
                        }
                    }
                }
            }
            for grid in &mut grids {
                fft_3d(grid, mesh, &forward_twiddles, &mut fft_line);
            }
            for left in 0..species.len() {
                let mut potential = vec![Complex::default(); count];
                for index in 0..count {
                    let [x, y, z] = grid_coordinates(index, mesh);
                    let correction = euler[0][x] * euler[1][y] * euler[2][z];
                    for right in 0..species.len() {
                        potential[index] = potential[index].add(
                            grids[right][index]
                                .scale(correction * kernels[left * species.len() + right][index]),
                        );
                    }
                }
                fft_3d(&mut potential, mesh, &reverse_twiddles, &mut fft_line);
                for atom in 0..atoms {
                    if atom_species[atom] != left {
                        continue;
                    }
                    let mut gathered = 0.0;
                    for third in 0..6 {
                        let z = modulo(bases[atom][2] + third as isize + 1, mesh[2]);
                        for second in 0..6 {
                            let y = modulo(bases[atom][1] + second as isize + 1, mesh[1]);
                            for first in 0..6 {
                                let x = modulo(bases[atom][0] + first as isize + 1, mesh[0]);
                                gathered += splines[atom][0][first]
                                    * splines[atom][1][second]
                                    * splines[atom][2][third]
                                    * potential[grid_index(x, y, z, mesh)].real;
                            }
                        }
                    }
                    energy -= 0.5 / volume
                        * lowrank.eigenvalues[term]
                        * factors[term * atoms + atom]
                        * gathered;
                }
            }
        }
        energy
    })
    .into_iter()
    .sum::<f64>();
    Ok(energy)
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "the fractional offset keeps the value within the validated mesh size"
)]
fn spline_weights(positions: &[f64], inverse: &[f64; 9], mesh: [usize; 3]) -> SplineData {
    let atoms = positions.len() / 3;
    let mut bases = vec![[0; 3]; atoms];
    let mut weights = vec![[[0.0; 6]; 3]; atoms];
    let mut derivatives = vec![[[0.0; 6]; 3]; atoms];
    for atom in 0..atoms {
        let fraction: [f64; 3] = std::array::from_fn(|row| {
            (0..3)
                .map(|column| inverse[row + 3 * column] * positions[3 * atom + column])
                .sum()
        });
        for direction in 0..3 {
            let coordinate =
                mesh[direction] as f64 * (fraction[direction] - fraction[direction].floor());
            bases[atom][direction] = coordinate.floor() as isize - 6;
            (weights[atom][direction], derivatives[atom][direction]) =
                bspline_derivative(coordinate - coordinate.floor());
        }
    }
    (bases, weights, derivatives)
}

fn bspline(offset: f64) -> [f64; 6] {
    bspline_derivative(offset).0
}

fn bspline_derivative(offset: f64) -> ([f64; 6], [f64; 6]) {
    let mut theta = [0.0; 6];
    let mut derivative = [0.0; 6];
    theta[0] = 1.0 - offset;
    theta[1] = offset;
    for order in 3..6 {
        let divisor = 1.0 / (order - 1) as f64;
        theta[order - 1] = divisor * offset * theta[order - 2];
        for index in 1..=order - 2 {
            let target = order - index - 1;
            theta[target] = divisor
                * ((offset + index as f64) * theta[target - 1]
                    + (order as f64 - index as f64 - offset) * theta[target]);
        }
        theta[0] *= divisor * (1.0 - offset);
    }
    derivative[0] = -theta[0];
    for index in 1..6 {
        derivative[index] = theta[index - 1] - theta[index];
    }
    let divisor = 0.2;
    theta[5] = divisor * offset * theta[4];
    for index in 1..=4 {
        let target = 5 - index;
        theta[target] = divisor
            * ((offset + index as f64) * theta[target - 1]
                + (6.0 - index as f64 - offset) * theta[target]);
    }
    theta[0] *= divisor * (1.0 - offset);
    (theta, derivative)
}

fn euler_factors(mesh: [usize; 3]) -> [Vec<f64>; 3] {
    let spline = bspline(0.0);
    std::array::from_fn(|direction| {
        (0..mesh[direction])
            .map(|wave| {
                let mut denominator = Complex::default();
                for (point, spline) in spline.iter().enumerate().take(5) {
                    let angle = 2.0 * PI * wave as f64 * point as f64 / mesh[direction] as f64;
                    denominator = denominator.add(
                        Complex {
                            real: angle.cos(),
                            imaginary: angle.sin(),
                        }
                        .scale(*spline),
                    );
                }
                1.0 / (denominator.real * denominator.real
                    + denominator.imaginary * denominator.imaginary)
            })
            .collect()
    })
}

fn fft_3d(
    grid: &mut [Complex],
    mesh: [usize; 3],
    twiddles: &[Vec<Complex>; 3],
    line: &mut [Complex],
) {
    for z in 0..mesh[2] {
        for y in 0..mesh[1] {
            for x in 0..mesh[0] {
                line[x] = grid[grid_index(x, y, z, mesh)];
            }
            fft(&mut line[..mesh[0]], &twiddles[0]);
            for x in 0..mesh[0] {
                grid[grid_index(x, y, z, mesh)] = line[x];
            }
        }
    }
    for z in 0..mesh[2] {
        for x in 0..mesh[0] {
            for y in 0..mesh[1] {
                line[y] = grid[grid_index(x, y, z, mesh)];
            }
            fft(&mut line[..mesh[1]], &twiddles[1]);
            for y in 0..mesh[1] {
                grid[grid_index(x, y, z, mesh)] = line[y];
            }
        }
    }
    for y in 0..mesh[1] {
        for x in 0..mesh[0] {
            for z in 0..mesh[2] {
                line[z] = grid[grid_index(x, y, z, mesh)];
            }
            fft(&mut line[..mesh[2]], &twiddles[2]);
            for z in 0..mesh[2] {
                grid[grid_index(x, y, z, mesh)] = line[z];
            }
        }
    }
}

fn fft_twiddles(size: usize, sign: i32) -> Vec<Complex> {
    let mut twiddles = Vec::with_capacity(size - 1);
    let mut half = 1;
    while half < size {
        let step = 2 * half;
        twiddles.extend((0..half).map(|offset| {
            let angle = sign as f64 * 2.0 * PI * offset as f64 / step as f64;
            Complex {
                real: angle.cos(),
                imaginary: angle.sin(),
            }
        }));
        half = step;
    }
    twiddles
}

fn fft(values: &mut [Complex], twiddles: &[Complex]) {
    let size = values.len();
    let mut target = 0;
    for source in 0..size {
        if target > source {
            values.swap(target, source);
        }
        let mut bit = size / 2;
        while bit >= 2 && target >= bit {
            target -= bit;
            bit /= 2;
        }
        target += bit;
    }
    let mut half = 1;
    let mut twiddle_start = 0;
    while half < size {
        let step = 2 * half;
        for offset in 0..half {
            let twiddle = twiddles[twiddle_start + offset];
            for index in (offset..size).step_by(step) {
                let value = twiddle.multiply(values[index + half]);
                values[index + half] = values[index].sub(value);
                values[index] = values[index].add(value);
            }
        }
        twiddle_start += half;
        half = step;
    }
}

fn grid_index(x: usize, y: usize, z: usize, mesh: [usize; 3]) -> usize {
    x + mesh[0] * (y + mesh[1] * z)
}

fn grid_coordinates(index: usize, mesh: [usize; 3]) -> [usize; 3] {
    [
        index % mesh[0],
        (index / mesh[0]) % mesh[1],
        index / (mesh[0] * mesh[1]),
    ]
}

fn modulo(value: isize, modulus: usize) -> usize {
    value.rem_euclid(modulus as isize) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_ewald_matches_cubic_reference() {
        let energy = periodic_energy(
            &[6, 8],
            &[0.4, 0.8, 1.2, 4.4, 2.8, 3.6],
            &[8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0],
            Damping::Rational {
                s6: 1.0,
                s8: 0.7875,
                a1: 0.4289,
                a2: 4.4407,
            },
            EwaldConfig {
                rank: 10_000,
                tolerance: 1.0e-4,
                kcut: 10.0,
                mesh: -1,
            },
        )
        .unwrap();
        assert!(
            (energy - -0.0028149228007689916).abs() < 1.0e-12,
            "{energy}"
        );
        let result = periodic_derivatives_partitioned(
            &[6, 8],
            &[0.4, 0.8, 1.2, 4.4, 2.8, 3.6],
            &[8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0],
            Damping::Rational {
                s6: 1.0,
                s8: 0.7875,
                a1: 0.4289,
                a2: 4.4407,
            },
            EwaldConfig {
                rank: 10_000,
                tolerance: 1.0e-4,
                kcut: 10.0,
                mesh: -1,
            },
            WorkPartition::SERIAL,
        )
        .unwrap();
        assert!((result.energy - energy).abs() < 1.0e-14);
        assert!((result.gradient[0] - -5.428996265343892e-20).abs() < 1.0e-13);
        assert!((result.gradient[1] - 2.939667997396428e-5).abs() < 1.0e-13);
        assert!((result.gradient[2] - 5.2457766011066856e-5).abs() < 1.0e-13);
        assert!((result.virial[0] - 0.003746104731458539).abs() < 1.0e-12);
        assert!((result.virial[4] - 0.003644616810547004).abs() < 1.0e-12);
        assert!((result.virial[5] - -0.0002457998363346487).abs() < 1.0e-12);
    }

    #[test]
    fn spme_mesh_rejects_overflow() {
        let lattice = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0];
        assert_eq!(spme_mesh(&lattice, 17, 10.0), Ok(([32; 3], 32768)));
        assert_eq!(spme_mesh(&lattice, 0, 10.0), Ok(([32; 3], 32768)));
        assert!(spme_mesh(&lattice, i32::MAX, 10.0).is_err());
        assert!(spme_mesh(&lattice, 0, f64::MAX).is_err());
        assert!(spme_mesh(&lattice, 0, 1e10).is_err());
    }

    #[test]
    fn spme_matches_cubic_reference() {
        let numbers = [6, 8];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
        let lattice = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0];
        let damping = Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        };
        let evaluate = |mesh| {
            periodic_energy(
                &numbers,
                &positions,
                &lattice,
                damping,
                EwaldConfig {
                    rank: 10_000,
                    tolerance: 1.0e-4,
                    kcut: 10.0,
                    mesh,
                },
            )
            .unwrap()
        };
        let fixed = evaluate(64);
        let automatic = evaluate(0);
        assert!((fixed - -0.0028149228007021015).abs() < 1.0e-11, "{fixed}");
        assert!(
            (automatic - -0.0028149227957823923).abs() < 1.0e-10,
            "{automatic}"
        );
        let derivatives = |mesh| {
            periodic_derivatives_partitioned(
                &numbers,
                &positions,
                &lattice,
                damping,
                EwaldConfig {
                    rank: 10_000,
                    tolerance: 1.0e-4,
                    kcut: 10.0,
                    mesh,
                },
                WorkPartition::SERIAL,
            )
            .unwrap()
        };
        let fixed = derivatives(64);
        assert!((fixed.energy - -0.0028149228007021015).abs() < 1.0e-11);
        assert!((fixed.gradient[0] - 1.1237229657839283e-12).abs() < 1.0e-10);
        assert!((fixed.gradient[1] - 2.9396680396406752e-5).abs() < 1.0e-10);
        assert!((fixed.gradient[2] - 5.2457765748044166e-5).abs() < 1.0e-10);
        assert!((fixed.virial[0] - 0.0037461047321046635).abs() < 1.0e-9);
        assert!((fixed.virial[4] - 0.003644616811269228).abs() < 1.0e-9);
        assert!((fixed.virial[5] - -0.0002457998363918225).abs() < 1.0e-9);
        let automatic = derivatives(0);
        assert!((automatic.energy - -0.0028149227957823923).abs() < 1.0e-10);
        assert!((automatic.gradient[1] - 2.939671231545548e-5).abs() < 1.0e-9);
        assert!((automatic.virial[0] - 0.0037461048190068438).abs() < 1.0e-8);
    }

    #[test]
    fn fourier_partitions_sum_to_serial() {
        let numbers = [6, 8];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
        let lattice = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0];
        let damping = Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        };
        for mesh in [-1, 32] {
            let config = EwaldConfig {
                rank: 10_000,
                tolerance: 1.0e-4,
                kcut: 6.0,
                mesh,
            };
            let serial = periodic_energy(&numbers, &positions, &lattice, damping, config).unwrap();
            let partitioned = (0..3)
                .map(|part| {
                    periodic_energy_partitioned(
                        &numbers,
                        &positions,
                        &lattice,
                        damping,
                        config,
                        WorkPartition::new(part, 3).unwrap(),
                    )
                    .unwrap()
                })
                .sum::<f64>();
            assert!((partitioned - serial).abs() < 1.0e-13);
            let serial = periodic_derivatives_partitioned(
                &numbers,
                &positions,
                &lattice,
                damping,
                config,
                WorkPartition::SERIAL,
            )
            .unwrap();
            let parts: Vec<_> = (0..3)
                .map(|part| {
                    periodic_derivatives_partitioned(
                        &numbers,
                        &positions,
                        &lattice,
                        damping,
                        config,
                        WorkPartition::new(part, 3).unwrap(),
                    )
                    .unwrap()
                })
                .collect();
            assert!(
                (parts.iter().map(|part| part.energy).sum::<f64>() - serial.energy).abs() < 1.0e-12
            );
            for index in 0..serial.gradient.len() {
                assert!(
                    (parts.iter().map(|part| part.gradient[index]).sum::<f64>()
                        - serial.gradient[index])
                        .abs()
                        < 1.0e-11
                );
            }
            for index in 0..9 {
                assert!(
                    (parts.iter().map(|part| part.virial[index]).sum::<f64>()
                        - serial.virial[index])
                        .abs()
                        < 1.0e-10
                );
            }
        }
    }

    #[test]
    fn zero_damping_rejects_unrepresentable_exponents() {
        let zero = |alpha| Damping::Zero {
            s6: 1.0,
            s8: 0.722,
            rs6: 1.217,
            rs8: 1.0,
            alpha,
        };
        for alpha in [4.0, 14.0, 16.0, 64.0] {
            assert!(fourier_terms(zero(alpha), 5, 7).is_ok(), "alpha {alpha}");
        }
        // Exponents that would overflow the i32 arithmetic or explode the pole sum.
        for alpha in [3.0, 65.0, 2147483648.0, 1.0e30, f64::MAX, -14.0] {
            assert!(fourier_terms(zero(alpha), 5, 7).is_err(), "alpha {alpha}");
        }
    }

    #[test]
    fn zero_damping_periodic_derivatives_match_reference() {
        let numbers = [6, 8];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
        let lattice = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0];
        let damping = Damping::Zero {
            s6: 1.0,
            s8: 0.722,
            rs6: 1.217,
            rs8: 1.0,
            alpha: 14.0,
        };
        for (mesh, reference) in [
            (
                -1,
                (
                    -0.002736274482112929,
                    0.00032704430213445804,
                    0.0003747814029760236,
                    0.0015928119650939646,
                    0.0031855499202944016,
                    -0.0009357420563380467,
                ),
            ),
            (
                64,
                (
                    -0.0027362534329724675,
                    0.00032700580544871955,
                    0.0003748431231891855,
                    0.00159185457126035,
                    0.0031844462953180097,
                    -0.0009355922365882624,
                ),
            ),
        ] {
            let result = periodic_derivatives_partitioned(
                &numbers,
                &positions,
                &lattice,
                damping,
                EwaldConfig {
                    rank: 10_000,
                    tolerance: 1.0e-4,
                    kcut: 10.0,
                    mesh,
                },
                WorkPartition::SERIAL,
            )
            .unwrap();
            assert!((result.energy - reference.0).abs() < 1.0e-10);
            assert!((result.gradient[1] - reference.1).abs() < 1.0e-9);
            assert!((result.gradient[2] - reference.2).abs() < 1.0e-9);
            assert!((result.virial[0] - reference.3).abs() < 1.0e-8);
            assert!((result.virial[4] - reference.4).abs() < 1.0e-8);
            assert!((result.virial[5] - reference.5).abs() < 1.0e-8);
        }
    }

    #[test]
    fn periodic_ghost_derivatives_match_reference() {
        let numbers = [6, 8];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
        let lattice = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0];
        let damping = Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        };
        for (mesh, reference) in [
            (
                -1,
                (
                    -0.0002548699546013089,
                    8.463014652590683e-8,
                    9.980536834115003e-8,
                    0.000464682300116663,
                ),
            ),
            (
                64,
                (
                    -0.0002548699545715037,
                    8.463014651600948e-8,
                    9.980536832947791e-8,
                    0.0004646822993399868,
                ),
            ),
        ] {
            let result = periodic_derivatives_partitioned_with_ghosts(
                &numbers,
                &positions,
                &lattice,
                damping,
                EwaldConfig {
                    rank: 10_000,
                    tolerance: 1.0e-4,
                    kcut: 10.0,
                    mesh,
                },
                &[true, false],
                WorkPartition::SERIAL,
            )
            .unwrap();
            assert!((result.energy - reference.0).abs() < 1.0e-10);
            assert!((result.gradient[1] - reference.1).abs() < 1.0e-10);
            assert!((result.gradient[2] - reference.2).abs() < 1.0e-10);
            assert!((result.virial[0] - reference.3).abs() < 1.0e-9);
        }
    }

    #[test]
    fn periodic_cn_cutoff_matches_reference() {
        let numbers = [6, 8];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
        let lattice = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0];
        let damping = Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        };
        for (mesh, reference) in [
            (
                -1,
                (
                    -0.0028160243304364476,
                    2.8620570755438185e-5,
                    5.154974134251637e-5,
                    0.0037543250300080327,
                ),
            ),
            (
                64,
                (
                    -0.002816024330369562,
                    2.8620571178024156e-5,
                    5.1549741079440304e-5,
                    0.0037543250306542285,
                ),
            ),
        ] {
            let result = periodic_derivatives_partitioned_with_cutoff(
                &numbers,
                &positions,
                &lattice,
                damping,
                EwaldConfig {
                    rank: 10_000,
                    tolerance: 1.0e-4,
                    kcut: 10.0,
                    mesh,
                },
                RealspaceCutoff {
                    cn: 2.0,
                    ..RealspaceCutoff::default()
                },
                &[],
                WorkPartition::SERIAL,
            )
            .unwrap();
            assert!((result.energy - reference.0).abs() < 1.0e-10);
            assert!((result.gradient[1] - reference.1).abs() < 1.0e-10);
            assert!((result.gradient[2] - reference.2).abs() < 1.0e-10);
            assert!((result.virial[0] - reference.3).abs() < 1.0e-9);
        }
    }
}
