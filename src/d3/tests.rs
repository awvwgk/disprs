use super::*;

use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::dual::Dual2;

pub fn energy(numbers: &[i32], positions: &[f64], damping: Damping) -> Result<f64, &'static str> {
    energy_partitioned(numbers, positions, damping, WorkPartition::SERIAL)
}

pub fn gradient(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
) -> Result<(f64, Vec<f64>, [f64; 9]), &'static str> {
    gradient_partitioned(numbers, positions, damping, WorkPartition::SERIAL)
}

pub fn hessian(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
) -> Result<(f64, Vec<f64>), &'static str> {
    hessian_partitioned(numbers, positions, damping, WorkPartition::SERIAL)
}

fn dual_hessian_partitioned_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>), &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, species, _, _) = prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let size = positions.len();
    let coordinates: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(index, &value)| Dual2::variable(value, size, index))
        .collect();
    let mut coordination = vec![Dual2::constant(0.0, size); numbers.len()];
    for first in 0..numbers.len() {
        for second in 0..first {
            let distance2 = dual_distance2(&coordinates, first, second);
            if distance2.value > cutoff.cn * cutoff.cn || distance2.value < 1.0e-12 {
                continue;
            }
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (((radius / distance2.sqrt() - 1.0) * -16.0).exp() + 1.0);
            coordination[first] += count.clone();
            coordination[second] += count;
        }
    }

    let weights: Vec<Vec<Dual2>> = elements
        .iter()
        .zip(coordination.iter())
        .map(|(&element, coordination)| {
            let count = reference_count(element);
            let raw: Vec<_> = (0..count)
                .map(|reference| {
                    let delta = coordination.clone()
                        - value(REFERENCE_CN, element * REFERENCES + reference);
                    (delta.clone() * delta * -4.0).exp()
                })
                .collect();
            let norm = raw
                .iter()
                .cloned()
                .fold(Dual2::constant(0.0, size), |sum, item| sum + item);
            raw.into_iter().map(|item| item / norm.clone()).collect()
        })
        .collect();

    let mut total = Dual2::constant(0.0, size);
    for first in 0..numbers.len() {
        for second in 0..first {
            if ghosts.get(first) == Some(&true)
                || ghosts.get(second) == Some(&true)
                || !partition.owns_pair(first, second)
            {
                continue;
            }
            let distance2 = dual_distance2(&coordinates, first, second);
            if distance2.value > cutoff.disp2 * cutoff.disp2 || distance2.value <= f64::EPSILON {
                continue;
            }
            let mut c6 = Dual2::constant(0.0, size);
            for first_reference in 0..reference_count(elements[first]) {
                for second_reference in 0..reference_count(elements[second]) {
                    c6 += weights[first][first_reference].clone()
                        * weights[second][second_reference].clone()
                        * reference_c6(
                            elements[first],
                            elements[second],
                            first_reference,
                            second_reference,
                        );
                }
            }
            let switch = dual_smooth_cutoff(distance2.clone(), cutoff.disp2, cutoff.width2);
            total += dual_kernel(
                damping,
                elements[first],
                elements[second],
                species[first],
                species[second],
                distance2,
                c6,
            ) * switch;
        }
    }
    Ok((total.value, total.hessian))
}

fn dual_distance2(coordinates: &[Dual2], first: usize, second: usize) -> Dual2 {
    (0..3)
        .map(|axis| coordinates[3 * first + axis].clone() - coordinates[3 * second + axis].clone())
        .map(|delta| delta.clone() * delta)
        .fold(Dual2::constant(0.0, coordinates.len()), |sum, item| {
            sum + item
        })
}

fn dual_smooth_cutoff(distance2: Dual2, cutoff: f64, width: f64) -> Dual2 {
    let size = distance2.gradient.len();
    let distance = distance2.sqrt();
    if width <= 0.0 || width >= cutoff || distance.value <= cutoff - width {
        Dual2::constant(1.0, size)
    } else if distance.value >= cutoff {
        Dual2::constant(0.0, size)
    } else {
        let x = (Dual2::constant(cutoff, size) - distance) / width;
        x.clone().powf(3.0)
            * (Dual2::constant(10.0, size) + x.clone() * (Dual2::constant(-15.0, size) + x * 6.0))
    }
}

fn dual_kernel(
    damping: Damping,
    first_element: usize,
    second_element: usize,
    first_species: usize,
    second_species: usize,
    u: Dual2,
    c6: Dual2,
) -> Dual2 {
    let size = u.gradient.len();
    let rrij = 3.0 * value(R4R2, first_element) * value(R4R2, second_element);
    let potential = match damping {
        Damping::Rational { s6, s8, a1, a2 } => {
            let r0 = a1 * rrij.sqrt() + a2;
            s6 / (u.clone().powf(3.0) + r0.powi(6)) + (s8 * rrij) / (u.powf(4.0) + r0.powi(8))
        }
        Damping::Zero {
            s6,
            s8,
            rs6,
            rs8,
            alpha,
        } => {
            let r0 = value(VDW_RADII, pair_index(first_element, second_element));
            let f6 = 1.0 / (((rs6 * r0) / u.clone().sqrt()).powf(alpha) * 6.0 + 1.0);
            let f8 = 1.0 / (((rs8 * r0) / u.clone().sqrt()).powf(alpha + 2.0) * 6.0 + 1.0);
            f6 * s6 / u.clone().powf(3.0) + f8 * (s8 * rrij) / u.powf(4.0)
        }
        Damping::ModifiedZero {
            s6,
            s8,
            rs6,
            rs8,
            alpha,
            beta,
        } => {
            let r0 = value(VDW_RADII, pair_index(first_element, second_element));
            let distance = u.clone().sqrt();
            let f6 = Dual2::constant(1.0, size)
                / ((distance.clone() / (rs6 * r0) + beta * r0).powf(-alpha) * 6.0 + 1.0);
            let f8 = Dual2::constant(1.0, size)
                / ((distance / (rs8 * r0) + beta * r0).powf(-alpha - 2.0) * 6.0 + 1.0);
            f6 * s6 / u.clone().powf(3.0) + f8 * (s8 * rrij) / u.powf(4.0)
        }
        Damping::OptimizedPower {
            s6,
            s8,
            a1,
            a2,
            beta,
        } => {
            let r0 = a1 * rrij.sqrt() + a2;
            let d6 = u.clone().powf(3.0) + u.clone().powf(-0.5 * beta) * r0.powf(6.0 + beta);
            let d8 = u.clone().powf(4.0) + u.powf(-0.5 * beta) * r0.powf(8.0 + beta);
            s6 / d6 + (s8 * rrij) / d8
        }
        Damping::Cso { s6, a1, a2, a3, a4 } => {
            let r0 = rrij.sqrt();
            let sigmoid = 1.0 / ((u.clone().sqrt() - a2 * r0).exp() + 1.0);
            let scale = sigmoid * a1 + s6;
            scale / (u.powf(3.0) + (a3 * r0 + a4).powi(6))
        }
        Damping::Z { s6, s8, a1 } => {
            let r0 = a1 / (first_species + second_species + 2) as f64;
            s6 / (u.clone().powf(3.0) + c6.clone() * r0)
                + (s8 * rrij) / (u.powf(4.0) + c6.clone() * (r0 * rrij))
        }
    };
    -(c6 * potential)
}

impl LocalDual {
    fn constant(value: f64) -> Self {
        Self {
            value,
            gradient: [0.0; 6],
            hessian: [0.0; 36],
        }
    }

    fn variable(value: f64, index: usize) -> Self {
        let mut result = Self::constant(value);
        result.gradient[index] = 1.0;
        result
    }

    fn powf(mut self, exponent: f64) -> Self {
        let value = self.value.powf(exponent);
        let first = exponent * self.value.powf(exponent - 1.0);
        let second = exponent * (exponent - 1.0) * self.value.powf(exponent - 2.0);
        for index in 0..36 {
            self.hessian[index] = first * self.hessian[index]
                + second * self.gradient[index / 6] * self.gradient[index % 6];
        }
        for item in &mut self.gradient {
            *item *= first;
        }
        self.value = value;
        self
    }

    fn sqrt(self) -> Self {
        self.powf(0.5)
    }
}

impl Add for LocalDual {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self {
        self.value += rhs.value;
        for index in 0..6 {
            self.gradient[index] += rhs.gradient[index];
        }
        for index in 0..36 {
            self.hessian[index] += rhs.hessian[index];
        }
        self
    }
}

impl Add<f64> for LocalDual {
    type Output = Self;
    fn add(mut self, rhs: f64) -> Self {
        self.value += rhs;
        self
    }
}

impl Neg for LocalDual {
    type Output = Self;
    fn neg(mut self) -> Self {
        self.value = -self.value;
        self.gradient.iter_mut().for_each(|value| *value = -*value);
        self.hessian.iter_mut().for_each(|value| *value = -*value);
        self
    }
}

impl Sub for LocalDual {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl Mul<f64> for LocalDual {
    type Output = Self;
    fn mul(mut self, rhs: f64) -> Self {
        self.value *= rhs;
        self.gradient.iter_mut().for_each(|value| *value *= rhs);
        self.hessian.iter_mut().for_each(|value| *value *= rhs);
        self
    }
}

impl Mul for LocalDual {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self {
        let left_value = self.value;
        for index in 0..36 {
            let row = index / 6;
            let column = index % 6;
            self.hessian[index] = self.hessian[index] * rhs.value
                + self.gradient[row] * rhs.gradient[column]
                + self.gradient[column] * rhs.gradient[row]
                + left_value * rhs.hessian[index];
        }
        for index in 0..6 {
            self.gradient[index] =
                self.gradient[index] * rhs.value + left_value * rhs.gradient[index];
        }
        self.value *= rhs.value;
        self
    }
}

impl Div for LocalDual {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self * rhs.powf(-1.0)
    }
}

impl Div<LocalDual> for f64 {
    type Output = LocalDual;
    fn div(self, rhs: LocalDual) -> LocalDual {
        rhs.powf(-1.0) * self
    }
}

fn local_atm_dual(
    distances: [f64; 3],
    c6: [f64; 3],
    radii: [f64; 3],
    atm: Atm,
    cutoff: RealspaceCutoff,
) -> LocalDual {
    let u: [LocalDual; 3] =
        std::array::from_fn(|index| LocalDual::variable(distances[index], index));
    let c: [LocalDual; 3] = std::array::from_fn(|index| LocalDual::variable(c6[index], 3 + index));
    let product = u[0].clone() * u[1].clone() * u[2].clone();
    let angular = (u[0].clone() + u[2].clone() - u[1].clone())
        * (u[0].clone() - u[2].clone() + u[1].clone())
        * (-u[0].clone() + u[2].clone() + u[1].clone())
        * 0.375
        / product.clone().powf(2.5)
        + 1.0 / product.clone().powf(1.5);
    let r0 = (4.0_f64 / 3.0).powi(3) * radii.iter().product::<f64>();
    let damping = 1.0 / ((r0 / product.sqrt()).powf(atm.alpha / 3.0) * 6.0 + 1.0);
    let switch = local_smooth_cutoff(u[0].clone(), cutoff.disp3, cutoff.width3)
        * local_smooth_cutoff(u[1].clone(), cutoff.disp3, cutoff.width3)
        * local_smooth_cutoff(u[2].clone(), cutoff.disp3, cutoff.width3);
    (c[0].clone() * c[1].clone() * c[2].clone()).sqrt() * angular * damping * atm.s9 * switch
}

fn local_smooth_cutoff(distance2: LocalDual, cutoff: f64, width: f64) -> LocalDual {
    let distance = distance2.sqrt();
    if width <= 0.0 || width >= cutoff || distance.value <= cutoff - width {
        LocalDual::constant(1.0)
    } else if distance.value >= cutoff {
        LocalDual::constant(0.0)
    } else {
        let x = (LocalDual::constant(cutoff) - distance) * width.recip();
        x.clone().powf(3.0)
            * (LocalDual::constant(10.0) + x.clone() * (LocalDual::constant(-15.0) + x * 6.0))
    }
}

fn dual_atm_derivatives_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<AtmDerivatives, &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, _, _, _) = prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let size = positions.len();
    let coordinates: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(index, &value)| Dual2::variable(value, size, index))
        .collect();
    let weights = dual_weights(&elements, &coordinates, cutoff.cn);
    let mut total = Dual2::constant(0.0, size);
    for first in 0..numbers.len() {
        for second in 0..first {
            if !partition.owns_pair(first, second) {
                continue;
            }
            for third in 0..second {
                if ghosts.get(first) == Some(&true)
                    || ghosts.get(second) == Some(&true)
                    || ghosts.get(third) == Some(&true)
                {
                    continue;
                }
                let u12 = dual_distance2(&coordinates, first, second);
                let u13 = dual_distance2(&coordinates, first, third);
                let u23 = dual_distance2(&coordinates, second, third);
                if [u12.value, u13.value, u23.value].iter().any(|&distance| {
                    !(f64::EPSILON..=cutoff.disp3 * cutoff.disp3).contains(&distance)
                }) {
                    continue;
                }
                let c12 = dual_c6(&elements, &weights, first, second, size);
                let c13 = dual_c6(&elements, &weights, first, third, size);
                let c23 = dual_c6(&elements, &weights, second, third, size);
                let product = u12.clone() * u13.clone() * u23.clone();
                let r0 = (4.0_f64 / 3.0).powi(3)
                    * value(VDW_RADII, pair_index(elements[first], elements[second]))
                    * value(VDW_RADII, pair_index(elements[first], elements[third]))
                    * value(VDW_RADII, pair_index(elements[second], elements[third]));
                let switch = dual_smooth_cutoff(u12.clone(), cutoff.disp3, cutoff.width3)
                    * dual_smooth_cutoff(u13.clone(), cutoff.disp3, cutoff.width3)
                    * dual_smooth_cutoff(u23.clone(), cutoff.disp3, cutoff.width3);
                let angular = (u12.clone() + u23.clone() - u13.clone())
                    * (u12.clone() - u23.clone() + u13.clone())
                    * (-u12 + u23 + u13)
                    * 0.375
                    / product.clone().powf(2.5)
                    + 1.0 / product.clone().powf(1.5);
                let damping = 1.0 / ((r0 / product.sqrt()).powf(atm.alpha / 3.0) * 6.0 + 1.0);
                total += (c12 * c13 * c23).sqrt() * angular * damping * atm.s9 * switch;
            }
        }
    }
    let mut virial = [0.0; 9];
    for atom in 0..numbers.len() {
        for axis in 0..3 {
            for other in 0..3 {
                virial[axis + 3 * other] +=
                    total.gradient[3 * atom + axis] * positions[3 * atom + other];
            }
        }
    }
    Ok((total.value, total.gradient, virial, total.hessian))
}

fn dual_weights(elements: &[usize], coordinates: &[Dual2], cn_cutoff: f64) -> Vec<Vec<Dual2>> {
    let size = coordinates.len();
    let mut coordination = vec![Dual2::constant(0.0, size); elements.len()];
    for first in 0..elements.len() {
        for second in 0..first {
            let distance2 = dual_distance2(coordinates, first, second);
            if distance2.value > cn_cutoff * cn_cutoff || distance2.value < 1.0e-12 {
                continue;
            }
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (((radius / distance2.sqrt() - 1.0) * -16.0).exp() + 1.0);
            coordination[first] += count.clone();
            coordination[second] += count;
        }
    }
    elements
        .iter()
        .zip(coordination)
        .map(|(&element, coordination)| {
            let raw: Vec<_> = (0..reference_count(element))
                .map(|reference| {
                    let delta = coordination.clone()
                        - value(REFERENCE_CN, element * REFERENCES + reference);
                    (delta.clone() * delta * -4.0).exp()
                })
                .collect();
            let norm = raw
                .iter()
                .cloned()
                .fold(Dual2::constant(0.0, size), |sum, item| sum + item);
            raw.into_iter().map(|item| item / norm.clone()).collect()
        })
        .collect()
}

fn dual_c6(
    elements: &[usize],
    weights: &[Vec<Dual2>],
    first: usize,
    second: usize,
    size: usize,
) -> Dual2 {
    let mut result = Dual2::constant(0.0, size);
    for first_reference in 0..reference_count(elements[first]) {
        for second_reference in 0..reference_count(elements[second]) {
            result += weights[first][first_reference].clone()
                * weights[second][second_reference].clone()
                * reference_c6(
                    elements[first],
                    elements[second],
                    first_reference,
                    second_reference,
                );
        }
    }
    result
}

#[test]
fn properties_follow_geometry_cutoffs_and_periodic_images() {
    let numbers = [6, 8, 1];
    let positions = [1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.5, 3.0, 1.5];
    let lattice = [6.0, 0.0, 0.0, 1.0, 7.0, 0.0, 0.5, 0.2, 8.0];
    let (cn, c6) = properties(&numbers, &positions, None, [false; 3], 7.0).unwrap();
    let (elements, _, expected_cn, weights) =
        prepare_with_cutoff(&numbers, &positions, 7.0).unwrap();
    assert_eq!(cn, expected_cn);
    for first in 0..3 {
        for second in 0..3 {
            assert!(
                (c6[first * 3 + second]
                    - atomic_c6(
                        elements[first],
                        elements[second],
                        &weights[first],
                        &weights[second],
                    ))
                .abs()
                    < 1e-12
            );
            assert_eq!(c6[first * 3 + second], c6[second * 3 + first]);
        }
    }
    assert_eq!(
        properties(&numbers, &positions, None, [false; 3], 0.0)
            .unwrap()
            .0,
        [0.0; 3]
    );
    for dimensions in 1..=3 {
        let periodic = std::array::from_fn(|axis| axis < dimensions);
        let reference = properties(&numbers, &positions, Some(&lattice), periodic, 7.0).unwrap();
        assert!(reference
            .0
            .iter()
            .zip(&cn)
            .any(|(periodic, molecular)| periodic > molecular));
        let mut shifted = positions;
        for axis in 0..3 {
            shifted[axis] += 3.0 * lattice[axis];
        }
        let translated = properties(&numbers, &shifted, Some(&lattice), periodic, 7.0).unwrap();
        for (original, moved) in reference
            .0
            .iter()
            .chain(&reference.1)
            .zip(translated.0.iter().chain(&translated.1))
        {
            assert!((original - moved).abs() < 1e-11);
        }
    }
    assert!(properties(&numbers, &positions[..8], None, [false; 3], 7.0).is_err());
    assert!(properties(&numbers, &positions, None, [true; 3], 7.0).is_err());
    assert!(properties(&numbers, &positions, None, [false; 3], f64::NAN).is_err());
}

#[test]
fn analytical_kernel_second_derivatives_match_dual() {
    let dampings = [
        Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        },
        Damping::Zero {
            s6: 1.0,
            s8: 0.722,
            rs6: 1.217,
            rs8: 1.0,
            alpha: 14.0,
        },
        Damping::ModifiedZero {
            s6: 1.0,
            s8: 0.5,
            rs6: 1.2,
            rs8: 1.1,
            alpha: 14.0,
            beta: 0.1,
        },
        Damping::OptimizedPower {
            s6: 1.0,
            s8: 0.5,
            a1: 0.4,
            a2: 4.5,
            beta: 2.0,
        },
        Damping::Cso {
            s6: 1.0,
            a1: 0.24,
            a2: 2.5,
            a3: 0.0,
            a4: 6.25,
        },
        Damping::Z {
            s6: 1.0,
            s8: 1.0,
            a1: 1.0,
        },
    ];
    for damping in dampings {
        let analytical = kernel_second_derivatives(damping, 5, 7, 0, 1, 25.0, 30.0);
        let dual = dual_kernel(
            damping,
            5,
            7,
            0,
            1,
            Dual2::variable(25.0, 2, 0),
            Dual2::variable(30.0, 2, 1),
        );
        let expected = [
            dual.value,
            dual.gradient[0],
            dual.gradient[1],
            dual.hessian[0],
            dual.hessian[1],
            dual.hessian[3],
        ];
        for (actual, expected) in analytical.iter().zip(expected) {
            assert!(
                (actual - expected).abs() < 1.0e-12,
                "{actual} != {expected}"
            );
        }
    }
}

#[test]
fn analytical_hessian_matches_dual_oracle() {
    let numbers = [6, 8, 7];
    let positions = [0.2, 0.4, 0.7, 4.8, 0.3, 0.5, 1.1, 4.2, 0.9];
    let dampings = [
        Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        },
        Damping::Zero {
            s6: 1.0,
            s8: 0.722,
            rs6: 1.217,
            rs8: 1.0,
            alpha: 14.0,
        },
        Damping::ModifiedZero {
            s6: 1.0,
            s8: 0.5,
            rs6: 1.2,
            rs8: 1.1,
            alpha: 14.0,
            beta: 0.1,
        },
        Damping::OptimizedPower {
            s6: 1.0,
            s8: 0.5,
            a1: 0.4,
            a2: 4.5,
            beta: 2.0,
        },
        Damping::Cso {
            s6: 1.0,
            a1: 0.24,
            a2: 2.5,
            a3: 0.0,
            a4: 6.25,
        },
        Damping::Z {
            s6: 1.0,
            s8: 1.0,
            a1: 1.0,
        },
    ];
    for damping in dampings {
        for (cutoff, ghosts) in [
            (RealspaceCutoff::default(), &[false, false, false][..]),
            (
                RealspaceCutoff {
                    disp2: 6.0,
                    width2: 2.0,
                    ..RealspaceCutoff::default()
                },
                &[true, false, false][..],
            ),
        ] {
            let analytical = analytical_hessian(
                &numbers,
                &positions,
                damping,
                cutoff,
                ghosts,
                WorkPartition::SERIAL,
            )
            .unwrap();
            let dual = dual_hessian_partitioned_with_cutoff(
                &numbers,
                &positions,
                damping,
                cutoff,
                ghosts,
                WorkPartition::SERIAL,
            )
            .unwrap();
            assert!((analytical.0 - dual.0).abs() < 1.0e-13);
            for (actual, expected) in analytical.1.iter().zip(dual.1) {
                assert!(
                    (actual - expected).abs() < 2.0e-11,
                    "{actual} != {expected}"
                );
            }
        }
    }
}

#[test]
fn local_atm_partials_match_dual() {
    let cutoff = RealspaceCutoff {
        disp3: 7.0,
        width3: 2.0,
        ..RealspaceCutoff::default()
    };
    for distances in [
        [9.0, 16.0, 25.0],
        [9.0, 9.0, 36.0],
        [30.0, 35.0, 40.0],
        [25.0, 36.0, 49.0],
    ] {
        let atm = Atm {
            s9: 0.7,
            alpha: 14.0,
        };
        let actual = local_atm_partials(distances, [12.0, 18.0, 24.0], [3.0; 3], atm, cutoff);
        let expected = local_atm_dual(distances, [12.0, 18.0, 24.0], [3.0; 3], atm, cutoff);
        for (actual, expected) in std::iter::once(actual.value)
            .chain(actual.gradient)
            .chain(actual.hessian)
            .zip(
                std::iter::once(expected.value)
                    .chain(expected.gradient)
                    .chain(expected.hessian),
            )
        {
            assert!(
                (actual - expected).abs() < 1.0e-14,
                "{actual} != {expected}"
            );
        }
    }
}

#[test]
fn analytical_atm_hessian_matches_dual_oracle() {
    let numbers = [6, 8, 7, 1];
    let positions = [0.2, 0.4, 0.7, 4.8, 0.3, 0.5, 1.1, 4.2, 0.9, 3.3, 3.6, 4.1];
    let atm = Atm {
        s9: 1.0,
        alpha: 14.0,
    };
    for (cutoff, ghosts) in [
        (
            RealspaceCutoff::default(),
            &[false, false, false, false][..],
        ),
        (
            RealspaceCutoff {
                disp3: 7.0,
                width3: 2.0,
                ..RealspaceCutoff::default()
            },
            &[true, false, false, false][..],
        ),
    ] {
        for partition in [WorkPartition::SERIAL, WorkPartition::new(1, 2).unwrap()] {
            let analytical = analytical_atm_derivatives(
                &numbers, &positions, atm, cutoff, ghosts, partition, true,
            )
            .unwrap();
            let dual = dual_atm_derivatives_with_cutoff(
                &numbers, &positions, atm, cutoff, ghosts, partition,
            )
            .unwrap();
            assert!((analytical.0 - dual.0).abs() < 1.0e-13);
            for (actual, expected) in analytical.1.iter().zip(dual.1) {
                assert!(
                    (actual - expected).abs() < 2.0e-12,
                    "{actual} != {expected}"
                );
            }
            for (actual, expected) in analytical.2.iter().zip(dual.2) {
                assert!(
                    (actual - expected).abs() < 2.0e-12,
                    "{actual} != {expected}"
                );
            }
            for (actual, expected) in analytical.3.iter().zip(dual.3) {
                assert!(
                    (actual - expected).abs() < 2.0e-10,
                    "{actual} != {expected}"
                );
            }
        }
    }
}

#[test]
fn matches_pbe_d3_bj_reference() {
    let energy = energy(
        &[6, 6],
        &[0.0, 0.0, 0.0, 6.0, 0.0, 0.0],
        Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        },
    )
    .unwrap();
    assert!((energy - -0.000_534_141_393_133_826_7).abs() < 1.0e-15);
}

#[test]
fn ghost_atoms_only_affect_coordination() {
    let numbers = [6, 8, 1];
    let positions = [0.0, 0.0, 0.0, 2.1, 0.3, 0.1, 1.0, 1.7, 0.4];
    let damping = Damping::Rational {
        s6: 1.0,
        s8: 0.7875,
        a1: 0.4289,
        a2: 4.4407,
    };
    let ghosts = [true, false, false];
    let (energy, gradient, _) = gradient_partitioned_with_ghosts(
        &numbers,
        &positions,
        damping,
        &ghosts,
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((energy - -0.0001337981636225901).abs() < 1.0e-15);
    for (actual, reference) in gradient.iter().zip([
        1.2298897985897666e-9,
        4.179935613941928e-10,
        1.1339746952941622e-10,
        4.493982480184129e-8,
        -5.871688528980564e-8,
        -1.2600459900417153e-8,
        -4.616971460043105e-8,
        5.829889172841145e-8,
        1.248706243088774e-8,
    ]) {
        assert!((actual - reference).abs() < 1.0e-15);
    }
    let pairwise = pairwise_partitioned_with_ghosts(
        &numbers,
        &positions,
        damping,
        &ghosts,
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!(pairwise[..3].iter().all(|&value| value == 0.0));
    assert!((pairwise.iter().sum::<f64>() - energy).abs() < 1.0e-15);
    let (hessian_energy, _) = hessian_partitioned_with_ghosts(
        &numbers,
        &positions,
        damping,
        &ghosts,
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((hessian_energy - energy).abs() < 1.0e-15);
    let atm = Atm {
        s9: 1.0,
        alpha: 14.0,
    };
    assert!(
        atm_pairwise_with_ghosts(&numbers, &positions, atm, &ghosts, WorkPartition::SERIAL,)
            .unwrap()
            .iter()
            .all(|&value| value == 0.0)
    );
    assert_eq!(
        atm_derivatives_with_ghosts(&numbers, &positions, atm, &ghosts, WorkPartition::SERIAL,)
            .unwrap()
            .0,
        0.0
    );
    assert!(energy_partitioned_with_ghosts(
        &numbers,
        &positions,
        damping,
        &[true],
        WorkPartition::SERIAL,
    )
    .is_err());
}

#[test]
fn smooth_two_body_cutoff_matches_reference() {
    let cutoff = RealspaceCutoff {
        disp2: 6.0,
        width2: 2.0,
        ..RealspaceCutoff::default()
    };
    let (energy, gradient, virial) = gradient_partitioned_with_cutoff(
        &[6, 6],
        &[0.0, 0.0, 0.0, 5.5, 0.0, 0.0],
        Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        },
        cutoff,
        &[],
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((energy - -6.405046062946415e-5).abs() < 1.0e-15);
    assert!((gradient[0] - -3.4189782676699354e-4).abs() < 1.0e-15);
    assert!((virial[0] - 1.8804380472184643e-3).abs() < 1.0e-15);
    let (hessian_energy, hessian) = hessian_partitioned_with_cutoff(
        &[6, 6],
        &[0.0, 0.0, 0.0, 5.5, 0.0, 0.0],
        Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        },
        cutoff,
        &[],
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((hessian_energy - energy).abs() < 1.0e-15);
    assert!((hessian[0] - -1.0205240365877704e-3).abs() < 1.0e-14);
    assert!((hessian[3] - 1.0205240365877704e-3).abs() < 1.0e-14);
    assert!((hessian[7] - 6.216324123036246e-5).abs() < 1.0e-14);
}

#[test]
fn smooth_atm_cutoff_matches_reference() {
    let numbers = [6, 6, 6];
    let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 2.5, 4.5, 0.0];
    let cutoff = RealspaceCutoff {
        disp3: 6.0,
        width3: 2.0,
        ..RealspaceCutoff::default()
    };
    let result = atm_derivatives_with_cutoff(
        &numbers,
        &positions,
        Atm {
            s9: 1.0,
            alpha: 16.0,
        },
        cutoff,
        &[],
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((result.0 - 6.630248136070531e-9).abs() < 1.0e-18);
    assert!((result.1[0] - 1.5761211815267338e-8).abs() < 1.0e-17);
    assert!((result.1[1] - 1.1709141438805946e-8).abs() < 1.0e-17);
    assert!((result.2[0] - -7.880605907633669e-8).abs() < 1.0e-17);
    let gradient = atm_gradient_with_cutoff(
        &numbers,
        &positions,
        Atm {
            s9: 1.0,
            alpha: 16.0,
        },
        cutoff,
        &[],
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((gradient.0 - result.0).abs() < 1.0e-18);
    for (actual, expected) in gradient.1.iter().zip(result.1.iter()) {
        assert!((actual - expected).abs() < 1.0e-17);
    }
    for (actual, expected) in gradient.2.iter().zip(result.2.iter()) {
        assert!((actual - expected).abs() < 1.0e-17);
    }
    let pairwise = atm_pairwise_with_cutoff(
        &numbers,
        &positions,
        Atm {
            s9: 1.0,
            alpha: 16.0,
        },
        cutoff,
        &[],
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((pairwise.iter().sum::<f64>() - result.0).abs() < 1.0e-18);
}

#[test]
fn zero_damping_is_finite_and_attractive() {
    let energy = energy(
        &[6, 6],
        &[0.0, 0.0, 0.0, 6.0, 0.0, 0.0],
        Damping::Zero {
            s6: 1.0,
            s8: 0.722,
            rs6: 1.217,
            rs8: 1.0,
            alpha: 14.0,
        },
    )
    .unwrap();
    assert!(energy.is_finite() && energy < 0.0);
}

#[test]
fn loads_named_pbe_parameters() {
    assert!(matches!(
        load_named("PBE", 1, false),
        Some(Damping::Rational {
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
            ..
        })
    ));
}

#[test]
fn matches_remaining_damping_references() {
    let cases = [
        (
            Damping::ModifiedZero {
                s6: 1.0,
                s8: 0.0,
                rs6: 2.340218,
                rs8: 1.0,
                alpha: 14.0,
                beta: 0.129434,
            },
            -0.0006540906317314853,
        ),
        (
            Damping::Rational {
                s6: 1.0,
                s8: 0.358940,
                a1: 0.012092,
                a2: 5.938951,
            },
            -0.000674950674541602,
        ),
        (
            Damping::OptimizedPower {
                s6: 0.91826,
                s8: 0.0,
                a1: 0.2,
                a2: 4.75,
                beta: 6.0,
            },
            -0.0005665429155128956,
        ),
        (
            Damping::Cso {
                s6: 1.0,
                a1: 0.24,
                a2: 2.5,
                a3: 0.0,
                a4: 6.25,
            },
            -0.0005716102446037327,
        ),
        (
            Damping::Z {
                s6: 1.0,
                s8: 1.0,
                a1: 1.0,
            },
            -0.0018926574597399033,
        ),
    ];
    for (damping, reference) in cases {
        let actual = energy(&[6, 6], &[0.0, 0.0, 0.0, 6.0, 0.0, 0.0], damping).unwrap();
        assert!(
            (actual - reference).abs() < 1.0e-14,
            "{actual} != {reference}"
        );
    }
}

#[test]
fn matches_gradient_references() {
    let cases = [
        (
            Damping::Zero {
                s6: 1.0,
                s8: 0.722,
                rs6: 1.217,
                rs8: 1.0,
                alpha: 14.0,
            },
            0.00010995853975013035,
            -0.0006597512385007821,
        ),
        (
            Damping::Rational {
                s6: 1.0,
                s8: 0.7875,
                a1: 0.4289,
                a2: 4.4407,
            },
            -0.00018488441198083488,
            0.0011093064718850092,
        ),
        (
            Damping::ModifiedZero {
                s6: 1.0,
                s8: 0.0,
                rs6: 2.340218,
                rs8: 1.0,
                alpha: 14.0,
                beta: 0.129434,
            },
            -0.0004262990601175771,
            0.0025577943607054624,
        ),
        (
            Damping::OptimizedPower {
                s6: 0.91826,
                s8: 0.0,
                a1: 0.2,
                a2: 4.75,
                beta: 6.0,
            },
            -9.919798358517386e-05,
            0.0005951879015110431,
        ),
        (
            Damping::Cso {
                s6: 1.0,
                a1: 0.24,
                a2: 2.5,
                a3: 0.0,
                a4: 6.25,
            },
            -0.00025099222300498694,
            0.0015059533380299215,
        ),
        (
            Damping::Z {
                s6: 1.0,
                s8: 1.0,
                a1: 1.0,
            },
            -0.002172534208460211,
            0.013035205250761267,
        ),
    ];
    for (damping, gradient_reference, virial_reference) in cases {
        let (_, gradient, virial) =
            gradient(&[6, 6], &[0.0, 0.0, 0.0, 6.0, 0.0, 0.0], damping).unwrap();
        assert!((gradient[0] - gradient_reference).abs() < 1.0e-14);
        assert!((gradient[3] + gradient_reference).abs() < 1.0e-14);
        assert!((virial[0] - virial_reference).abs() < 1.0e-13);
    }
}

#[test]
fn matches_analytical_hessian_reference() {
    let (energy, hessian) = hessian(
        &[6, 6],
        &[0.0, 0.0, 0.0, 6.0, 0.0, 0.0],
        Damping::Rational {
            s6: 1.0,
            s8: 0.7875,
            a1: 0.4289,
            a2: 4.4407,
        },
    )
    .unwrap();
    let longitudinal = 4.9042518143140056e-05;
    let transverse = 3.081406866347248e-05;
    let expected = [
        longitudinal,
        0.0,
        0.0,
        -longitudinal,
        0.0,
        0.0,
        0.0,
        transverse,
        0.0,
        0.0,
        -transverse,
        0.0,
        0.0,
        0.0,
        transverse,
        0.0,
        0.0,
        -transverse,
        -longitudinal,
        0.0,
        0.0,
        longitudinal,
        0.0,
        0.0,
        0.0,
        -transverse,
        0.0,
        0.0,
        transverse,
        0.0,
        0.0,
        0.0,
        -transverse,
        0.0,
        0.0,
        transverse,
    ];
    assert!((energy - -0.0005341413931338267).abs() < 1.0e-15);
    for (actual, expected) in hessian.iter().zip(expected) {
        assert!(
            (actual - expected).abs() < 1.0e-13,
            "{actual} != {expected}"
        );
    }
    for row in 0..6 {
        assert!(hessian[row * 6..(row + 1) * 6].iter().sum::<f64>().abs() < 1.0e-15);
    }
}

#[test]
fn realspace_partitions_sum_to_serial() {
    let numbers = [6, 8, 1];
    let positions = [0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 1.0, 3.0, 2.0];
    let damping = Damping::Rational {
        s6: 1.0,
        s8: 0.7875,
        a1: 0.4289,
        a2: 4.4407,
    };
    let serial = gradient(&numbers, &positions, damping).unwrap();
    let serial_hessian = hessian(&numbers, &positions, damping).unwrap().1;
    let parts: Vec<_> = (0..3)
        .map(|part| {
            gradient_partitioned(
                &numbers,
                &positions,
                damping,
                WorkPartition::new(part, 3).unwrap(),
            )
            .unwrap()
        })
        .collect();
    assert!((parts.iter().map(|part| part.0).sum::<f64>() - serial.0).abs() < 1.0e-15);
    for index in 0..positions.len() {
        assert!(
            (parts.iter().map(|part| part.1[index]).sum::<f64>() - serial.1[index]).abs() < 1.0e-14
        );
    }
    for index in 0..9 {
        assert!(
            (parts.iter().map(|part| part.2[index]).sum::<f64>() - serial.2[index]).abs() < 1.0e-14
        );
    }
    for (index, expected) in serial_hessian.iter().enumerate() {
        let sum = (0..3)
            .map(|part| {
                hessian_partitioned(
                    &numbers,
                    &positions,
                    damping,
                    WorkPartition::new(part, 3).unwrap(),
                )
                .unwrap()
                .1[index]
            })
            .sum::<f64>();
        assert!((sum - expected).abs() < 1.0e-13);
    }
}

#[test]
fn matches_atm_reference() {
    let numbers = [6, 8, 7];
    let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 4.0, 0.0];
    let atm = Atm {
        s9: 1.0,
        alpha: 16.0,
    };
    let pairwise = atm_pairwise(&numbers, &positions, atm, WorkPartition::SERIAL).unwrap();
    let scalar = atm_energy_with_cutoff(
        &numbers,
        &positions,
        atm,
        RealspaceCutoff::default(),
        &[],
        WorkPartition::SERIAL,
    )
    .unwrap();
    let (energy, gradient, virial, hessian) =
        atm_derivatives(&numbers, &positions, atm, WorkPartition::SERIAL).unwrap();
    assert!((energy - 1.7772815945555145e-7).abs() < 1.0e-16, "{energy}");
    assert!((scalar - energy).abs() < 1.0e-16);
    assert!((pairwise.iter().sum::<f64>() - energy).abs() < 1.0e-16);
    assert!((pairwise[1] - 2.962135990925684e-8).abs() < 1.0e-17);
    let gradient_reference = [
        -1.2861635008461366e-7,
        -1.1895399868174997e-7,
        0.0,
        1.3131318927472858e-7,
        -3.243533886057104e-8,
        0.0,
        -2.696839190128469e-9,
        1.5138933754231423e-7,
        0.0,
    ];
    for (actual, expected) in gradient.iter().zip(gradient_reference) {
        assert!((actual - expected).abs() < 1.0e-17);
    }
    let first_order = atm_gradient_with_cutoff(
        &numbers,
        &positions,
        atm,
        RealspaceCutoff::default(),
        &[],
        WorkPartition::SERIAL,
    )
    .unwrap();
    assert!((first_order.0 - energy).abs() < 1.0e-16);
    for (actual, expected) in first_order.1.iter().zip(gradient_reference) {
        assert!((actual - expected).abs() < 1.0e-17);
    }
    for (actual, expected) in first_order.2.iter().zip(virial) {
        assert!((actual - expected).abs() < 1.0e-16);
    }
    assert!(
        (virial[0] - 6.538691071835551e-7).abs() < 1.0e-16,
        "{} {}",
        virial[0],
        virial[4]
    );
    assert!((virial[4] - 6.055573501692569e-7).abs() < 1.0e-16);
    assert!((hessian[0] - 6.002548300879269e-8).abs() < 1.0e-17);
    assert!((hessian[1] - 4.646237216554474e-8).abs() < 1.0e-17);
    assert!((hessian[3] - -7.187609143240512e-8).abs() < 1.0e-17);
    let partitioned = (0..2)
        .map(|part| {
            atm_derivatives(
                &numbers,
                &positions,
                atm,
                WorkPartition::new(part, 2).unwrap(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    assert!((partitioned.iter().map(|value| value.0).sum::<f64>() - energy).abs() < 1.0e-16);
    for (index, expected) in gradient.iter().enumerate() {
        assert!(
            (partitioned.iter().map(|value| value.1[index]).sum::<f64>() - expected).abs()
                < 1.0e-16
        );
    }
}
