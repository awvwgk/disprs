//! D3/D3S two-body dispersion, optional ATM three-body terms, and gCP.
//!
//! # Inputs
//! Unless stated otherwise, `numbers` contains `N > 0` atomic numbers in
//! `1..=103`; `positions` contains exactly `3*N` finite Cartesian coordinates
//! in bohr, ordered by atom. Avoid coincident atoms. D3S supports `1..=94`.
//! `damping` supplies finite, physically valid [`Damping`] parameters;
//! `atm` supplies [`Atm`] parameters, independently of two-body damping.
//! `cutoff` supplies nonnegative finite radii and widths in bohr.
//! `ghosts` is empty or has `N` entries: `true` removes that atom from explicit
//! dispersion terms, but not from coordination numbers or their derivatives.
//! `partition` assigns interaction terms to a valid [`WorkPartition`]; sum all
//! partitions to obtain the total. Molecular routines check shapes and elements,
//! but do not uniformly reject nonfinite coordinates or invalid parameters.
//!
//! `lattice` stores three lattice vectors as consecutive triples (column-major
//! 3-by-3 matrix), in bohr. `periodic` selects the active vectors, which must be
//! finite and linearly independent. An optional lattice is required if any
//! periodic flag is true. `cn_cutoff` is a finite, nonnegative CN radius in bohr.
//! Hard cutoff boundaries must not be crossed when interpreting derivatives.
//!
//! # Example
//! ```
//! use disprs::d3::{energy_partitioned, load_named, WorkPartition};
//! let numbers = [8, 1, 1];
//! let positions = [0.0, 0.0, 0.0, 1.8, 0.0, 0.0, -0.45, 1.74, 0.0];
//! let damping = load_named("pbe", 1, false).unwrap();
//! let energy = energy_partitioned(&numbers, &positions, damping, WorkPartition::SERIAL)?;
//! assert!(energy.is_finite());
//! # Ok::<(), &'static str>(())
//! ```

use std::convert::TryInto;

use crate::geometry::{displacement, squared_distance};
use crate::parameters::{lookup, parameter};
pub use crate::WorkPartition;

const ELEMENTS: usize = 103;
const REFERENCES: usize = 7;
const PAIRS: usize = ELEMENTS * (ELEMENTS + 1) / 2;
const COVALENT_RADII: &[u8; ELEMENTS * 8] = include_bytes!("../assets/covalent_radii.bin");
const R4R2: &[u8; ELEMENTS * 8] = include_bytes!("../assets/r4r2.bin");
const VDW_RADII: &[u8; PAIRS * 8] = include_bytes!("../assets/vdw_radii.bin");
const REFERENCE_CN: &[u8; ELEMENTS * REFERENCES * 8] = include_bytes!("../assets/reference_cn.bin");
const REFERENCE_C6: &[u8; PAIRS * REFERENCES * REFERENCES * 8] =
    include_bytes!("../assets/reference_c6.bin");
pub(crate) mod fourier;
mod gcp;
mod realspace;
mod smooth;
pub use self::smooth::Model;
#[cfg(test)]
mod tests;
pub use self::fourier::{
    periodic_derivatives_partitioned, periodic_derivatives_partitioned_with_cutoff,
    periodic_derivatives_partitioned_with_ghosts, periodic_energy_partitioned,
    periodic_energy_partitioned_with_cutoff, periodic_energy_partitioned_with_ghosts, EwaldConfig,
    PeriodicResult,
};
pub use self::gcp::{
    derivatives as gcp_derivatives, energy as gcp_energy, hessian as gcp_hessian, load as load_gcp,
    Gcp, GcpCutoff, GcpResult,
};
pub use self::realspace::{periodic_realspace, periodic_realspace_hessian};
type Prepared = (Vec<usize>, Vec<usize>, Vec<f64>, Vec<[f64; REFERENCES]>);
type AtmDerivatives = (f64, Vec<f64>, [f64; 9], Vec<f64>);
type AtmGradient = (f64, Vec<f64>, [f64; 9]);

/// Axilrod-Teller-Muto three-body parameters; values are caller-validated.
#[derive(Clone, Copy)]
pub struct Atm {
    /// Finite dimensionless three-body scale; zero disables the contribution.
    pub s9: f64,
    /// Finite positive damping exponent, usually 16.
    pub alpha: f64,
}

/// Real-space summation ranges in bohr; all fields must be finite and nonnegative.
///
/// Defaults: CN 40, pairs 60, triplets 40, no switching. For smooth truncation
/// use `0 < width < radius`; D3 treats `width >= radius` as no switching.
#[derive(Clone, Copy)]
pub struct RealspaceCutoff {
    /// Coordination-number radius.
    pub cn: f64,
    /// Two-body interaction radius.
    pub disp2: f64,
    /// Three-body edge radius.
    pub disp3: f64,
    /// Two-body quintic switching interval; zero gives a hard cutoff.
    pub width2: f64,
    /// Three-body edge switching interval; zero gives a hard cutoff.
    pub width3: f64,
}

impl Default for RealspaceCutoff {
    fn default() -> Self {
        Self {
            cn: 40.0,
            disp2: 60.0,
            disp3: 40.0,
            width2: 0.0,
            width3: 0.0,
        }
    }
}

/// Two-body damping family and parameters, preferably obtained with [`load_named`].
///
/// All fields must be finite. Scales may be signed to represent fitted values;
/// radius expressions and power bases must remain positive. These preconditions
/// are not uniformly checked by evaluators. Length parameters use bohr.
#[derive(Clone, Copy)]
pub enum Damping {
    /// Original zero damping with positive `rs6`, `rs8`, and `alpha`.
    Zero {
        /// Dimensionless C6 scale.
        s6: f64,
        /// Dimensionless C8 scale.
        s8: f64,
        /// Dimensionless C6 reference-radius multiplier.
        rs6: f64,
        /// Dimensionless C8 reference-radius multiplier.
        rs8: f64,
        /// C6 damping exponent; C8 uses `alpha + 2`.
        alpha: f64,
    },
    /// Becke-Johnson damping; `a1*sqrt(C8/C6) + a2` must be positive.
    Rational {
        /// Dimensionless C6 scale.
        s6: f64,
        /// Dimensionless C8 scale.
        s8: f64,
        /// Dimensionless radius slope.
        a1: f64,
        /// Radius offset in bohr.
        a2: f64,
    },
    /// Modified zero damping; positive radii/exponent and positive shifted bases.
    ModifiedZero {
        /// Dimensionless C6 scale.
        s6: f64,
        /// Dimensionless C8 scale.
        s8: f64,
        /// Dimensionless C6 reference-radius multiplier; positive.
        rs6: f64,
        /// Dimensionless C8 reference-radius multiplier; positive.
        rs8: f64,
        /// C6 damping exponent; positive, with C8 using `alpha + 2`.
        alpha: f64,
        /// Inverse-bohr shift coefficient in `R/(rs*RvdW) + beta*RvdW`.
        beta: f64,
    },
    /// Optimized-power damping; positive BJ radius and normally nonnegative `beta`.
    OptimizedPower {
        /// Dimensionless C6 scale.
        s6: f64,
        /// Dimensionless C8 scale.
        s8: f64,
        /// Dimensionless BJ radius slope.
        a1: f64,
        /// BJ radius offset in bohr.
        a2: f64,
        /// Additional dimensionless damping power.
        beta: f64,
    },
    /// C6-only sigmoid scaling with a positive `a3*sqrt(C8/C6) + a4` radius.
    Cso {
        /// Long-range dimensionless C6 scale.
        s6: f64,
        /// Dimensionless sigmoid amplitude.
        a1: f64,
        /// Dimensionless sigmoid midpoint multiplier of `sqrt(C8/C6)`.
        a2: f64,
        /// Dimensionless damping-radius slope.
        a3: f64,
        /// Damping-radius offset in bohr.
        a4: f64,
    },
    /// C6-dependent Z damping; not available from the named-parameter table.
    Z {
        /// Dimensionless C6 scale.
        s6: f64,
        /// Dimensionless C8 scale.
        s8: f64,
        /// Positive coefficient in the C6-dependent denominator (atomic units).
        a1: f64,
    },
}

/// Load fitted two-body parameters for functional `method`.
///
/// `damping`: 0 zero, 1 BJ, 2 modified zero, 3 modified BJ, 4 optimized power,
/// 5 CSO. Unknown names/selectors return `None`. `_atm` is ignored; add [`Atm`]
/// separately. Method spelling is normalized by the parameter lookup.
pub fn load_named(method: &str, damping: i32, _atm: bool) -> Option<Damping> {
    if !(0..=5).contains(&damping) {
        return None;
    }
    let line = lookup(
        method,
        match damping {
            0 => "d3.zero",
            1 => "d3.bj",
            2 => "d3.zerom",
            3 => "d3.bjm",
            4 => "d3.op",
            5 => "d3.cso",
            _ => unreachable!(),
        },
    )?;
    Some(match damping {
        0 => Damping::Zero {
            s6: parameter(line, "s6", 1.0),
            s8: parameter(line, "s8", 0.0),
            rs6: parameter(line, "rs6", 1.0),
            rs8: parameter(line, "rs8", 1.0),
            alpha: parameter(line, "alp", 14.0),
        },
        1 | 3 => Damping::Rational {
            s6: parameter(line, "s6", 1.0),
            s8: parameter(line, "s8", 0.0),
            a1: parameter(line, "a1", 0.0),
            a2: parameter(line, "a2", 0.0),
        },
        2 => Damping::ModifiedZero {
            s6: parameter(line, "s6", 1.0),
            s8: parameter(line, "s8", 0.0),
            rs6: parameter(line, "rs6", 1.0),
            rs8: parameter(line, "rs8", 1.0),
            alpha: parameter(line, "alp", 14.0),
            beta: parameter(line, "bet", 0.0),
        },
        4 => Damping::OptimizedPower {
            s6: parameter(line, "s6", 0.0),
            s8: parameter(line, "s8", 0.0),
            a1: parameter(line, "a1", 0.0),
            a2: parameter(line, "a2", 0.0),
            beta: parameter(line, "bet", 0.0),
        },
        5 => Damping::Cso {
            s6: parameter(line, "s6", 1.0),
            a1: parameter(line, "a1", 0.0),
            a2: parameter(line, "a2", 2.5),
            a3: parameter(line, "rs6", 0.0),
            a4: parameter(line, "rs8", 6.25),
        },
        _ => unreachable!(),
    })
}

pub(super) fn value(data: &[u8], index: usize) -> f64 {
    f64::from_le_bytes(data[index * 8..(index + 1) * 8].try_into().unwrap())
}

fn pair_index(first: usize, second: usize) -> usize {
    let (low, high) = if first <= second {
        (first, second)
    } else {
        (second, first)
    };
    low + high * (high + 1) / 2
}

pub(super) fn reference_count(element: usize) -> usize {
    (0..REFERENCES)
        .take_while(|&reference| value(REFERENCE_CN, element * REFERENCES + reference) >= 0.0)
        .count()
}

pub(super) fn weights(element: usize, coordination: f64) -> [f64; REFERENCES] {
    let count = reference_count(element);
    let mut result = [0.0; REFERENCES];
    let mut norm = 0.0;
    for (reference, weight) in result.iter_mut().enumerate().take(count) {
        let delta = coordination - value(REFERENCE_CN, element * REFERENCES + reference);
        *weight = (-4.0 * delta * delta).exp();
        norm += *weight;
    }
    if norm.is_finite() && norm > 0.0 {
        for weight in &mut result[..count] {
            *weight /= norm;
        }
    } else {
        result[count - 1] = 1.0;
    }
    result
}

fn weight_derivatives(
    element: usize,
    coordination: f64,
    weights: &[f64; REFERENCES],
) -> [f64; REFERENCES] {
    let count = reference_count(element);
    let mean = weights[..count]
        .iter()
        .enumerate()
        .map(|(reference, weight)| {
            weight * 8.0 * (value(REFERENCE_CN, element * REFERENCES + reference) - coordination)
        })
        .sum::<f64>();
    let mut derivatives = [0.0; REFERENCES];
    for reference in 0..count {
        derivatives[reference] = weights[reference]
            * (8.0 * (value(REFERENCE_CN, element * REFERENCES + reference) - coordination) - mean);
    }
    derivatives
}

fn weight_second_derivatives(
    element: usize,
    coordination: f64,
    weights: &[f64; REFERENCES],
) -> [f64; REFERENCES] {
    let count = reference_count(element);
    let slopes: [f64; REFERENCES] = std::array::from_fn(|reference| {
        8.0 * (value(REFERENCE_CN, element * REFERENCES + reference) - coordination)
    });
    let mean = weights[..count]
        .iter()
        .enumerate()
        .map(|(reference, weight)| weight * slopes[reference])
        .sum::<f64>();
    let variance = weights[..count]
        .iter()
        .enumerate()
        .map(|(reference, weight)| weight * (slopes[reference] - mean).powi(2))
        .sum::<f64>();
    let mut result = [0.0; REFERENCES];
    for reference in 0..count {
        result[reference] = weights[reference] * ((slopes[reference] - mean).powi(2) - variance);
    }
    result
}

pub(super) fn reference_c6(
    first_element: usize,
    second_element: usize,
    first_reference: usize,
    second_reference: usize,
) -> f64 {
    let pair = pair_index(first_element, second_element);
    let index = if first_element > second_element {
        first_reference + REFERENCES * second_reference
    } else {
        second_reference + REFERENCES * first_reference
    };
    value(REFERENCE_C6, pair * REFERENCES * REFERENCES + index)
}

fn atomic_c6(
    first_element: usize,
    second_element: usize,
    first_weights: &[f64; REFERENCES],
    second_weights: &[f64; REFERENCES],
) -> f64 {
    let mut result = 0.0;
    for (first_reference, &first_weight) in first_weights
        .iter()
        .enumerate()
        .take(reference_count(first_element))
    {
        for (second_reference, &second_weight) in second_weights
            .iter()
            .enumerate()
            .take(reference_count(second_element))
        {
            result += first_weight
                * second_weight
                * reference_c6(
                    first_element,
                    second_element,
                    first_reference,
                    second_reference,
                );
        }
    }
    result
}

fn atomic_c6_derivatives(
    first_element: usize,
    second_element: usize,
    first_weights: &[f64; REFERENCES],
    second_weights: &[f64; REFERENCES],
    first_derivatives: &[f64; REFERENCES],
    second_derivatives: &[f64; REFERENCES],
) -> (f64, f64, f64) {
    let mut c6 = 0.0;
    let mut first_derivative = 0.0;
    let mut second_derivative = 0.0;
    for first_reference in 0..reference_count(first_element) {
        for second_reference in 0..reference_count(second_element) {
            let reference = reference_c6(
                first_element,
                second_element,
                first_reference,
                second_reference,
            );
            c6 += first_weights[first_reference] * second_weights[second_reference] * reference;
            first_derivative +=
                first_derivatives[first_reference] * second_weights[second_reference] * reference;
            second_derivative +=
                first_weights[first_reference] * second_derivatives[second_reference] * reference;
        }
    }
    (c6, first_derivative, second_derivative)
}

#[allow(clippy::too_many_arguments)]
fn atomic_c6_second_derivatives(
    first_element: usize,
    second_element: usize,
    first_weights: &[f64; REFERENCES],
    second_weights: &[f64; REFERENCES],
    first_derivatives: &[f64; REFERENCES],
    second_derivatives: &[f64; REFERENCES],
    first_second_derivatives: &[f64; REFERENCES],
    second_second_derivatives: &[f64; REFERENCES],
) -> [f64; 6] {
    let mut result = [0.0; 6];
    for first_reference in 0..reference_count(first_element) {
        for second_reference in 0..reference_count(second_element) {
            let reference = reference_c6(
                first_element,
                second_element,
                first_reference,
                second_reference,
            );
            let first_weight = first_weights[first_reference];
            let second_weight = second_weights[second_reference];
            result[0] += first_weight * second_weight * reference;
            result[1] += first_derivatives[first_reference] * second_weight * reference;
            result[2] += first_weight * second_derivatives[second_reference] * reference;
            result[3] += first_second_derivatives[first_reference] * second_weight * reference;
            result[4] += first_derivatives[first_reference]
                * second_derivatives[second_reference]
                * reference;
            result[5] += first_weight * second_second_derivatives[second_reference] * reference;
        }
    }
    result
}

fn inverse_sum_power(u: f64, exponent: i32, constant: f64) -> [f64; 3] {
    let value = 1.0 / (u.powi(exponent) + constant);
    let first = -exponent as f64 * u.powi(exponent - 1) * value * value;
    let second = -(exponent * (exponent - 1)) as f64 * u.powi(exponent - 2) * value * value
        + 2.0 * (exponent * exponent) as f64 * u.powi(2 * exponent - 2) * value.powi(3);
    [value, first, second]
}

fn damped_inverse_power(u: f64, radius: f64, exponent: f64, power: i32) -> [f64; 3] {
    let half_exponent = 0.5 * exponent;
    let ratio = (radius / u.sqrt()).powf(exponent);
    let ratio_first = -half_exponent * ratio / u;
    let ratio_second = half_exponent * (half_exponent + 1.0) * ratio / (u * u);
    let damping = 1.0 / (1.0 + 6.0 * ratio);
    let damping_first = -6.0 * ratio_first * damping * damping;
    let damping_second = -6.0 * ratio_second * damping * damping
        + 72.0 * ratio_first * ratio_first * damping.powi(3);
    let inverse = u.powi(-power);
    let inverse_first = -power as f64 * inverse / u;
    let inverse_second = (power * (power + 1)) as f64 * inverse / (u * u);
    [
        damping * inverse,
        damping_first * inverse + damping * inverse_first,
        damping_second * inverse + 2.0 * damping_first * inverse_first + damping * inverse_second,
    ]
}

fn kernel_second_derivatives(
    damping: Damping,
    first_element: usize,
    second_element: usize,
    first_species: usize,
    second_species: usize,
    u: f64,
    c6: f64,
) -> [f64; 6] {
    let rrij = 3.0 * value(R4R2, first_element) * value(R4R2, second_element);
    if let Damping::Z { s6, s8, a1 } = damping {
        let radius = a1 / (first_species + second_species + 2) as f64;
        let b8 = radius * rrij;
        let t6 = 1.0 / (u.powi(3) + radius * c6);
        let t8 = 1.0 / (u.powi(4) + b8 * c6);
        let t6u = -3.0 * u.powi(2) * t6.powi(2);
        let t8u = -4.0 * u.powi(3) * t8.powi(2);
        let t6uu = -6.0 * u * t6.powi(2) + 18.0 * u.powi(4) * t6.powi(3);
        let t8uu = -12.0 * u.powi(2) * t8.powi(2) + 32.0 * u.powi(6) * t8.powi(3);
        let t6c = -radius * t6.powi(2);
        let t8c = -b8 * t8.powi(2);
        let t6uc = 6.0 * u.powi(2) * radius * t6.powi(3);
        let t8uc = 8.0 * u.powi(3) * b8 * t8.powi(3);
        let t6cc = 2.0 * radius * radius * t6.powi(3);
        let t8cc = 2.0 * b8 * b8 * t8.powi(3);
        let p = s6 * t6 + s8 * rrij * t8;
        let pu = s6 * t6u + s8 * rrij * t8u;
        let puu = s6 * t6uu + s8 * rrij * t8uu;
        let pc = s6 * t6c + s8 * rrij * t8c;
        let puc = s6 * t6uc + s8 * rrij * t8uc;
        let pcc = s6 * t6cc + s8 * rrij * t8cc;
        return [
            -c6 * p,
            -c6 * pu,
            -p - c6 * pc,
            -c6 * puu,
            -pu - c6 * puc,
            -2.0 * pc - c6 * pcc,
        ];
    }
    let (potential, first, second) = match damping {
        Damping::Rational { s6, s8, a1, a2 } => {
            let radius = a1 * rrij.sqrt() + a2;
            let sixth = inverse_sum_power(u, 3, radius.powi(6));
            let eighth = inverse_sum_power(u, 4, radius.powi(8));
            (
                s6 * sixth[0] + s8 * rrij * eighth[0],
                s6 * sixth[1] + s8 * rrij * eighth[1],
                s6 * sixth[2] + s8 * rrij * eighth[2],
            )
        }
        Damping::Zero {
            s6,
            s8,
            rs6,
            rs8,
            alpha,
        } => {
            let radius = value(VDW_RADII, pair_index(first_element, second_element));
            let sixth = damped_inverse_power(u, rs6 * radius, alpha, 3);
            let eighth = damped_inverse_power(u, rs8 * radius, alpha + 2.0, 4);
            (
                s6 * sixth[0] + s8 * rrij * eighth[0],
                s6 * sixth[1] + s8 * rrij * eighth[1],
                s6 * sixth[2] + s8 * rrij * eighth[2],
            )
        }
        Damping::ModifiedZero {
            s6,
            s8,
            rs6,
            rs8,
            alpha,
            beta,
        } => {
            let radius = value(VDW_RADII, pair_index(first_element, second_element));
            let term = |scale: f64, exponent: f64, power: i32| {
                let root = u.sqrt();
                let w = root / (scale * radius) + beta * radius;
                let w_first = 1.0 / (2.0 * scale * radius * root);
                let w_second = -1.0 / (4.0 * scale * radius * u * root);
                let ratio = w.powf(-exponent);
                let ratio_first = -exponent * w.powf(-exponent - 1.0) * w_first;
                let ratio_second =
                    exponent * (exponent + 1.0) * w.powf(-exponent - 2.0) * w_first * w_first
                        - exponent * w.powf(-exponent - 1.0) * w_second;
                let damping = 1.0 / (1.0 + 6.0 * ratio);
                let damping_first = -6.0 * ratio_first * damping * damping;
                let damping_second = -6.0 * ratio_second * damping * damping
                    + 72.0 * ratio_first * ratio_first * damping.powi(3);
                let inverse = u.powi(-power);
                let inverse_first = -power as f64 * inverse / u;
                let inverse_second = (power * (power + 1)) as f64 * inverse / (u * u);
                [
                    damping * inverse,
                    damping_first * inverse + damping * inverse_first,
                    damping_second * inverse
                        + 2.0 * damping_first * inverse_first
                        + damping * inverse_second,
                ]
            };
            let sixth = term(rs6, alpha, 3);
            let eighth = term(rs8, alpha + 2.0, 4);
            (
                s6 * sixth[0] + s8 * rrij * eighth[0],
                s6 * sixth[1] + s8 * rrij * eighth[1],
                s6 * sixth[2] + s8 * rrij * eighth[2],
            )
        }
        Damping::OptimizedPower {
            s6,
            s8,
            a1,
            a2,
            beta,
        } => {
            let radius = a1 * rrij.sqrt() + a2;
            let q = 0.5 * beta;
            let term = |power: i32| {
                let constant = radius.powf(2.0 * power as f64 + beta);
                let denominator = u.powi(power) + constant * u.powf(-q);
                let first = power as f64 * u.powi(power - 1) - q * constant * u.powf(-q - 1.0);
                let second = (power * (power - 1)) as f64 * u.powi(power - 2)
                    + q * (q + 1.0) * constant * u.powf(-q - 2.0);
                let inverse = denominator.recip();
                [
                    inverse,
                    -first * inverse * inverse,
                    -second * inverse * inverse + 2.0 * first * first * inverse.powi(3),
                ]
            };
            let sixth = term(3);
            let eighth = term(4);
            (
                s6 * sixth[0] + s8 * rrij * eighth[0],
                s6 * sixth[1] + s8 * rrij * eighth[1],
                s6 * sixth[2] + s8 * rrij * eighth[2],
            )
        }
        Damping::Cso { s6, a1, a2, a3, a4 } => {
            let root = u.sqrt();
            let radius = rrij.sqrt();
            let sigmoid = 1.0 / (1.0 + (root - a2 * radius).exp());
            let scale = s6 + a1 * sigmoid;
            let scale_r = -a1 * sigmoid * (1.0 - sigmoid);
            let scale_rr = a1 * sigmoid * (1.0 - sigmoid) * (1.0 - 2.0 * sigmoid);
            let scale_u = scale_r / (2.0 * root);
            let scale_uu = scale_rr / (4.0 * u) - scale_r / (4.0 * u * root);
            let inverse = inverse_sum_power(u, 3, (a3 * radius + a4).powi(6));
            (
                scale * inverse[0],
                scale_u * inverse[0] + scale * inverse[1],
                scale_uu * inverse[0] + 2.0 * scale_u * inverse[1] + scale * inverse[2],
            )
        }
        Damping::Z { .. } => unreachable!(),
    };
    [
        -c6 * potential,
        -c6 * first,
        -potential,
        -c6 * second,
        -first,
        0.0,
    ]
}

fn kernel(
    damping: Damping,
    first_element: usize,
    second_element: usize,
    first_species: usize,
    second_species: usize,
    u: f64,
    c6: f64,
) -> (f64, f64, f64) {
    let rrij = 3.0 * value(R4R2, first_element) * value(R4R2, second_element);
    let (p, derivative, c6_derivative) = match damping {
        Damping::Rational { s6, s8, a1, a2 } => {
            let r0 = a1 * rrij.sqrt() + a2;
            let t6 = 1.0 / (u.powi(3) + r0.powi(6));
            let t8 = 1.0 / (u.powi(4) + r0.powi(8));
            (
                s6 * t6 + s8 * rrij * t8,
                s6 * -3.0 * u.powi(2) * t6.powi(2) + s8 * rrij * -4.0 * u.powi(3) * t8.powi(2),
                0.0,
            )
        }
        Damping::Zero {
            s6,
            s8,
            rs6,
            rs8,
            alpha,
        } => {
            let r0 = value(VDW_RADII, pair_index(first_element, second_element));
            let t6 = (rs6 * r0 / u.sqrt()).powf(alpha);
            let t8 = (rs8 * r0 / u.sqrt()).powf(alpha + 2.0);
            let f6 = 1.0 / (1.0 + 6.0 * t6);
            let f8 = 1.0 / (1.0 + 6.0 * t8);
            let f6d = 3.0 * alpha * t6 * f6 * f6 / u;
            let f8d = 3.0 * (alpha + 2.0) * t8 * f8 * f8 / u;
            (
                s6 * f6 / u.powi(3) + s8 * rrij * f8 / u.powi(4),
                s6 * (f6d / u.powi(3) - 3.0 * f6 / u.powi(4))
                    + s8 * rrij * (f8d / u.powi(4) - 4.0 * f8 / u.powi(5)),
                0.0,
            )
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
            let r = u.sqrt();
            let w6 = r / (rs6 * r0) + beta * r0;
            let w8 = r / (rs8 * r0) + beta * r0;
            let t6 = w6.powf(-alpha);
            let t8 = w8.powf(-alpha - 2.0);
            let t6d = -alpha * w6.powf(-alpha - 1.0) / (2.0 * rs6 * r0 * r);
            let t8d = -(alpha + 2.0) * w8.powf(-alpha - 3.0) / (2.0 * rs8 * r0 * r);
            let f6 = 1.0 / (1.0 + 6.0 * t6);
            let f8 = 1.0 / (1.0 + 6.0 * t8);
            let f6d = -6.0 * t6d * f6 * f6;
            let f8d = -6.0 * t8d * f8 * f8;
            (
                s6 * f6 / u.powi(3) + s8 * rrij * f8 / u.powi(4),
                s6 * (f6d / u.powi(3) - 3.0 * f6 / u.powi(4))
                    + s8 * rrij * (f8d / u.powi(4) - 4.0 * f8 / u.powi(5)),
                0.0,
            )
        }
        Damping::OptimizedPower {
            s6,
            s8,
            a1,
            a2,
            beta,
        } => {
            let r0 = a1 * rrij.sqrt() + a2;
            let q = 0.5 * beta;
            let d6 = u.powi(3) + r0.powf(6.0 + beta) * u.powf(-q);
            let d8 = u.powi(4) + r0.powf(8.0 + beta) * u.powf(-q);
            let d6d = 3.0 * u.powi(2) - q * r0.powf(6.0 + beta) * u.powf(-q - 1.0);
            let d8d = 4.0 * u.powi(3) - q * r0.powf(8.0 + beta) * u.powf(-q - 1.0);
            (
                s6 / d6 + s8 * rrij / d8,
                -s6 * d6d / d6.powi(2) - s8 * rrij * d8d / d8.powi(2),
                0.0,
            )
        }
        Damping::Cso { s6, a1, a2, a3, a4 } => {
            let r = u.sqrt();
            let r0 = rrij.sqrt();
            let sigmoid = 1.0 / (1.0 + (r - a2 * r0).exp());
            let scale = s6 + a1 * sigmoid;
            let scale_derivative = -a1 * sigmoid * (1.0 - sigmoid) / (2.0 * r);
            let t6 = 1.0 / (u.powi(3) + (a3 * r0 + a4).powi(6));
            (
                scale * t6,
                scale_derivative * t6 - 3.0 * scale * u.powi(2) * t6.powi(2),
                0.0,
            )
        }
        Damping::Z { s6, s8, a1 } => {
            let r0 = a1 / (first_species + second_species + 2) as f64;
            let t6 = 1.0 / (u.powi(3) + r0 * c6);
            let t8 = 1.0 / (u.powi(4) + r0 * rrij * c6);
            let p = s6 * t6 + s8 * rrij * t8;
            let pu = -3.0 * s6 * u.powi(2) * t6.powi(2) - 4.0 * s8 * rrij * u.powi(3) * t8.powi(2);
            let pc = -s6 * r0 * t6.powi(2) - s8 * rrij.powi(2) * r0 * t8.powi(2);
            (p, pu, pc)
        }
    };
    (-c6 * p, -c6 * derivative, -p - c6 * c6_derivative)
}

/// Molecular two-body energy in hartree using default cutoffs and no ghosts.
/// All arguments follow the [module input contract](#inputs); excludes ATM.
pub fn energy_partitioned(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    energy_partitioned_with_ghosts(numbers, positions, damping, &[], partition)
}

/// Two-body energy with a ghost mask; uses default cutoffs and excludes ATM.
/// All arguments follow the [module input contract](#inputs).
pub fn energy_partitioned_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    energy_partitioned_with_cutoff(
        numbers,
        positions,
        damping,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

/// Two-body energy with explicit cutoffs, ghost mask, and work partition.
/// All arguments follow the [module input contract](#inputs); excludes ATM.
pub fn energy_partitioned_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    Ok(pair_energies::<false>(numbers, positions, damping, cutoff, ghosts, partition)?.0)
}

/// Return two-body `(energy, gradient, virial)` with default cutoffs and no ghosts.
/// Gradient has `3*N` entries (hartree/bohr); virial uses `row + 3*column`
/// (hartree). All inputs follow the [module contract](#inputs); excludes ATM.
pub fn gradient_partitioned(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>, [f64; 9]), &'static str> {
    gradient_partitioned_with_ghosts(numbers, positions, damping, &[], partition)
}

/// As [`gradient_partitioned`], with an empty or `N`-entry ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn gradient_partitioned_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>, [f64; 9]), &'static str> {
    gradient_partitioned_with_cutoff(
        numbers,
        positions,
        damping,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

/// As [`gradient_partitioned`], with explicit cutoffs and ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn gradient_partitioned_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>, [f64; 9]), &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, species, coordination, weights) =
        prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let derivatives: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .zip(weights.iter())
        .map(|((&element, &coordination), weights)| {
            weight_derivatives(element, coordination, weights)
        })
        .collect();
    let contributions = crate::parallel::map(numbers.len(), 256, |start, stride| {
        let mut energy = 0.0;
        let mut gradient = vec![0.0; positions.len()];
        let mut virial = [0.0; 9];
        let mut dedcn = vec![0.0; numbers.len()];

        for first in (start..numbers.len()).step_by(stride) {
            for second in 0..first {
                if ghosts.get(first) == Some(&true)
                    || ghosts.get(second) == Some(&true)
                    || !partition.owns_pair(first, second)
                {
                    continue;
                }
                let vector = displacement(positions, first, second);
                let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
                if distance2 > cutoff.disp2 * cutoff.disp2 || distance2 <= f64::EPSILON {
                    continue;
                }
                let (c6, dc6_first, dc6_second) = atomic_c6_derivatives(
                    elements[first],
                    elements[second],
                    &weights[first],
                    &weights[second],
                    &derivatives[first],
                    &derivatives[second],
                );
                let (pair_energy, dedr2, dedc6) = kernel(
                    damping,
                    elements[first],
                    elements[second],
                    species[first],
                    species[second],
                    distance2,
                    c6,
                );
                let (switch, switch_derivative) =
                    smooth_cutoff(distance2, cutoff.disp2, cutoff.width2);
                let switched_derivative = switch * dedr2 + switch_derivative * pair_energy;
                energy += switch * pair_energy;
                dedcn[first] += switch * dedc6 * dc6_first;
                dedcn[second] += switch * dedc6 * dc6_second;
                for axis in 0..3 {
                    let component = 2.0 * switched_derivative * vector[axis];
                    gradient[3 * first + axis] += component;
                    gradient[3 * second + axis] -= component;
                    for other in 0..3 {
                        virial[axis + 3 * other] += component * vector[other];
                    }
                }
            }
        }

        (energy, gradient, virial, dedcn)
    });
    let mut contributions = contributions.into_iter();
    let (mut energy, mut gradient, mut virial, mut dedcn) = contributions.next().unwrap();
    for (local_energy, local_gradient, local_virial, local_dedcn) in contributions {
        energy += local_energy;
        for (total, value) in gradient.iter_mut().zip(local_gradient) {
            *total += value;
        }
        for (total, value) in virial.iter_mut().zip(local_virial) {
            *total += value;
        }
        for (total, value) in dedcn.iter_mut().zip(local_dedcn) {
            *total += value;
        }
    }
    let contributions = crate::parallel::map(numbers.len(), 256, |start, stride| {
        let mut gradient = vec![0.0; positions.len()];
        let mut virial = [0.0; 9];
        for first in (start..numbers.len()).step_by(stride) {
            for second in 0..first {
                let vector = displacement(positions, first, second);
                let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
                if distance2 > cutoff.cn * cutoff.cn || distance2 < 1.0e-12 {
                    continue;
                }
                let distance = distance2.sqrt();
                let radius = value(COVALENT_RADII, elements[first])
                    + value(COVALENT_RADII, elements[second]);
                let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
                let derivative = -count * (1.0 - count) * 16.0 * radius / distance2;
                for axis in 0..3 {
                    let component =
                        derivative * vector[axis] / distance * (dedcn[first] + dedcn[second]);
                    gradient[3 * first + axis] += component;
                    gradient[3 * second + axis] -= component;
                    for other in 0..3 {
                        virial[axis + 3 * other] += component * vector[other];
                    }
                }
            }
        }
        (gradient, virial)
    });
    for (local_gradient, local_virial) in contributions {
        for (total, value) in gradient.iter_mut().zip(local_gradient) {
            *total += value;
        }
        for (total, value) in virial.iter_mut().zip(local_virial) {
            *total += value;
        }
    }
    Ok((energy, gradient, virial))
}

/// Return two-body energy and row-major `(3*N, 3*N)` Cartesian Hessian.
/// Uses default cutoffs and no ghosts; excludes ATM. Hessian units are
/// hartree/bohr squared. Inputs follow the [module contract](#inputs).
pub fn hessian_partitioned(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>), &'static str> {
    hessian_partitioned_with_ghosts(numbers, positions, damping, &[], partition)
}

/// As [`hessian_partitioned`], with an empty or `N`-entry ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn hessian_partitioned_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>), &'static str> {
    hessian_partitioned_with_cutoff(
        numbers,
        positions,
        damping,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

/// As [`hessian_partitioned`], with explicit cutoffs and ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn hessian_partitioned_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>), &'static str> {
    analytical_hessian(numbers, positions, damping, cutoff, ghosts, partition)
}

#[allow(clippy::needless_range_loop)]
fn analytical_hessian(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>), &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, species, coordination, weights) =
        prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let atoms = numbers.len();
    let size = positions.len();
    let first_weights: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .zip(weights.iter())
        .map(|((&element, &coordination), weights)| {
            weight_derivatives(element, coordination, weights)
        })
        .collect();
    let second_weights: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .zip(weights.iter())
        .map(|((&element, &coordination), weights)| {
            weight_second_derivatives(element, coordination, weights)
        })
        .collect();
    let mut cn_jacobian = vec![0.0; size * atoms];
    for first in 0..atoms {
        for second in 0..first {
            let vector = displacement(positions, first, second);
            let u = vector.iter().map(|value| value * value).sum::<f64>();
            if u > cutoff.cn * cutoff.cn || u < 1.0e-12 {
                continue;
            }
            let distance = u.sqrt();
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
            let radial = -count * (1.0 - count) * 16.0 * radius / u;
            for (axis, &coordinate) in vector.iter().enumerate() {
                let derivative = radial * coordinate / distance;
                let first_coordinate = 3 * first + axis;
                let second_coordinate = 3 * second + axis;
                cn_jacobian[first_coordinate * atoms + first] += derivative;
                cn_jacobian[second_coordinate * atoms + first] -= derivative;
                cn_jacobian[first_coordinate * atoms + second] += derivative;
                cn_jacobian[second_coordinate * atoms + second] -= derivative;
            }
        }
    }
    let mut energy = 0.0;
    let mut hessian = vec![0.0; size * size];
    let mut dedcn = vec![0.0; atoms];
    let mut mixed = vec![0.0; size * atoms];
    let mut cn_hessian = vec![0.0; atoms * atoms];
    for first in 0..atoms {
        for second in 0..first {
            if ghosts.get(first) == Some(&true)
                || ghosts.get(second) == Some(&true)
                || !partition.owns_pair(first, second)
            {
                continue;
            }
            let vector = displacement(positions, first, second);
            let u = vector.iter().map(|value| value * value).sum::<f64>();
            if u > cutoff.disp2 * cutoff.disp2 || u <= f64::EPSILON {
                continue;
            }
            let c6 = atomic_c6_second_derivatives(
                elements[first],
                elements[second],
                &weights[first],
                &weights[second],
                &first_weights[first],
                &first_weights[second],
                &second_weights[first],
                &second_weights[second],
            );
            let kernel = kernel_second_derivatives(
                damping,
                elements[first],
                elements[second],
                species[first],
                species[second],
                u,
                c6[0],
            );
            let switch = smooth_cutoff_second(u, cutoff.disp2, cutoff.width2);
            let pair = [
                switch[0] * kernel[0],
                switch[1] * kernel[0] + switch[0] * kernel[1],
                switch[0] * kernel[2],
                switch[2] * kernel[0] + 2.0 * switch[1] * kernel[1] + switch[0] * kernel[3],
                switch[1] * kernel[2] + switch[0] * kernel[4],
                switch[0] * kernel[5],
            ];
            energy += pair[0];
            let direct_radial = 2.0 * pair[1];
            let direct_cartesian = 4.0 * pair[3];
            for row in 0..3 {
                for column in 0..3 {
                    let block = direct_cartesian * vector[row] * vector[column]
                        + if row == column { direct_radial } else { 0.0 };
                    add_hessian_block(&mut hessian, size, first, second, row, column, block);
                }
            }
            dedcn[first] += pair[2] * c6[1];
            dedcn[second] += pair[2] * c6[2];
            cn_hessian[first * atoms + first] += pair[5] * c6[1] * c6[1] + pair[2] * c6[3];
            cn_hessian[second * atoms + second] += pair[5] * c6[2] * c6[2] + pair[2] * c6[5];
            let cross = pair[5] * c6[1] * c6[2] + pair[2] * c6[4];
            cn_hessian[first * atoms + second] += cross;
            cn_hessian[second * atoms + first] += cross;
            for axis in 0..3 {
                let coordinate = 2.0 * vector[axis];
                let first_coordinate = 3 * first + axis;
                let second_coordinate = 3 * second + axis;
                mixed[first_coordinate * atoms + first] += coordinate * pair[4] * c6[1];
                mixed[first_coordinate * atoms + second] += coordinate * pair[4] * c6[2];
                mixed[second_coordinate * atoms + first] -= coordinate * pair[4] * c6[1];
                mixed[second_coordinate * atoms + second] -= coordinate * pair[4] * c6[2];
            }
        }
    }
    let mut jacobian_cn_hessian = vec![0.0; size * atoms];
    for coordinate in 0..size {
        for second in 0..atoms {
            jacobian_cn_hessian[coordinate * atoms + second] = (0..atoms)
                .map(|first| {
                    cn_jacobian[coordinate * atoms + first] * cn_hessian[first * atoms + second]
                })
                .sum();
        }
    }
    for row in 0..size {
        for column in 0..size {
            let chain = (0..atoms)
                .map(|atom| {
                    mixed[row * atoms + atom] * cn_jacobian[column * atoms + atom]
                        + cn_jacobian[row * atoms + atom] * mixed[column * atoms + atom]
                        + jacobian_cn_hessian[row * atoms + atom]
                            * cn_jacobian[column * atoms + atom]
                })
                .sum::<f64>();
            hessian[row * size + column] += chain;
        }
    }
    for first in 0..atoms {
        for second in 0..first {
            let vector = displacement(positions, first, second);
            let u = vector.iter().map(|value| value * value).sum::<f64>();
            if u > cutoff.cn * cutoff.cn || u < 1.0e-12 {
                continue;
            }
            let distance = u.sqrt();
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
            let radial_first = -count * (1.0 - count) * 16.0 * radius / u;
            let radial_second =
                count * (1.0 - count) * (1.0 - 2.0 * count) * (16.0 * radius).powi(2) / u.powi(2)
                    + 2.0 * count * (1.0 - count) * 16.0 * radius / (u * distance);
            let radial = radial_first / distance;
            let cartesian = (radial_second - radial) / u;
            let scale = dedcn[first] + dedcn[second];
            for row in 0..3 {
                for column in 0..3 {
                    let block = scale
                        * (cartesian * vector[row] * vector[column]
                            + if row == column { radial } else { 0.0 });
                    add_hessian_block(&mut hessian, size, first, second, row, column, block);
                }
            }
        }
    }
    Ok((energy, hessian))
}

fn add_hessian_block(
    hessian: &mut [f64],
    size: usize,
    first: usize,
    second: usize,
    row: usize,
    column: usize,
    block: f64,
) {
    let first_row = 3 * first + row;
    let second_row = 3 * second + row;
    let first_column = 3 * first + column;
    let second_column = 3 * second + column;
    hessian[first_row * size + first_column] += block;
    hessian[second_row * size + second_column] += block;
    hessian[first_row * size + second_column] -= block;
    hessian[second_row * size + first_column] -= block;
}

fn smooth_cutoff_second(u: f64, cutoff: f64, width: f64) -> [f64; 3] {
    let distance = u.sqrt();
    if width <= 0.0 || width >= cutoff || distance <= cutoff - width {
        [1.0, 0.0, 0.0]
    } else if distance >= cutoff {
        [0.0, 0.0, 0.0]
    } else {
        let x = (cutoff - distance) / width;
        let value = x.powi(3) * (10.0 - 15.0 * x + 6.0 * x * x);
        let radial = -30.0 * x * x * (1.0 - x).powi(2) / width;
        let radial_second = 60.0 * x * (1.0 - x) * (1.0 - 2.0 * x) / (width * width);
        [
            value,
            radial / (2.0 * distance),
            (radial_second - radial / distance) / (4.0 * u),
        ]
    }
}

/// ATM-only symmetric `N*N` energy matrix; its entries sum to the ATM energy.
/// Uses default cutoffs and no ghosts. Inputs follow the [module contract](#inputs).
pub fn atm_pairwise(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    partition: WorkPartition,
) -> Result<Vec<f64>, &'static str> {
    atm_pairwise_with_ghosts(numbers, positions, atm, &[], partition)
}

/// ATM-only energy in hartree, with explicit cutoffs, ghosts, and partition.
/// All arguments follow the [module input contract](#inputs).
pub fn atm_energy_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, _, _, weights) = prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let (distances, c6) = atm_pair_data(positions, &elements, &weights);
    let mut energy = 0.0;
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
                energy += atm_energy_triplet(
                    numbers.len(),
                    &distances,
                    &c6,
                    &elements,
                    first,
                    second,
                    third,
                    atm,
                    cutoff,
                );
            }
        }
    }
    Ok(energy)
}

/// As [`atm_pairwise`], with an empty or `N`-entry ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn atm_pairwise_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<Vec<f64>, &'static str> {
    atm_pairwise_with_cutoff(
        numbers,
        positions,
        atm,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

/// As [`atm_pairwise`], with explicit cutoffs and ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn atm_pairwise_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<Vec<f64>, &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, _, _, weights) = prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let (distances, c6) = atm_pair_data(positions, &elements, &weights);
    let mut pairwise = vec![0.0; numbers.len() * numbers.len()];
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
                let energy = atm_energy_triplet(
                    numbers.len(),
                    &distances,
                    &c6,
                    &elements,
                    first,
                    second,
                    third,
                    atm,
                    cutoff,
                );
                for (left, right) in [
                    (first, second),
                    (second, first),
                    (first, third),
                    (third, first),
                    (second, third),
                    (third, second),
                ] {
                    pairwise[left * numbers.len() + right] += energy / 6.0;
                }
            }
        }
    }
    Ok(pairwise)
}

#[allow(clippy::needless_range_loop)]
fn atm_pair_data(
    positions: &[f64],
    elements: &[usize],
    weights: &[[f64; REFERENCES]],
) -> (Vec<f64>, Vec<f64>) {
    let atoms = elements.len();
    let mut distances = vec![0.0; atoms * atoms];
    let mut c6 = vec![0.0; atoms * atoms];
    for first in 0..atoms {
        for second in 0..first {
            let index = first * atoms + second;
            let reverse = second * atoms + first;
            distances[index] = squared_distance(positions, first, second);
            distances[reverse] = distances[index];
            c6[index] = atomic_c6(
                elements[first],
                elements[second],
                &weights[first],
                &weights[second],
            );
            c6[reverse] = c6[index];
        }
    }
    (distances, c6)
}

#[allow(clippy::too_many_arguments)]
fn atm_energy_triplet(
    atoms: usize,
    distances: &[f64],
    c6: &[f64],
    elements: &[usize],
    first: usize,
    second: usize,
    third: usize,
    atm: Atm,
    cutoff: RealspaceCutoff,
) -> f64 {
    let u12 = distances[first * atoms + second];
    let u13 = distances[first * atoms + third];
    let u23 = distances[second * atoms + third];
    if [u12, u13, u23]
        .iter()
        .any(|&distance| !(f64::EPSILON..=cutoff.disp3 * cutoff.disp3).contains(&distance))
    {
        return 0.0;
    }
    let c12 = c6[first * atoms + second];
    let c13 = c6[first * atoms + third];
    let c23 = c6[second * atoms + third];
    let r0 = (4.0_f64 / 3.0).powi(3)
        * value(VDW_RADII, pair_index(elements[first], elements[second]))
        * value(VDW_RADII, pair_index(elements[first], elements[third]))
        * value(VDW_RADII, pair_index(elements[second], elements[third]));
    let product = u12 * u13 * u23;
    let root = product.sqrt();
    let inverse = product.recip();
    let inverse_root = root.recip();
    let angular = 0.375
        * (u12 + u23 - u13)
        * (u12 - u23 + u13)
        * (-u12 + u23 + u13)
        * inverse
        * inverse
        * inverse_root
        + inverse * inverse_root;
    let damping = 1.0 / (1.0 + 6.0 * (r0 / root).powf(atm.alpha / 3.0));
    let switch = smooth_cutoff(u12, cutoff.disp3, cutoff.width3).0
        * smooth_cutoff(u13, cutoff.disp3, cutoff.width3).0
        * smooth_cutoff(u23, cutoff.disp3, cutoff.width3).0;
    atm.s9 * (c12 * c13 * c23).abs().sqrt() * angular * damping * switch
}

/// ATM energy, `3*N` gradient, column-major virial, and symmetric `N*N` pair matrix.
/// Uses default cutoffs and no ghosts. Units are hartree and bohr;
/// all arguments follow the [module input contract](#inputs).
pub fn atm_derivatives(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    partition: WorkPartition,
) -> Result<AtmDerivatives, &'static str> {
    atm_derivatives_with_ghosts(numbers, positions, atm, &[], partition)
}

/// As [`atm_derivatives`], with an empty or `N`-entry ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn atm_derivatives_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<AtmDerivatives, &'static str> {
    atm_derivatives_with_cutoff(
        numbers,
        positions,
        atm,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

/// As [`atm_derivatives`], with explicit cutoffs and ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn atm_derivatives_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<AtmDerivatives, &'static str> {
    analytical_atm_derivatives(numbers, positions, atm, cutoff, ghosts, partition, true)
}

/// ATM energy and row-major `(3*N, 3*N)` Cartesian Hessian (hartree/bohr squared).
/// All arguments follow the [module input contract](#inputs).
pub fn atm_hessian_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<f64>), &'static str> {
    let result =
        analytical_atm_derivatives(numbers, positions, atm, cutoff, ghosts, partition, false)?;
    Ok((result.0, result.3))
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::needless_range_loop)]
fn analytical_atm_derivatives(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
    gradients: bool,
) -> Result<AtmDerivatives, &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, _, coordination, weights) = prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let first_weights: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .zip(weights.iter())
        .map(|((&element, &coordination), weights)| {
            weight_derivatives(element, coordination, weights)
        })
        .collect();
    let second_weights: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .zip(weights.iter())
        .map(|((&element, &coordination), weights)| {
            weight_second_derivatives(element, coordination, weights)
        })
        .collect();
    let atoms = numbers.len();
    let size = positions.len();
    let cn_jacobian = coordination_jacobian(&elements, positions, cutoff.cn);
    let mut pair_c6 = vec![[0.0; 6]; atoms * atoms];
    for first in 0..atoms {
        for second in 0..first {
            let values = atomic_c6_second_derivatives(
                elements[first],
                elements[second],
                &weights[first],
                &weights[second],
                &first_weights[first],
                &first_weights[second],
                &second_weights[first],
                &second_weights[second],
            );
            pair_c6[first * atoms + second] = values;
            pair_c6[second * atoms + first] = [
                values[0], values[2], values[1], values[5], values[4], values[3],
            ];
        }
    }
    let mut energy = 0.0;
    let mut gradient = vec![0.0; size];
    let mut hessian = vec![0.0; size * size];
    let mut dedcn = vec![0.0; atoms];
    let mut mixed = vec![0.0; size * atoms];
    let mut cn_hessian = vec![0.0; atoms * atoms];
    for first in 0..atoms {
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
                let triplet = [first, second, third];
                let pairs = [(first, second), (first, third), (second, third)];
                let distances = pairs.map(|(left, right)| squared_distance(positions, left, right));
                if distances.iter().any(|&distance| {
                    !(f64::EPSILON..=cutoff.disp3 * cutoff.disp3).contains(&distance)
                }) {
                    continue;
                }
                let c6 = pairs.map(|(left, right)| pair_c6[left * atoms + right]);
                let local = local_atm_partials(
                    distances,
                    c6.map(|value| value[0]),
                    pairs.map(|(left, right)| {
                        value(VDW_RADII, pair_index(elements[left], elements[right]))
                    }),
                    atm,
                    cutoff,
                );
                energy += local.value;
                let mut u_jacobian = [[0.0; 3]; 9];
                for (pair, &(left, right)) in pairs.iter().enumerate() {
                    let vector = displacement(positions, left, right);
                    for (local_atom, &atom) in triplet.iter().enumerate() {
                        for axis in 0..3 {
                            u_jacobian[3 * local_atom + axis][pair] = if atom == left {
                                2.0 * vector[axis]
                            } else if atom == right {
                                -2.0 * vector[axis]
                            } else {
                                0.0
                            };
                        }
                    }
                    for row in 0..3 {
                        for column in 0..3 {
                            let block = if row == column {
                                2.0 * local.gradient[pair]
                            } else {
                                0.0
                            };
                            add_hessian_block(&mut hessian, size, left, right, row, column, block);
                        }
                    }
                }
                for local_row in 0..9 {
                    let global_row = 3 * triplet[local_row / 3] + local_row % 3;
                    let contracted: [f64; 3] = std::array::from_fn(|right| {
                        (0..3)
                            .map(|left| {
                                u_jacobian[local_row][left] * local.hessian[left * 6 + right]
                            })
                            .sum()
                    });
                    if gradients {
                        gradient[global_row] += (0..3)
                            .map(|pair| local.gradient[pair] * u_jacobian[local_row][pair])
                            .sum::<f64>();
                    }
                    for local_column in 0..9 {
                        let global_column = 3 * triplet[local_column / 3] + local_column % 3;
                        let value: f64 = (0..3)
                            .map(|right| contracted[right] * u_jacobian[local_column][right])
                            .sum();
                        hessian[global_row * size + global_column] += value;
                    }
                }
                let mut c_jacobian = [[0.0; 3]; 3];
                let mut c_second = [[[0.0; 3]; 3]; 3];
                for (pair, &(left, right)) in pairs.iter().enumerate() {
                    let left_local = triplet.iter().position(|&atom| atom == left).unwrap();
                    let right_local = triplet.iter().position(|&atom| atom == right).unwrap();
                    c_jacobian[pair][left_local] = c6[pair][1];
                    c_jacobian[pair][right_local] = c6[pair][2];
                    c_second[pair][left_local][left_local] = c6[pair][3];
                    c_second[pair][left_local][right_local] = c6[pair][4];
                    c_second[pair][right_local][left_local] = c6[pair][4];
                    c_second[pair][right_local][right_local] = c6[pair][5];
                }
                for local_atom in 0..3 {
                    let atom = triplet[local_atom];
                    dedcn[atom] += (0..3)
                        .map(|pair| local.gradient[3 + pair] * c_jacobian[pair][local_atom])
                        .sum::<f64>();
                    for local_other in 0..3 {
                        let other = triplet[local_other];
                        let mut value = 0.0;
                        for left in 0..3 {
                            for right in 0..3 {
                                value += local.hessian[(3 + left) * 6 + 3 + right]
                                    * c_jacobian[left][local_atom]
                                    * c_jacobian[right][local_other];
                            }
                            value +=
                                local.gradient[3 + left] * c_second[left][local_atom][local_other];
                        }
                        cn_hessian[atom * atoms + other] += value;
                    }
                }
                for local_coordinate in 0..9 {
                    let coordinate = 3 * triplet[local_coordinate / 3] + local_coordinate % 3;
                    let contracted: [f64; 3] = std::array::from_fn(|right| {
                        (0..3)
                            .map(|left| {
                                u_jacobian[local_coordinate][left]
                                    * local.hessian[left * 6 + 3 + right]
                            })
                            .sum()
                    });
                    for local_atom in 0..3 {
                        let atom = triplet[local_atom];
                        for right in 0..3 {
                            mixed[coordinate * atoms + atom] +=
                                contracted[right] * c_jacobian[right][local_atom];
                        }
                    }
                }
            }
        }
    }
    contract_cn_hessian(
        &elements,
        positions,
        cutoff.cn,
        &dedcn,
        &mixed,
        &cn_hessian,
        &cn_jacobian,
        &mut gradient,
        &mut hessian,
        gradients,
    );
    let mut virial = [0.0; 9];
    if gradients {
        for atom in 0..atoms {
            for axis in 0..3 {
                for other in 0..3 {
                    virial[axis + 3 * other] +=
                        gradient[3 * atom + axis] * positions[3 * atom + other];
                }
            }
        }
    }
    Ok((energy, gradient, virial, hessian))
}

#[derive(Clone)]
struct LocalDual {
    value: f64,
    gradient: [f64; 6],
    hessian: [f64; 36],
}

fn local_atm_partials(
    distances: [f64; 3],
    c6: [f64; 3],
    radii: [f64; 3],
    atm: Atm,
    cutoff: RealspaceCutoff,
) -> LocalDual {
    let inverse = distances.map(f64::recip);
    let product = distances.iter().product::<f64>();
    let inverse_cube = 1.0 / (product * product.sqrt());
    let inverse_fifth = inverse_cube / product;
    let sum = distances.iter().sum::<f64>();
    let factors = distances.map(|distance| sum - 2.0 * distance);
    let numerator = factors.iter().product::<f64>();
    let signs: [[f64; 3]; 3] = std::array::from_fn(|factor| {
        std::array::from_fn(|axis| if factor == axis { -1.0 } else { 1.0 })
    });
    let numerator_first: [f64; 3] = std::array::from_fn(|axis| {
        (0..3)
            .map(|factor| {
                signs[factor][axis] * factors[(factor + 1) % 3] * factors[(factor + 2) % 3]
            })
            .sum()
    });
    let angular = 0.375 * numerator * inverse_fifth + inverse_cube;
    let angular_first: [f64; 3] = std::array::from_fn(|axis| {
        0.375 * inverse_fifth * (numerator_first[axis] - 2.5 * numerator * inverse[axis])
            - 1.5 * inverse_cube * inverse[axis]
    });
    let radius = (4.0_f64 / 3.0).powi(3) * radii.iter().product::<f64>();
    let damping = 1.0 / (1.0 + 6.0 * (radius / product.sqrt()).powf(atm.alpha / 3.0));
    let exponent = atm.alpha / 6.0;
    let damping_first = exponent * damping * (1.0 - damping);
    let damping_second = exponent * damping_first * (1.0 - 2.0 * damping);
    let switches =
        distances.map(|distance| smooth_cutoff_second(distance, cutoff.disp3, cutoff.width3));
    let switch = switches.iter().map(|value| value[0]).product::<f64>();
    let switch_first: [f64; 3] = std::array::from_fn(|axis| {
        switches[axis][1] * switches[(axis + 1) % 3][0] * switches[(axis + 2) % 3][0]
    });
    let coefficient = atm.s9 * c6.iter().product::<f64>().sqrt();
    let radial = angular * damping;
    let radial_first: [f64; 3] = std::array::from_fn(|axis| {
        angular_first[axis] * damping + angular * damping_first * inverse[axis]
    });
    let mut result = LocalDual {
        value: coefficient * radial * switch,
        gradient: [0.0; 6],
        hessian: [0.0; 36],
    };
    for row in 0..3 {
        result.gradient[row] =
            coefficient * (radial_first[row] * switch + radial * switch_first[row]);
        result.gradient[3 + row] = 0.5 * result.value / c6[row];
        for column in 0..3 {
            let diagonal = if row == column { 1.0 } else { 0.0 };
            let mut numerator_second = 0.0;
            for first in 0..3 {
                for second in 0..3 {
                    if first != second {
                        numerator_second +=
                            signs[first][row] * signs[second][column] * factors[3 - first - second];
                    }
                }
            }
            let angular_second = 0.375
                * inverse_fifth
                * (numerator_second
                    - 2.5
                        * (numerator_first[row] * inverse[column]
                            + numerator_first[column] * inverse[row])
                    + numerator * (6.25 + 2.5 * diagonal) * inverse[row] * inverse[column])
                + inverse_cube * (2.25 + 1.5 * diagonal) * inverse[row] * inverse[column];
            let radial_second = angular_second * damping
                + damping_first
                    * (angular_first[row] * inverse[column] + angular_first[column] * inverse[row])
                + angular
                    * (damping_second - diagonal * damping_first)
                    * inverse[row]
                    * inverse[column];
            let switch_second = if row == column {
                switches[row][2] * switches[(row + 1) % 3][0] * switches[(row + 2) % 3][0]
            } else {
                switches[row][1] * switches[column][1] * switches[3 - row - column][0]
            };
            result.hessian[row * 6 + column] = coefficient
                * (radial_second * switch
                    + radial_first[row] * switch_first[column]
                    + radial_first[column] * switch_first[row]
                    + radial * switch_second);
            let mixed = 0.5 * result.gradient[row] / c6[column];
            result.hessian[row * 6 + 3 + column] = mixed;
            result.hessian[(3 + column) * 6 + row] = mixed;
            result.hessian[(3 + row) * 6 + 3 + column] =
                result.value * (0.25 - 0.5 * diagonal) / (c6[row] * c6[column]);
        }
    }
    result
}

fn coordination_jacobian(elements: &[usize], positions: &[f64], cutoff: f64) -> Vec<f64> {
    let atoms = elements.len();
    let mut jacobian = vec![0.0; 3 * atoms * atoms];
    for first in 0..atoms {
        for second in 0..first {
            let vector = displacement(positions, first, second);
            let u = squared_distance(positions, first, second);
            if u > cutoff * cutoff || u < 1.0e-12 {
                continue;
            }
            let distance = u.sqrt();
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
            let radial = -count * (1.0 - count) * 16.0 * radius / u;
            for (axis, &coordinate) in vector.iter().enumerate() {
                let derivative = radial * coordinate / distance;
                let first_coordinate = 3 * first + axis;
                let second_coordinate = 3 * second + axis;
                jacobian[first_coordinate * atoms + first] += derivative;
                jacobian[second_coordinate * atoms + first] -= derivative;
                jacobian[first_coordinate * atoms + second] += derivative;
                jacobian[second_coordinate * atoms + second] -= derivative;
            }
        }
    }
    jacobian
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::needless_range_loop)]
fn contract_cn_hessian(
    elements: &[usize],
    positions: &[f64],
    cutoff: f64,
    dedcn: &[f64],
    mixed: &[f64],
    cn_hessian: &[f64],
    jacobian: &[f64],
    gradient: &mut [f64],
    hessian: &mut [f64],
    gradients: bool,
) {
    let atoms = elements.len();
    let size = positions.len();
    if gradients {
        for coordinate in 0..size {
            gradient[coordinate] += (0..atoms)
                .map(|atom| dedcn[atom] * jacobian[coordinate * atoms + atom])
                .sum::<f64>();
        }
    }
    let mut jq = vec![0.0; size * atoms];
    for coordinate in 0..size {
        for second in 0..atoms {
            jq[coordinate * atoms + second] = (0..atoms)
                .map(|first| {
                    jacobian[coordinate * atoms + first] * cn_hessian[first * atoms + second]
                })
                .sum();
        }
    }
    for row in 0..size {
        for column in 0..size {
            hessian[row * size + column] += (0..atoms)
                .map(|atom| {
                    mixed[row * atoms + atom] * jacobian[column * atoms + atom]
                        + jacobian[row * atoms + atom] * mixed[column * atoms + atom]
                        + jq[row * atoms + atom] * jacobian[column * atoms + atom]
                })
                .sum::<f64>();
        }
    }
    for first in 0..atoms {
        for second in 0..first {
            let vector = displacement(positions, first, second);
            let u = squared_distance(positions, first, second);
            if u > cutoff * cutoff || u < 1.0e-12 {
                continue;
            }
            let distance = u.sqrt();
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
            let radial_first = -count * (1.0 - count) * 16.0 * radius / u;
            let radial_second =
                count * (1.0 - count) * (1.0 - 2.0 * count) * (16.0 * radius).powi(2) / u.powi(2)
                    + 2.0 * count * (1.0 - count) * 16.0 * radius / (u * distance);
            let radial = radial_first / distance;
            let cartesian = (radial_second - radial) / u;
            let scale = dedcn[first] + dedcn[second];
            for row in 0..3 {
                for column in 0..3 {
                    let block = scale
                        * (cartesian * vector[row] * vector[column]
                            + if row == column { radial } else { 0.0 });
                    add_hessian_block(hessian, size, first, second, row, column, block);
                }
            }
        }
    }
}

/// ATM energy, `3*N` gradient, and column-major virial without allocating pair output.
/// All arguments follow the [module input contract](#inputs).
pub fn atm_gradient_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    atm: Atm,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<AtmGradient, &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, _, coordination, weights) = prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let derivatives: Vec<_> = elements
        .iter()
        .zip(coordination.iter())
        .zip(weights.iter())
        .map(|((&element, &coordination), weights)| {
            weight_derivatives(element, coordination, weights)
        })
        .collect();
    let atoms = numbers.len();
    let mut distance = vec![0.0; atoms * atoms];
    let mut c6 = vec![0.0; atoms * atoms];
    let mut dc6 = vec![0.0; atoms * atoms];
    for first in 0..atoms {
        for second in 0..first {
            let ij = first * atoms + second;
            let ji = second * atoms + first;
            distance[ij] = squared_distance(positions, first, second);
            distance[ji] = distance[ij];
            let values = atomic_c6_derivatives(
                elements[first],
                elements[second],
                &weights[first],
                &weights[second],
                &derivatives[first],
                &derivatives[second],
            );
            c6[ij] = values.0;
            c6[ji] = values.0;
            dc6[ij] = values.1;
            dc6[ji] = values.2;
        }
    }
    let mut energy = 0.0;
    let mut gradient = vec![0.0; positions.len()];
    let mut virial = [0.0; 9];
    let mut dedcn = vec![0.0; atoms];
    for first in 0..atoms {
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
                let pairs = [(first, second), (first, third), (second, third)];
                let u = pairs.map(|(left, right)| distance[left * atoms + right]);
                if u.iter().any(|&distance| {
                    !(f64::EPSILON..=cutoff.disp3 * cutoff.disp3).contains(&distance)
                }) {
                    continue;
                }
                let pair_c6 = pairs.map(|(left, right)| c6[left * atoms + right]);
                let radii = pairs.map(|(left, right)| {
                    value(VDW_RADII, pair_index(elements[left], elements[right]))
                });
                let (value, dedu, dedc6) = atm_triplet_derivatives(u, pair_c6, radii, atm, cutoff);
                energy += value;
                for (index, &(left, right)) in pairs.iter().enumerate() {
                    let vector = displacement(positions, left, right);
                    for axis in 0..3 {
                        let component = 2.0 * dedu[index] * vector[axis];
                        gradient[3 * left + axis] += component;
                        gradient[3 * right + axis] -= component;
                        for other in 0..3 {
                            virial[axis + 3 * other] += component * vector[other];
                        }
                    }
                    dedcn[left] += dedc6[index] * dc6[left * atoms + right];
                    dedcn[right] += dedc6[index] * dc6[right * atoms + left];
                }
            }
        }
    }
    for first in 0..atoms {
        for second in 0..first {
            let vector = displacement(positions, first, second);
            let distance2 = squared_distance(positions, first, second);
            if distance2 > cutoff.cn * cutoff.cn || distance2 < 1.0e-12 {
                continue;
            }
            let distance = distance2.sqrt();
            let radius =
                value(COVALENT_RADII, elements[first]) + value(COVALENT_RADII, elements[second]);
            let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
            let scale = -count * (1.0 - count) * 16.0 * radius / (distance2 * distance)
                * (dedcn[first] + dedcn[second]);
            for axis in 0..3 {
                let component = scale * vector[axis];
                gradient[3 * first + axis] += component;
                gradient[3 * second + axis] -= component;
                for other in 0..3 {
                    virial[axis + 3 * other] += component * vector[other];
                }
            }
        }
    }
    Ok((energy, gradient, virial))
}

fn atm_triplet_derivatives(
    u: [f64; 3],
    c6: [f64; 3],
    radii: [f64; 3],
    atm: Atm,
    cutoff: RealspaceCutoff,
) -> (f64, [f64; 3], [f64; 3]) {
    let product = u.iter().product::<f64>();
    let factors = [u[0] + u[2] - u[1], u[0] - u[2] + u[1], -u[0] + u[2] + u[1]];
    let numerator = factors.iter().product::<f64>();
    let inverse = product.recip();
    let inverse_root = product.sqrt().recip();
    let inverse_3_2 = inverse * inverse_root;
    let inverse_5_2 = inverse_3_2 * inverse;
    let inverse_7_2 = inverse_5_2 * inverse;
    let angular = 0.375 * numerator * inverse_5_2 + inverse_3_2;
    let r0 = (4.0_f64 / 3.0).powi(3) * radii.iter().product::<f64>();
    let damping = 1.0 / (1.0 + 6.0 * (r0 / product.sqrt()).powf(atm.alpha / 3.0));
    let switches = u.map(|value| smooth_cutoff(value, cutoff.disp3, cutoff.width3));
    let switch = switches.iter().map(|item| item.0).product::<f64>();
    let coefficient = c6.iter().product::<f64>().abs().sqrt();
    let energy = atm.s9 * coefficient * angular * damping * switch;
    let mut dedu = [0.0; 3];
    for index in 0..3 {
        let dproduct = product / u[index];
        let dnumerator = match index {
            0 => factors[1] * factors[2] + factors[0] * factors[2] - factors[0] * factors[1],
            1 => -factors[1] * factors[2] + factors[0] * factors[2] + factors[0] * factors[1],
            _ => factors[1] * factors[2] - factors[0] * factors[2] + factors[0] * factors[1],
        };
        let dangular = 0.375
            * (dnumerator * inverse_5_2 - 2.5 * numerator * inverse_7_2 * dproduct)
            - 1.5 * inverse_5_2 * dproduct;
        let ddamping = damping * (1.0 - damping) * atm.alpha / (6.0 * u[index]);
        let dswitch = switches[index].1 * switches[(index + 1) % 3].0 * switches[(index + 2) % 3].0;
        dedu[index] = atm.s9
            * coefficient
            * (dangular * damping * switch
                + angular * ddamping * switch
                + angular * damping * dswitch);
    }
    let dedc6 = std::array::from_fn(|index| energy / (2.0 * c6[index]));
    (energy, dedu, dedc6)
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "ELEMENTS is the compile-time constant 103"
)]
fn prepare_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    cn_cutoff: f64,
) -> Result<Prepared, &'static str> {
    if positions.len() != numbers.len() * 3 {
        return Err("positions must contain three coordinates per atom");
    }
    let elements: Vec<_> = numbers
        .iter()
        .map(|&number| {
            if (1..=ELEMENTS as i32).contains(&number) {
                Ok(number as usize - 1)
            } else {
                Err("D3 supports atomic numbers 1 through 103")
            }
        })
        .collect::<Result<_, _>>()?;
    let mut unique_elements = Vec::new();
    let species = elements
        .iter()
        .map(|element| {
            unique_elements
                .iter()
                .position(|known| known == element)
                .unwrap_or_else(|| {
                    unique_elements.push(*element);
                    unique_elements.len() - 1
                })
        })
        .collect();
    let contributions = crate::parallel::map(numbers.len(), 256, |start, stride| {
        let mut coordination = vec![0.0; numbers.len()];
        for first in (start..numbers.len()).step_by(stride) {
            for second in 0..first {
                let vector = displacement(positions, first, second);
                let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
                if distance2 > cn_cutoff * cn_cutoff || distance2 < 1.0e-12 {
                    continue;
                }
                let radius = value(COVALENT_RADII, elements[first])
                    + value(COVALENT_RADII, elements[second]);
                let count = 1.0 / (1.0 + (-16.0 * (radius / distance2.sqrt() - 1.0)).exp());
                coordination[first] += count;
                coordination[second] += count;
            }
        }
        coordination
    });
    let mut contributions = contributions.into_iter();
    let mut coordination = contributions.next().unwrap();
    for local in contributions {
        for (total, value) in coordination.iter_mut().zip(local) {
            *total += value;
        }
    }
    let weights = elements
        .iter()
        .zip(coordination.iter())
        .map(|(&element, &coordination)| weights(element, coordination))
        .collect();
    Ok((elements, species, coordination, weights))
}

/// Return `N` coordination numbers and the symmetric row-major `N*N` C6 matrix.
/// C6 units are hartree bohr^6. All arguments follow the
/// [module input contract](#inputs); supports zero to three periodic axes.
pub fn properties(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    cn_cutoff: f64,
) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    properties_with_model(numbers, positions, lattice, periodic, cn_cutoff, Model::D3)
}

fn properties_with_model(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    cn_cutoff: f64,
    model: Model,
) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    model.validate(numbers)?;
    if numbers.is_empty() || positions.iter().any(|value| !value.is_finite()) {
        return Err("D3 requires a nonempty finite structure");
    }
    if !cn_cutoff.is_finite() || cn_cutoff < 0.0 {
        return Err("D3 CN cutoff must be finite and nonnegative");
    }
    let (elements, _, coordination, weights) = if periodic.iter().any(|&active| active) {
        realspace::prepare(
            numbers,
            positions,
            lattice.ok_or("periodic D3 requires a lattice")?,
            periodic,
            cn_cutoff,
            cn_cutoff,
        )?
        .0
    } else {
        prepare_with_cutoff(numbers, positions, cn_cutoff)?
    };
    let atoms = numbers.len();
    let mut c6 = vec![0.0; atoms * atoms];
    for first in 0..atoms {
        for second in 0..=first {
            let coefficient = if model == Model::D3S {
                smooth::pair_coefficients(
                    elements[first],
                    elements[second],
                    coordination[first],
                    coordination[second],
                )[0]
            } else {
                atomic_c6(
                    elements[first],
                    elements[second],
                    &weights[first],
                    &weights[second],
                )
            };
            c6[first * atoms + second] = coefficient;
            c6[second * atoms + first] = coefficient;
        }
    }
    Ok((coordination, c6))
}

/// D3 properties and analytical Jacobians, flattened property-major.
pub struct PropertyResponse {
    /// `N` dimensionless coordination numbers.
    pub coordination: Vec<f64>,
    /// Row-major `N*N` C6 coefficients in hartree bohr^6.
    pub c6: Vec<f64>,
    /// CN Jacobian: `N` rows of `3*N` Cartesian derivatives (1/bohr).
    pub coordination_cartesian: Vec<f64>,
    /// CN strain Jacobian: `N` rows of 9 derivatives, `row + 3*column`.
    pub coordination_strain: Vec<f64>,
    /// C6 Jacobian: `N*N` rows of `3*N` derivatives (hartree bohr^5).
    pub c6_cartesian: Vec<f64>,
    /// C6 strain Jacobian: `N*N` rows of 9 derivatives (hartree bohr^6).
    pub c6_strain: Vec<f64>,
}

/// Inputs follow [`properties`]; returns CN/C6 values and their geometry responses.
/// Property-major derivatives: coordinate = 3*atom+axis, strain = row+3*column.
/// Strain transforms both coordinates and lattice; hard cutoffs must not be crossed.
/// The full Cartesian C6 Jacobian requires 3*N^3 elements.
pub fn property_response(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    cn_cutoff: f64,
) -> Result<PropertyResponse, &'static str> {
    property_response_with_model(numbers, positions, lattice, periodic, cn_cutoff, Model::D3)
}

fn property_response_with_model(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    cn_cutoff: f64,
    model: Model,
) -> Result<PropertyResponse, &'static str> {
    let response_size = numbers
        .len()
        .checked_pow(3)
        .and_then(|size| size.checked_mul(3))
        .filter(|&size| size <= isize::MAX as usize / std::mem::size_of::<f64>())
        .ok_or("D3 property response dimensions overflow")?;
    let (coordination, c6) =
        properties_with_model(numbers, positions, lattice, periodic, cn_cutoff, model)?;
    let atoms = numbers.len();
    let coordinates = positions.len();
    let elements: Vec<_> = numbers.iter().map(|&number| number as usize - 1).collect();
    let mut result = PropertyResponse {
        coordination,
        c6,
        coordination_cartesian: vec![0.0; atoms * coordinates],
        coordination_strain: vec![0.0; atoms * 9],
        c6_cartesian: vec![0.0; response_size],
        c6_strain: vec![0.0; atoms * atoms * 9],
    };
    if periodic.iter().any(|&active| active) {
        let (_, neighbors) = realspace::prepare(
            numbers,
            positions,
            lattice.unwrap(),
            periodic,
            cn_cutoff,
            cn_cutoff,
        )?;
        for (first, neighbors) in neighbors.iter().enumerate() {
            for &(second, vector, distance2) in neighbors {
                let distance = distance2.sqrt();
                let radius = value(COVALENT_RADII, elements[first])
                    + value(COVALENT_RADII, elements[second]);
                let count = 1.0 / (1.0 + (-16.0 * (radius / distance - 1.0)).exp());
                let radial = -count * (1.0 - count) * 16.0 * radius / (distance2 * distance);
                for axis in 0..3 {
                    let derivative = radial * vector[axis];
                    result.coordination_cartesian[first * coordinates + 3 * first + axis] +=
                        derivative;
                    result.coordination_cartesian[first * coordinates + 3 * second + axis] -=
                        derivative;
                    for (column, &coordinate) in vector.iter().enumerate() {
                        result.coordination_strain[first * 9 + axis + 3 * column] +=
                            derivative * coordinate;
                    }
                }
            }
        }
    } else {
        let jacobian = coordination_jacobian(&elements, positions, cn_cutoff);
        for atom in 0..atoms {
            for coordinate in 0..coordinates {
                let derivative = jacobian[coordinate * atoms + atom];
                result.coordination_cartesian[atom * coordinates + coordinate] = derivative;
                for column in 0..3 {
                    result.coordination_strain[atom * 9 + coordinate % 3 + 3 * column] +=
                        derivative * positions[coordinate / 3 * 3 + column];
                }
            }
        }
    }
    let weights: Vec<_> = elements
        .iter()
        .zip(&result.coordination)
        .map(|(&element, &cn)| weights(element, cn))
        .collect();
    let derivatives: Vec<_> = elements
        .iter()
        .zip(&result.coordination)
        .zip(&weights)
        .map(|((&element, &cn), weights)| weight_derivatives(element, cn, weights))
        .collect();
    for first in 0..atoms {
        for second in 0..atoms {
            let (_, left, right) = if model == Model::D3S {
                let partial = smooth::pair_coefficients(
                    elements[first],
                    elements[second],
                    result.coordination[first],
                    result.coordination[second],
                );
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
            for (output, input, size) in [
                (
                    &mut result.c6_cartesian,
                    &result.coordination_cartesian,
                    coordinates,
                ),
                (&mut result.c6_strain, &result.coordination_strain, 9),
            ] {
                for component in 0..size {
                    output[(first * atoms + second) * size + component] = left
                        * input[first * size + component]
                        + right * input[second * size + component];
                }
            }
        }
    }
    Ok(result)
}

#[test]
fn property_responses_match_differences() {
    let numbers = [6, 8, 1];
    let positions = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, -1.0, 2.0, 0.5];
    let lattice = [7.0, 0.0, 0.0, 1.0, 8.0, 0.0, 0.5, 0.2, 9.0];
    let step = 1.0e-5;
    for periodic in [[false; 3], [true, false, false], [true; 3]] {
        let response =
            property_response(&numbers, &positions, Some(&lattice), periodic, 12.0).unwrap();
        for component in 0..18 {
            let evaluate = |sign: f64| {
                let mut xyz = positions;
                let mut cell = lattice;
                if component < 9 {
                    xyz[component] += sign * step;
                } else {
                    let row = (component - 9) % 3;
                    let column = (component - 9) / 3;
                    for atom in 0..3 {
                        xyz[3 * atom + row] += sign * step * positions[3 * atom + column];
                        cell[3 * atom + row] += sign * step * lattice[3 * atom + column];
                    }
                }
                properties(&numbers, &xyz, Some(&cell), periodic, 12.0).unwrap()
            };
            let plus = evaluate(1.0);
            let minus = evaluate(-1.0);
            for (plus, minus, cartesian, strain) in [
                (
                    &plus.0,
                    &minus.0,
                    &response.coordination_cartesian,
                    &response.coordination_strain,
                ),
                (
                    &plus.1,
                    &minus.1,
                    &response.c6_cartesian,
                    &response.c6_strain,
                ),
            ] {
                for property in 0..plus.len() {
                    let analytical = if component < 9 {
                        cartesian[property * 9 + component]
                    } else {
                        strain[property * 9 + component - 9]
                    };
                    assert!(
                        (analytical - (plus[property] - minus[property]) / (2.0 * step)).abs()
                            < 2.0e-6
                    );
                }
            }
        }
    }
}

/// Symmetric row-major `N*N` two-body energy matrix in hartree; sum for total energy.
/// Each off-diagonal entry holds half a pair energy. Uses default cutoffs,
/// no ghosts, and no ATM. Inputs follow the [module contract](#inputs).
pub fn pairwise_partitioned(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    partition: WorkPartition,
) -> Result<Vec<f64>, &'static str> {
    pairwise_partitioned_with_ghosts(numbers, positions, damping, &[], partition)
}

/// As [`pairwise_partitioned`], with an empty or `N`-entry ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn pairwise_partitioned_with_ghosts(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<Vec<f64>, &'static str> {
    pairwise_partitioned_with_cutoff(
        numbers,
        positions,
        damping,
        RealspaceCutoff::default(),
        ghosts,
        partition,
    )
}

/// As [`pairwise_partitioned`], with explicit cutoffs and ghost mask.
/// All arguments follow the [module input contract](#inputs).
pub fn pairwise_partitioned_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<Vec<f64>, &'static str> {
    let mut pairwise = vec![0.0; numbers.len() * numbers.len()];
    for (first, second, energy) in
        pair_energies::<true>(numbers, positions, damping, cutoff, ghosts, partition)?.1
    {
        pairwise[first * numbers.len() + second] = 0.5 * energy;
        pairwise[second * numbers.len() + first] = 0.5 * energy;
    }
    Ok(pairwise)
}

#[allow(clippy::type_complexity)]
fn pair_energies<const PAIRWISE: bool>(
    numbers: &[i32],
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
) -> Result<(f64, Vec<(usize, usize, f64)>), &'static str> {
    validate_ghosts(numbers, ghosts)?;
    let (elements, species, _, weights) = prepare_with_cutoff(numbers, positions, cutoff.cn)?;
    let mut unique_elements = Vec::new();
    for &element in &elements {
        if !unique_elements.contains(&element) {
            unique_elements.push(element);
        }
    }
    let contracted: Vec<Vec<[f64; REFERENCES]>> = elements
        .iter()
        .enumerate()
        .map(|(atom, &element)| {
            unique_elements
                .iter()
                .map(|&other| {
                    let mut result = [0.0; REFERENCES];
                    for (other_reference, total) in
                        result.iter_mut().enumerate().take(reference_count(other))
                    {
                        for (reference, &weight) in weights[atom]
                            .iter()
                            .enumerate()
                            .take(reference_count(element))
                        {
                            *total +=
                                weight * reference_c6(element, other, reference, other_reference);
                        }
                    }
                    result
                })
                .collect()
        })
        .collect();
    let contributions = crate::parallel::map(numbers.len(), 256, |start, stride| {
        pair_energy_range::<PAIRWISE>(
            (start..numbers.len()).step_by(stride),
            positions,
            damping,
            cutoff,
            ghosts,
            partition,
            &elements,
            &species,
            &weights,
            &contracted,
        )
    });
    let mut energy = 0.0;
    let mut pairs = Vec::new();
    for (contribution, values) in contributions {
        energy += contribution;
        pairs.extend(values);
    }
    Ok((energy, pairs))
}

#[allow(clippy::too_many_arguments)]
fn pair_energy_range<const PAIRWISE: bool>(
    range: impl Iterator<Item = usize>,
    positions: &[f64],
    damping: Damping,
    cutoff: RealspaceCutoff,
    ghosts: &[bool],
    partition: WorkPartition,
    elements: &[usize],
    species: &[usize],
    weights: &[[f64; REFERENCES]],
    contracted: &[Vec<[f64; REFERENCES]>],
) -> (f64, Vec<(usize, usize, f64)>) {
    let mut result = Vec::new();
    let mut total = 0.0;
    for first in range {
        for second in 0..first {
            if ghosts.get(first) == Some(&true)
                || ghosts.get(second) == Some(&true)
                || !partition.owns_pair(first, second)
            {
                continue;
            }
            let distance2 = squared_distance(positions, first, second);
            if distance2 > cutoff.disp2 * cutoff.disp2 || distance2 <= f64::EPSILON {
                continue;
            }
            let c6 = contracted[first][species[second]]
                .iter()
                .zip(&weights[second])
                .take(reference_count(elements[second]))
                .map(|(left, right)| left * right)
                .sum();
            let energy = kernel(
                damping,
                elements[first],
                elements[second],
                species[first],
                species[second],
                distance2,
                c6,
            )
            .0;
            let energy = smooth_cutoff(distance2, cutoff.disp2, cutoff.width2).0 * energy;
            if PAIRWISE {
                result.push((first, second, energy));
            } else {
                total += energy;
            }
        }
    }
    (total, result)
}

fn smooth_cutoff(distance2: f64, cutoff: f64, width: f64) -> (f64, f64) {
    if width >= cutoff {
        (1.0, 0.0)
    } else {
        crate::geometry::smooth_cutoff(distance2, cutoff, width)
    }
}

fn validate_ghosts(numbers: &[i32], ghosts: &[bool]) -> Result<(), &'static str> {
    if ghosts.is_empty() || ghosts.len() == numbers.len() {
        Ok(())
    } else {
        Err("ghost mask must contain one entry per atom")
    }
}
