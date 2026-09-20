use super::{pair_index, value, WorkPartition, VDW_RADII};
use crate::geometry::{lattice_repetitions, lattice_translation};

const ELEMENTS: usize = 36;
const TABLES: usize = 23;
const TABLE_STRIDE: usize = 2 * ELEMENTS;
const GCP_TABLES: &[u8; (TABLES * TABLE_STRIDE + ELEMENTS) * 8] =
    include_bytes!("../../assets/gcp_tables.bin");

#[derive(Clone)]
pub struct Gcp {
    numbers: Vec<i32>,
    pub(crate) effective: Vec<usize>,
    pub(crate) emiss: Vec<f64>,
    pub(crate) virtuals: Vec<f64>,
    pub(crate) slater: Vec<f64>,
    pub(crate) rvdw: Option<Vec<f64>>,
    pub(crate) rvdw_srb: Option<Vec<f64>>,
    eta: f64,
    eta_spec: f64,
    pub(crate) sigma: f64,
    pub(crate) alpha: f64,
    pub(crate) beta: f64,
    pub(crate) damp: bool,
    pub(crate) srb: bool,
    pub(crate) base: bool,
    pub(crate) rscal: f64,
    pub(crate) qscal: f64,
    pub(crate) dmp_scal: f64,
    pub(crate) dmp_exp: f64,
}

#[derive(Clone, Copy)]
pub struct GcpCutoff {
    pub gcp: f64,
    pub srb: f64,
}

impl Default for GcpCutoff {
    fn default() -> Self {
        Self {
            gcp: 60.0,
            srb: 60.0,
        }
    }
}

pub struct GcpResult {
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub virial: [f64; 9],
    pub hessian: Vec<f64>,
}

#[derive(Clone, Copy, PartialEq)]
enum Method {
    Unknown,
    Hf,
    Dft,
    Gga,
    B3lyp,
    Blyp,
    Pbe,
    Tpss,
    Pw6b95,
    Hf3c,
    Pbeh3c,
    Hse3c,
    B973c,
    B3pbe3c,
    R2scan3c,
}

#[derive(Clone, Copy)]
enum Basis {
    Sv = 1,
    SvP = 2,
    Svx = 3,
    Svp = 4,
    Minis = 6,
    G631 = 7,
    Tz = 8,
    Def1tzvp = 9,
    Ccdz = 10,
    Accdz = 11,
    Pobtz = 12,
    Minix = 13,
    Gcore = 14,
    TwoG = 15,
    Dzp = 16,
    Dz = 17,
    Msvp = 18,
    Lanl = 19,
    Pbeh3c = 20,
    Def2mtzvpp = 21,
    Def2mtzvp = 22,
}

pub fn load(numbers: &[i32], method: Option<&str>, basis: Option<&str>) -> Option<Gcp> {
    let method_id = method.map(method_id).unwrap_or(Method::Unknown);
    let basis_id = basis
        .and_then(basis_id)
        .or_else(|| method.and_then(basis_id));
    if method_id == Method::Unknown && basis_id.is_none() {
        return None;
    }
    let effective: Vec<_> = numbers
        .iter()
        .map(|&number| effective_number(number).filter(|&number| number < ELEMENTS))
        .collect::<Option<_>>()?;
    let mut sigma = 0.0;
    let mut eta = 0.0;
    let mut eta_spec = 0.0;
    let mut alpha = 0.0;
    let mut beta = 0.0;
    let mut damp = false;
    let mut base = false;
    let mut rscal = 0.0;
    let mut qscal = 0.0;
    if let Some(basis) = basis_id {
        match basis {
            Basis::Sv => match method_id {
                Method::Hf => (sigma, eta, alpha, beta) = (0.1724, 1.2804, 0.8568, 1.2342),
                Method::B3lyp | Method::Pw6b95 => {
                    (sigma, eta, alpha, beta) = (0.4048, 1.1626, 0.8652, 1.2375)
                }
                Method::Gga | Method::Tpss | Method::Blyp => {
                    (sigma, eta, alpha, beta) = (0.2727, 1.4022, 0.8055, 1.3)
                }
                _ => {}
            },
            Basis::SvP => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.1373, 1.4271, 0.8141, 1.2760);
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.2424, 1.2371, 0.6076, 1.4078);
                }
            }
            Basis::Svx if is_dft(method_id) => {
                (sigma, eta, alpha, beta) = (0.1861, 1.3200, 0.6171, 1.4019)
            }
            Basis::Svp => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.2054, 1.3157, 0.8136, 1.2572);
                } else if method_id == Method::Tpss {
                    (sigma, eta, alpha, beta) = (0.6647, 1.3306, 1.0792, 1.1651);
                } else if method_id == Method::Pw6b95 {
                    (sigma, eta, alpha, beta) = (0.3098, 1.2373, 0.6896, 1.3347);
                } else if is_hybrid(method_id) {
                    (sigma, eta, alpha, beta) = (0.2990, 1.2605, 0.6438, 1.3694);
                } else if is_gga(method_id) {
                    (sigma, eta, alpha, beta) = (0.6823, 1.2491, 0.8225, 1.2811);
                }
            }
            Basis::Minis => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.1290, 1.1526, 1.1549, 1.1763);
                } else if method_id == Method::Tpss {
                    (sigma, eta, alpha, beta) = (0.22982, 1.35401, 1.47633, 1.11300);
                } else if method_id == Method::Pw6b95 {
                    (sigma, eta, alpha, beta) = (0.21054, 1.25458, 1.35003, 1.14061);
                } else if is_gga(method_id) {
                    (sigma, eta, alpha, beta) = (0.1566, 1.0271, 1.0732, 1.1968);
                } else if is_hybrid(method_id) {
                    (sigma, eta, alpha, beta) = (0.2059, 0.9722, 1.1961, 1.1456);
                }
            }
            Basis::G631 => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.2048, 1.5652, 0.9447, 1.2100);
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.3405, 1.6127, 0.8589, 1.2830);
                }
            }
            Basis::Tz => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.3127, 1.9914, 1.0216, 1.2833);
                } else if is_hybrid(method_id) {
                    (sigma, eta, alpha, beta) = (0.2905, 2.2495, 0.8120, 1.4412);
                } else if is_gga(method_id) {
                    (sigma, eta, alpha, beta) = (0.1182, 1.0631, 1.0510, 1.1287);
                }
            }
            Basis::Def1tzvp => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.2600, 2.2448, 0.7998, 1.4381);
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.2393, 2.2247, 0.8185, 1.4298);
                }
            }
            Basis::Ccdz => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.4416, 1.5185, 0.6902, 1.3713);
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.5383, 1.6482, 0.6230, 1.4523);
                }
            }
            Basis::Accdz => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.0748, 0.0663, 0.3811, 1.0155);
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.1465, 0.0500, 0.6003, 0.8761);
                }
            }
            Basis::Pobtz if is_dft(method_id) => {
                (sigma, eta, alpha, beta) = (0.1300, 1.3743, 0.4792, 1.3962)
            }
            Basis::Minix => {
                if matches!(method_id, Method::Hf | Method::Hf3c) {
                    (sigma, eta, alpha, beta) = (0.1290, 1.1526, 1.1549, 1.1763);
                    base = method_id == Method::Hf3c;
                    if base {
                        (rscal, qscal) = (0.7, 0.03);
                    }
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.2059, 0.9722, 1.1961, 1.1456);
                }
            }
            Basis::TwoG if method_id == Method::Hf => {
                (sigma, eta, alpha, beta) = (0.2461, 1.1616, 0.7335, 1.4709)
            }
            Basis::Dzp => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.1443, 1.4547, 0.3711, 1.6300);
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.2687, 1.4634, 0.3513, 1.6880);
                }
            }
            Basis::Dz => {
                if method_id == Method::Hf {
                    (sigma, eta, alpha, beta) = (0.1059, 1.4554, 0.3711, 1.6342);
                } else if is_dft(method_id) {
                    (sigma, eta, alpha, beta) = (0.2687, 1.4634, 0.3513, 1.6880);
                }
            }
            Basis::Lanl if is_dft(method_id) => {
                (sigma, eta, alpha, beta) = (0.3405, 1.6127, 0.8589, 1.2830)
            }
            Basis::Def2mtzvp if method_id == Method::B3pbe3c => {
                (sigma, eta, alpha, beta) = (1.0, 2.98561, 0.3011, 2.4405)
            }
            Basis::Pbeh3c if method_id == Method::Pbeh3c => {
                (sigma, eta, alpha, beta) = (1.0, 1.32492, 0.27649, 1.95600);
                damp = true;
            }
            Basis::Pbeh3c if method_id == Method::Hse3c => {
                (sigma, eta, alpha, beta) = (1.0, 1.32378, 0.28314, 1.94527);
                damp = true;
            }
            Basis::Def2mtzvpp if method_id == Method::R2scan3c => {
                (sigma, eta, eta_spec, alpha, beta) = (1.0, 1.3150, 1.15, 0.9410, 1.4636);
                damp = true;
            }
            _ => {}
        }
    }
    let table = basis_id;
    let emiss = effective
        .iter()
        .map(|&element| table.map_or(0.0, |table| table_value(table, 0, element)))
        .collect();
    let virtuals = effective
        .iter()
        .map(|&element| {
            let basis_functions = table.map_or(0.0, |table| table_value(table, 1, element));
            if matches!(table, Some(Basis::Def2mtzvpp)) {
                match element + 1 {
                    6 => 3.0,
                    7 | 8 => 0.5,
                    _ => 1.0,
                }
            } else {
                basis_functions - 0.5 * number_of_electrons(element + 1, false) as f64
            }
        })
        .collect();
    let slater = effective
        .iter()
        .map(|&element| {
            let scale = if eta_spec > 0.0 && element + 1 > 10 {
                eta * eta_spec
            } else {
                eta
            };
            scale * value(GCP_TABLES, TABLES * TABLE_STRIDE + element)
        })
        .collect();
    let srb = method_id == Method::B973c;
    if srb {
        (rscal, qscal) = (10.0, 0.08);
    }
    Some(Gcp {
        numbers: numbers.to_vec(),
        effective,
        emiss,
        virtuals,
        slater,
        rvdw: None,
        rvdw_srb: None,
        eta,
        eta_spec,
        sigma,
        alpha,
        beta,
        damp,
        srb,
        base,
        rscal,
        qscal,
        dmp_scal: 4.0,
        dmp_exp: 6.0,
    })
}

impl Gcp {
    pub fn controls(&self) -> (f64, bool, bool) {
        (self.eta, self.base, self.srb)
    }

    pub fn set_controls(&mut self, eta: f64, base: bool, srb: bool) -> Result<(), &'static str> {
        if !eta.is_finite() || eta < 0.0 {
            return Err("gCP eta must be finite and nonnegative");
        }
        if eta != self.eta {
            let slater: Vec<_> = self
                .effective
                .iter()
                .map(|&element| {
                    let scale = if self.eta_spec > 0.0 && element + 1 > 10 {
                        eta * self.eta_spec
                    } else {
                        eta
                    };
                    scale * value(GCP_TABLES, TABLES * TABLE_STRIDE + element)
                })
                .collect();
            if slater.iter().any(|value| !value.is_finite()) {
                return Err("gCP Slater exponents overflow");
            }
            self.slater = slater;
        }
        self.eta = eta;
        self.base = base;
        self.srb = srb;
        Ok(())
    }

    pub(crate) fn radius(&self, first: usize, second: usize, srb: bool) -> f64 {
        let radii = if srb { &self.rvdw_srb } else { &self.rvdw };
        if let Some(radii) = radii {
            return radii[first * self.numbers.len() + second];
        }
        let elements = if srb {
            [
                self.numbers[first] as usize - 1,
                self.numbers[second] as usize - 1,
            ]
        } else {
            [self.effective[first], self.effective[second]]
        };
        value(VDW_RADII, pair_index(elements[0], elements[1]))
    }

    pub(crate) fn validate_parameters(&self) -> Result<(), &'static str> {
        let count = self.numbers.len();
        if [
            self.sigma,
            self.alpha,
            self.beta,
            self.dmp_scal,
            self.dmp_exp,
            self.rscal,
            self.qscal,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
            || self.effective.len() != count
            || self.effective.iter().any(|&value| value >= ELEMENTS)
            || [&self.emiss, &self.virtuals, &self.slater]
                .iter()
                .any(|values| {
                    values.len() != count || values.iter().any(|value| !value.is_finite())
                })
            || self.slater.iter().any(|&value| value < 0.0)
        {
            return Err("invalid gCP parameters");
        }
        for radii in [&self.rvdw, &self.rvdw_srb]
            .iter()
            .filter_map(|radii| radii.as_ref())
        {
            if radii.len() != count * count
                || radii
                    .iter()
                    .any(|value| !value.is_finite() || *value <= 0.0)
                || (0..count).any(|first| {
                    (0..first).any(|second| {
                        radii[first * count + second] != radii[second * count + first]
                    })
                })
            {
                return Err("gCP radii must be finite, positive, and symmetric");
            }
        }
        Ok(())
    }
}

pub fn energy(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    param: &Gcp,
    cutoff: GcpCutoff,
    partition: WorkPartition,
) -> Result<f64, &'static str> {
    Ok(evaluate::<0>(
        numbers, positions, lattice, periodic, param, cutoff, partition,
    )?
    .energy)
}

#[allow(clippy::too_many_arguments)]
pub fn derivatives(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    param: &Gcp,
    cutoff: GcpCutoff,
    partition: WorkPartition,
) -> Result<GcpResult, &'static str> {
    evaluate::<1>(
        numbers, positions, lattice, periodic, param, cutoff, partition,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn hessian(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    param: &Gcp,
    cutoff: GcpCutoff,
    partition: WorkPartition,
) -> Result<GcpResult, &'static str> {
    evaluate::<2>(
        numbers, positions, lattice, periodic, param, cutoff, partition,
    )
}

#[allow(clippy::too_many_arguments)]
fn evaluate<const ORDER: usize>(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    param: &Gcp,
    cutoff: GcpCutoff,
    partition: WorkPartition,
) -> Result<GcpResult, &'static str> {
    let order = ORDER;
    if positions.len() != 3 * numbers.len() || numbers != param.numbers {
        return Err("gCP structure does not match its parameters");
    }
    let size = positions.len();
    let mut result = GcpResult {
        energy: 0.0,
        gradient: if order >= 1 {
            vec![0.0; size]
        } else {
            Vec::new()
        },
        virial: [0.0; 9],
        hessian: if order >= 2 {
            vec![0.0; size * size]
        } else {
            Vec::new()
        },
    };
    if param.sigma != 0.0 && param.slater.iter().all(|value| *value > 0.0) {
        accumulate::<ORDER>(
            &mut result,
            positions,
            lattice,
            periodic,
            cutoff.gcp,
            partition,
            |first, second, distance, same| {
                standard_pair::<ORDER>(param, first, second, distance, same)
            },
        );
    }
    if param.srb {
        accumulate::<ORDER>(
            &mut result,
            positions,
            lattice,
            periodic,
            cutoff.srb,
            partition,
            |first, second, distance, same| srb_pair(param, first, second, distance, same, false),
        );
    }
    if param.base {
        accumulate::<ORDER>(
            &mut result,
            positions,
            lattice,
            periodic,
            cutoff.srb,
            partition,
            |first, second, distance, same| srb_pair(param, first, second, distance, same, true),
        );
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn accumulate<const ORDER: usize>(
    result: &mut GcpResult,
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    cutoff: f64,
    partition: WorkPartition,
    radial: impl Fn(usize, usize, f64, bool) -> (f64, f64, f64) + Sync,
) {
    let order = ORDER;
    let atoms = positions.len() / 3;
    let translations = lattice_points(lattice, periodic, cutoff);
    let minimum = if order >= 2 { usize::MAX } else { 96 };
    let partials = crate::parallel::map(atoms, minimum, |start, stride| {
        let mut result = GcpResult {
            energy: 0.0,
            gradient: vec![0.0; result.gradient.len()],
            virial: [0.0; 9],
            hessian: vec![0.0; result.hessian.len()],
        };
        for first in (start..atoms).step_by(stride) {
            for second in 0..=first {
                if !partition.owns_pair(first, second) {
                    continue;
                }
                for translation in &translations {
                    let vector: [f64; 3] = std::array::from_fn(|axis| {
                        positions[3 * first + axis]
                            - positions[3 * second + axis]
                            - translation[axis]
                    });
                    let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
                    if distance2 > cutoff * cutoff || distance2 < f64::EPSILON {
                        continue;
                    }
                    let distance = distance2.sqrt();
                    let same = first == second;
                    let (energy, first_derivative, second_derivative) =
                        radial(first, second, distance, same);
                    result.energy += energy;
                    if order >= 1 {
                        let scale = first_derivative / distance;
                        for (axis, &coordinate) in vector.iter().enumerate() {
                            let component = scale * coordinate;
                            if !same {
                                result.gradient[3 * first + axis] += component;
                                result.gradient[3 * second + axis] -= component;
                            }
                            for (other, &other_coordinate) in vector.iter().enumerate() {
                                result.virial[axis + 3 * other] += component * other_coordinate;
                            }
                        }
                    }
                    if order >= 2 && !same {
                        add_pair_hessian(
                            &mut result.hessian,
                            first,
                            second,
                            vector,
                            distance,
                            first_derivative,
                            second_derivative,
                        );
                    }
                }
            }
        }
        result
    });
    for partial in partials {
        result.energy += partial.energy;
        for (total, value) in result.gradient.iter_mut().zip(partial.gradient) {
            *total += value;
        }
        for (total, value) in result.virial.iter_mut().zip(partial.virial) {
            *total += value;
        }
        for (total, value) in result.hessian.iter_mut().zip(partial.hessian) {
            *total += value;
        }
    }
}

fn standard_pair<const ORDER: usize>(
    param: &Gcp,
    first: usize,
    second: usize,
    r: f64,
    same: bool,
) -> (f64, f64, f64) {
    let first_virtual = inverse_sqrt(param.virtuals[first]);
    let second_virtual = inverse_sqrt(param.virtuals[second]);
    let mut scale =
        param.sigma * (param.emiss[first] * second_virtual + param.emiss[second] * first_virtual);
    if same {
        scale *= 0.5;
    }
    let (overlap, overlap1, overlap2) = overlap::<ORDER>(
        r,
        param.effective[first] + 1,
        param.effective[second] + 1,
        param.slater[first],
        param.slater[second],
    );
    let root = overlap.sqrt();
    let inverse = 1.0 / root;
    let inverse1 = -0.5 * overlap1 / (overlap * root);
    let inverse2 = (0.75 * overlap1 * overlap1 / overlap - 0.5 * overlap2) / (overlap * root);
    let power = r.powf(param.beta);
    let argument1 = param.alpha * param.beta * power / r;
    let argument2 = argument1 * (param.beta - 1.0) / r;
    let exponential = (-param.alpha * power).exp();
    let exponential1 = -argument1 * exponential;
    let exponential2 = (argument1 * argument1 - argument2) * exponential;
    let bsse = exponential * inverse;
    let bsse1 = exponential1 * inverse + exponential * inverse1;
    let bsse2 = exponential2 * inverse + 2.0 * exponential1 * inverse1 + exponential * inverse2;
    let (damping, damping1, damping2) = if param.damp {
        let radius = param.radius(first, second, false);
        let scaled = r / radius;
        let power = param.dmp_scal * scaled.powf(param.dmp_exp);
        let power1 = param.dmp_exp * power / r;
        let power2 = (param.dmp_exp - 1.0) * power1 / r;
        (
            power / (1.0 + power),
            power1 / (1.0 + power).powi(2),
            power2 / (1.0 + power).powi(2) - 2.0 * power1 * power1 / (1.0 + power).powi(3),
        )
    } else {
        (1.0, 0.0, 0.0)
    };
    (
        scale * bsse * damping,
        scale * (bsse1 * damping + bsse * damping1),
        scale * (bsse2 * damping + 2.0 * bsse1 * damping1 + bsse * damping2),
    )
}

fn srb_pair(
    param: &Gcp,
    first: usize,
    second: usize,
    r: f64,
    same: bool,
    base: bool,
) -> (f64, f64, f64) {
    let (first_element, second_element, radius_exponent, charge_exponent) = if base {
        (param.effective[first], param.effective[second], 0.75, 1.5)
    } else {
        (
            param.numbers[first] as usize - 1,
            param.numbers[second] as usize - 1,
            -1.0,
            0.5,
        )
    };
    let radius = param.radius(first, second, !base);
    let rate = param.rscal * radius.powf(radius_exponent);
    let charge = -(((first_element + 1) * (second_element + 1)) as f64).powf(charge_exponent);
    let mut scale = param.qscal * charge;
    if same {
        scale *= 0.5;
    }
    let exponential = (-rate * r).exp();
    (
        scale * exponential,
        -scale * rate * exponential,
        scale * rate * rate * exponential,
    )
}

fn overlap<const ORDER: usize>(
    r: f64,
    first: usize,
    second: usize,
    mut za: f64,
    mut zb: f64,
) -> (f64, f64, f64) {
    let order = ORDER;
    let first_shell = shell(first);
    let second_shell = shell(second);
    let combination = first_shell * second_shell;
    if matches!(combination, 2 | 3 | 6) && first_shell >= second_shell {
        std::mem::swap(&mut za, &mut zb);
    }
    let (m, weights, pa, qb, norm): (i32, &[f64], &[usize], &[usize], f64) = match combination {
        1 => (
            3,
            &[1.0, -1.0],
            &[2, 0],
            &[0, 2],
            0.25 * (za * zb).powf(1.5),
        ),
        2 => (
            4,
            &[1.0, -1.0, 1.0, -1.0],
            &[3, 0, 2, 1],
            &[0, 3, 1, 2],
            (1.0_f64 / 3.0).sqrt() * (za.powi(3) * zb.powi(5)).sqrt() * 0.125,
        ),
        3 => (
            5,
            &[1.0, -1.0, 2.0, -2.0],
            &[4, 0, 3, 1],
            &[0, 4, 1, 3],
            (za.powi(3) * zb.powi(7) / 7.5).sqrt() * 0.0625 / 3.0_f64.sqrt(),
        ),
        4 => (
            5,
            &[1.0, 1.0, -2.0],
            &[4, 0, 2],
            &[0, 4, 2],
            (za * zb).powf(2.5) * 0.0625 / 3.0,
        ),
        6 => (
            6,
            &[1.0, 1.0, -2.0, -2.0, 1.0, 1.0],
            &[5, 4, 3, 2, 1, 0],
            &[0, 1, 2, 3, 4, 5],
            (za.powi(5) * zb.powi(7) / 7.5).sqrt() * 0.03125 / 3.0,
        ),
        9 => (
            7,
            &[1.0, -3.0, 3.0, -1.0],
            &[6, 4, 2, 0],
            &[0, 2, 4, 6],
            (za * zb).powf(3.5) / 1440.0,
        ),
        _ => unreachable!(),
    };
    let ha = 0.5 * (za + zb);
    let hb = 0.5 * (zb - za);
    let count = m as usize + order;
    let a = aaux(ha * r, count);
    let b = if (za - zb).abs() < 0.1 {
        bint(hb * r, count)
    } else {
        baux(hb * r, count)
    };
    let mut f0 = 0.0;
    let mut f1 = 0.0;
    let mut f2 = 0.0;
    for index in 0..weights.len() {
        let (p, q, weight) = (pa[index], qb[index], weights[index]);
        f0 += weight * a[p] * b[q];
        if order >= 1 {
            f1 -= weight * (ha * a[p + 1] * b[q] + hb * a[p] * b[q + 1]);
        }
        if order >= 2 {
            f2 += weight
                * (ha * ha * a[p + 2] * b[q]
                    + 2.0 * ha * hb * a[p + 1] * b[q + 1]
                    + hb * hb * a[p] * b[q + 2]);
        }
    }
    (
        norm * r.powi(m) * f0,
        norm * (m as f64 * r.powi(m - 1) * f0 + r.powi(m) * f1),
        norm * ((m * (m - 1)) as f64 * r.powi(m - 2) * f0
            + 2.0 * m as f64 * r.powi(m - 1) * f1
            + r.powi(m) * f2),
    )
}

fn aaux(x: f64, count: usize) -> [f64; 9] {
    let mut result = [0.0; 9];
    let exponential = (-x).exp();
    result[0] = exponential / x;
    for order in 1..count {
        result[order] = (order as f64 * result[order - 1] + exponential) / x;
    }
    result
}

fn baux(x: f64, count: usize) -> [f64; 9] {
    let mut result = [0.0; 9];
    let positive_exponential = x.exp();
    let negative_exponential = (-x).exp();
    for (order, value) in result.iter_mut().take(count).enumerate() {
        let mut term = 1.0 / x;
        let mut sign = if order % 2 == 0 { 1.0 } else { -1.0 };
        let mut positive = sign * term;
        let mut negative = term;
        for index in 1..=order {
            term *= (order - index + 1) as f64 / x;
            sign = -sign;
            positive += sign * term;
            negative += term;
        }
        *value = positive_exponential * positive - negative_exponential * negative;
    }
    result
}

fn bint(x: f64, count: usize) -> [f64; 9] {
    let mut result = [0.0; 9];
    if x.abs() < 1.0e-6 {
        for order in (0..count).step_by(2) {
            result[order] = 2.0 / (order + 1) as f64;
        }
        return result;
    }
    let mut powers = [0.0; 13];
    powers[0] = 1.0;
    for index in 1..powers.len() {
        powers[index] = powers[index - 1] * -x / index as f64;
    }
    for (order, value) in result.iter_mut().take(count).enumerate() {
        for index in (order % 2..powers.len()).step_by(2) {
            *value += 2.0 * powers[index] / (order + index + 1) as f64;
        }
    }
    result
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "the Hessian is square, so its length is an exactly representable perfect square"
)]
fn add_pair_hessian(
    hessian: &mut [f64],
    first: usize,
    second: usize,
    vector: [f64; 3],
    distance: f64,
    derivative: f64,
    second_derivative: f64,
) {
    let size = (hessian.len() as f64).sqrt() as usize;
    let radial = derivative / distance;
    let cartesian = (second_derivative - radial) / (distance * distance);
    for row in 0..3 {
        for column in 0..3 {
            let block =
                cartesian * vector[row] * vector[column] + if row == column { radial } else { 0.0 };
            let first_row = 3 * first + row;
            let second_row = 3 * second + row;
            let first_column = 3 * first + column;
            let second_column = 3 * second + column;
            hessian[first_row * size + first_column] += block;
            hessian[second_row * size + second_column] += block;
            hessian[first_row * size + second_column] -= block;
            hessian[second_row * size + first_column] -= block;
        }
    }
}

fn lattice_points(lattice: Option<&[f64; 9]>, periodic: [bool; 3], cutoff: f64) -> Vec<[f64; 3]> {
    let Some(lattice) = lattice.filter(|_| periodic.iter().any(|value| *value)) else {
        return vec![[0.0; 3]];
    };
    let repetitions = lattice_repetitions(lattice, periodic, cutoff);
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

fn table_value(table: Basis, field: usize, element: usize) -> f64 {
    value(
        GCP_TABLES,
        table as usize * TABLE_STRIDE + field * ELEMENTS + element,
    )
}

fn inverse_sqrt(value: f64) -> f64 {
    if value >= 0.5 {
        value.sqrt().recip()
    } else {
        0.0
    }
}

fn shell(number: usize) -> usize {
    match number {
        1..=2 => 1,
        3..=10 => 2,
        _ => 3,
    }
}

fn effective_number(number: i32) -> Option<usize> {
    let effective = match number {
        1..=36 => number,
        37..=54 => number - 18,
        55..=57 => number - 36,
        58..=71 | 90..=94 => 21,
        72..=89 => number - 50,
        _ => return None,
    };
    Some(effective as usize - 1)
}

fn number_of_electrons(number: usize, valence: bool) -> usize {
    if valence {
        match number {
            5..=10 => number - 2,
            11..=18 => number - 10,
            _ => number,
        }
    } else {
        number
    }
}

fn method_id(method: &str) -> Method {
    match method {
        "hf" => Method::Hf,
        "dft" => Method::Dft,
        "gga" => Method::Gga,
        "b3lyp" => Method::B3lyp,
        "blyp" => Method::Blyp,
        "pbe" => Method::Pbe,
        "tpss" => Method::Tpss,
        "pw6b95" => Method::Pw6b95,
        "hf3c" => Method::Hf3c,
        "pbeh3c" => Method::Pbeh3c,
        "hse3c" => Method::Hse3c,
        "b973c" => Method::B973c,
        "b3pbe3c" => Method::B3pbe3c,
        "r2scan3c" => Method::R2scan3c,
        _ => Method::Unknown,
    }
}

fn basis_id(basis: &str) -> Option<Basis> {
    Some(match basis {
        "sv" => Basis::Sv,
        "sv(p)" | "def2sv(p)" | "sv_p" | "def2sv_p" => Basis::SvP,
        "svx" => Basis::Svx,
        "svp" => Basis::Svp,
        "minis" => Basis::Minis,
        "631gd" | "631gs" => Basis::G631,
        "tz" | "def2tzvp" => Basis::Tz,
        "deftzvp" | "def1tzvp" => Basis::Def1tzvp,
        "ccdz" | "ccpvdz" => Basis::Ccdz,
        "accdz" | "augccpvdz" | "accpvdz" => Basis::Accdz,
        "pobtz" | "pobtzvp" => Basis::Pobtz,
        "minix" | "hf3c" => Basis::Minix,
        "gcore" => Basis::Gcore,
        "2g" | "twog" | "fitg" => Basis::TwoG,
        "dzp" => Basis::Dzp,
        "dz" => Basis::Dz,
        "msvp" | "def2msvp" => Basis::Msvp,
        "lanl" => Basis::Lanl,
        "pbeh3c" | "hse3c" => Basis::Pbeh3c,
        "mtzvp" | "def2mtzvp" => Basis::Def2mtzvp,
        "mtzvpp" | "def2mtzvpp" | "r2scan3c" => Basis::Def2mtzvpp,
        _ => return None,
    })
}

fn is_dft(method: Method) -> bool {
    matches!(
        method,
        Method::Dft | Method::B3lyp | Method::Blyp | Method::Gga | Method::Tpss | Method::Pw6b95
    )
}

fn is_hybrid(method: Method) -> bool {
    matches!(method, Method::B3lyp | Method::Pw6b95)
}

fn is_gga(method: Method) -> bool {
    matches!(method, Method::Gga | Method::Tpss | Method::Blyp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn controls_preserve_scaling_and_reject_invalid_eta() {
        let mut param = load(&[6, 17], Some("r2scan3c"), None).unwrap();
        let original = param.slater.clone();
        let (eta, base, srb) = param.controls();
        param.set_controls(2.0 * eta, !base, !srb).unwrap();
        assert_eq!(
            param.slater,
            original.iter().map(|value| 2.0 * value).collect::<Vec<_>>()
        );
        assert_eq!(param.controls(), (2.0 * eta, !base, !srb));
        param.set_controls(eta, base, srb).unwrap();
        assert_eq!(param.slater, original);
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(param.set_controls(bad, !base, !srb).is_err());
            assert_eq!(param.controls(), (eta, base, srb));
            assert_eq!(param.slater, original);
        }
        param.set_controls(0.0, false, false).unwrap();
        assert!(param.slater.iter().all(|&value| value == 0.0));
    }

    #[test]
    fn matches_reference_families() {
        let numbers = [6, 8, 7];
        let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 4.0, 0.0];
        let cases = [
            (
                "hf",
                Some("minix"),
                0.0013746606180596486,
                0.0005738507309585009,
            ),
            ("pbeh3c", None, 0.0013367957792382655, 0.0006857278480993941),
            ("hse3c", None, 0.001279680707661919, 0.0006580382057369658),
            (
                "b973c",
                None,
                -0.00018565968877661711,
                -0.00011603708836002195,
            ),
            ("hf3c", None, 0.0007812866164643338, 7.292761771734242e-5),
            (
                "r2scan3c",
                None,
                0.00044649249961611307,
                0.0001714116831929529,
            ),
            (
                "hf",
                Some("tz"),
                0.0001482559736207828,
                3.1186208539109405e-5,
            ),
            (
                "b3lyp",
                Some("631gd"),
                0.0011378185466096387,
                0.00024897565392941997,
            ),
        ];
        for (method, basis, expected_energy, expected_gradient) in cases {
            let param = load(&numbers, Some(method), basis).unwrap();
            let result = hessian(
                &numbers,
                &positions,
                None,
                [false; 3],
                &param,
                GcpCutoff::default(),
                WorkPartition::SERIAL,
            )
            .unwrap();
            assert!(
                (result.energy - expected_energy).abs() < 1.0e-13,
                "{method}: {}",
                result.energy
            );
            assert!(
                (result.gradient[0] - expected_gradient).abs() < 5.0e-13,
                "{method}: {}",
                result.gradient[0]
            );
        }
    }

    #[test]
    fn matches_derivative_hessian_and_cutoff_references() {
        let numbers = [6, 8, 7];
        let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 4.0, 0.0];
        let param = load(&numbers, Some("pbeh3c"), None).unwrap();
        let result = hessian(
            &numbers,
            &positions,
            None,
            [false; 3],
            &param,
            GcpCutoff::default(),
            WorkPartition::SERIAL,
        )
        .unwrap();
        assert!((result.energy - 0.0013367957792382655).abs() < 1.0e-13);
        assert!((result.gradient[0] - 0.0006857278480993941).abs() < 5.0e-13);
        assert!((result.gradient[1] - 0.0003944996529004371).abs() < 5.0e-13);
        assert!((result.virial[0] - -0.004347551721006959).abs() < 5.0e-13);
        assert!((result.virial[4] - -0.0028914107450121744).abs() < 5.0e-13);
        assert!((result.hessian[0] - 0.0004491912463559838).abs() < 5.0e-12);
        assert!((result.hessian[1] - -8.858055566709091e-5).abs() < 5.0e-12);
        assert!((result.hessian[3] - -0.0005699612984978659).abs() < 5.0e-12);

        let cutoff = hessian(
            &numbers,
            &positions,
            None,
            [false; 3],
            &param,
            GcpCutoff { gcp: 4.5, srb: 4.5 },
            WorkPartition::SERIAL,
        )
        .unwrap();
        assert!((cutoff.energy - 0.0006577241271364464).abs() < 1.0e-13);
        assert!((cutoff.gradient[0] - 9.862491322510927e-5).abs() < 5.0e-13);
        assert!((cutoff.virial[4] - -0.0015779986116017484).abs() < 5.0e-13);
        assert!((cutoff.hessian[0] - -0.0001207700521418821).abs() < 5.0e-12);
    }

    #[test]
    fn matches_periodic_reference_and_partitions() {
        let numbers = [6, 8, 7];
        let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 4.0, 0.0];
        let lattice = [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0];
        let param = load(&numbers, Some("pbeh3c"), None).unwrap();
        let serial = hessian(
            &numbers,
            &positions,
            Some(&lattice),
            [true; 3],
            &param,
            GcpCutoff::default(),
            WorkPartition::SERIAL,
        )
        .unwrap();
        assert!((serial.energy - 0.001805303821670533).abs() < 1.0e-13);
        assert!((serial.gradient[0] - 0.00011199508985802613).abs() < 5.0e-13);
        assert!((serial.gradient[1] - 0.0003142674657973423).abs() < 5.0e-13);
        assert!((serial.virial[0] - -0.007487015291475585).abs() < 5.0e-12);
        assert!((serial.hessian[0] - 0.0010101320772209767).abs() < 5.0e-12);
        assert!((serial.hessian[1] - -0.00011465243638749777).abs() < 5.0e-12);
        assert!((serial.hessian[3] - -0.0011399227050866401).abs() < 5.0e-12);
        let parts: Vec<_> = (0..2)
            .map(|part| {
                hessian(
                    &numbers,
                    &positions,
                    Some(&lattice),
                    [true; 3],
                    &param,
                    GcpCutoff::default(),
                    WorkPartition::new(part, 2).unwrap(),
                )
                .unwrap()
            })
            .collect();
        assert!(
            (parts.iter().map(|part| part.energy).sum::<f64>() - serial.energy).abs() < 1.0e-14
        );
        for index in 0..serial.gradient.len() {
            assert!(
                (parts.iter().map(|part| part.gradient[index]).sum::<f64>()
                    - serial.gradient[index])
                    .abs()
                    < 1.0e-13
            );
        }
        for index in 0..serial.hessian.len() {
            assert!(
                (parts.iter().map(|part| part.hessian[index]).sum::<f64>() - serial.hessian[index])
                    .abs()
                    < 1.0e-12
            );
        }
    }

    #[test]
    fn threaded_pairs_match_serial_hessian_path() {
        let numbers: Vec<_> = (0..96).map(|index| [6, 8, 7, 1][index % 4]).collect();
        let positions: Vec<_> = (0..96)
            .flat_map(|index| {
                [
                    3.1 * (index % 5) as f64,
                    3.2 * ((index / 5) % 5) as f64,
                    3.3 * (index / 25) as f64,
                ]
            })
            .collect();
        let lattice = [20.0, 0.0, 0.0, 0.7, 21.0, 0.0, 0.3, 0.5, 22.0];
        let cutoff = GcpCutoff {
            gcp: 12.0,
            srb: 12.0,
        };
        for periodic in [[false; 3], [true; 3]] {
            for method in ["pbeh3c", "hf3c", "b973c"] {
                let param = load(&numbers, Some(method), None).unwrap();
                let reference = hessian(
                    &numbers,
                    &positions,
                    Some(&lattice),
                    periodic,
                    &param,
                    cutoff,
                    WorkPartition::SERIAL,
                )
                .unwrap();
                let mut total = vec![0.0; 1 + positions.len() + 9];
                for part in 0..3 {
                    let partition = WorkPartition::new(part, 3).unwrap();
                    let actual = derivatives(
                        &numbers,
                        &positions,
                        Some(&lattice),
                        periodic,
                        &param,
                        cutoff,
                        partition,
                    )
                    .unwrap();
                    let scalar = energy(
                        &numbers,
                        &positions,
                        Some(&lattice),
                        periodic,
                        &param,
                        cutoff,
                        partition,
                    )
                    .unwrap();
                    assert!((scalar - actual.energy).abs() < 1e-12);
                    for (total, value) in total.iter_mut().zip(
                        std::iter::once(actual.energy)
                            .chain(actual.gradient)
                            .chain(actual.virial),
                    ) {
                        *total += value;
                    }
                }
                for (actual, expected) in total.into_iter().zip(
                    std::iter::once(reference.energy)
                        .chain(reference.gradient)
                        .chain(reference.virial),
                ) {
                    assert!((actual - expected).abs() < 1e-11);
                }
            }
        }
    }

    #[test]
    fn matches_self_partial_and_skew_periodic_references() {
        let lattice = [8.0, 0.0, 0.0, 1.0, 9.0, 0.0, 0.5, 0.7, 10.0];
        let cutoff = GcpCutoff {
            gcp: 12.0,
            srb: 12.0,
        };
        let self_param = load(&[6], Some("pbeh3c"), None).unwrap();
        let self_result = derivatives(
            &[6],
            &[0.4, 0.8, 1.2],
            Some(&lattice),
            [true; 3],
            &self_param,
            cutoff,
            WorkPartition::SERIAL,
        )
        .unwrap();
        assert!((self_result.energy - 1.5432984867825415e-7).abs() < 1.0e-16);
        assert!(self_result
            .gradient
            .iter()
            .all(|value| value.abs() < 1.0e-20));
        assert!((self_result.virial[0] + 3.6585164221279463e-6).abs() < 1.0e-15);

        let numbers = [6, 8];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
        let param = load(&numbers, Some("pbeh3c"), None).unwrap();
        let partial = derivatives(
            &numbers,
            &positions,
            Some(&lattice),
            [true, false, false],
            &param,
            cutoff,
            WorkPartition::SERIAL,
        )
        .unwrap();
        assert!((partial.energy - 0.0007326131208736937).abs() < 1.0e-13);
        assert!((partial.gradient[1] - 0.00042872230585806807).abs() < 1.0e-13);
        assert!((partial.virial[0] + 0.003484094274295485).abs() < 1.0e-12);

        let skew = derivatives(
            &numbers,
            &positions,
            Some(&lattice),
            [true; 3],
            &param,
            cutoff,
            WorkPartition::SERIAL,
        )
        .unwrap();
        assert!((skew.energy - 0.000733389638788683).abs() < 1.0e-13);
        assert!((skew.gradient[0] - 6.214839320259188e-7).abs() < 1.0e-13);
        assert!((skew.virial[0] + 0.003486778083367054).abs() < 1.0e-12);
    }

    #[test]
    fn derivatives_match_finite_differences() {
        let numbers = [6, 8, 7];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6, 2.1, 5.0, 1.7];
        let lattice = [8.0, 0.0, 0.0, 1.0, 9.0, 0.0, 0.5, 0.7, 10.0];
        let cutoff = GcpCutoff {
            gcp: 12.0,
            srb: 12.0,
        };
        for method in ["pbeh3c", "hf3c", "b973c"] {
            let param = load(&numbers, Some(method), None).unwrap();
            let result = derivatives(
                &numbers,
                &positions,
                Some(&lattice),
                [true; 3],
                &param,
                cutoff,
                WorkPartition::SERIAL,
            )
            .unwrap();
            let step = 1.0e-5;
            for coordinate in 0..positions.len() {
                let mut plus = positions;
                let mut minus = positions;
                plus[coordinate] += step;
                minus[coordinate] -= step;
                let finite_difference = (energy(
                    &numbers,
                    &plus,
                    Some(&lattice),
                    [true; 3],
                    &param,
                    cutoff,
                    WorkPartition::SERIAL,
                )
                .unwrap()
                    - energy(
                        &numbers,
                        &minus,
                        Some(&lattice),
                        [true; 3],
                        &param,
                        cutoff,
                        WorkPartition::SERIAL,
                    )
                    .unwrap())
                    / (2.0 * step);
                assert!(
                    (finite_difference - result.gradient[coordinate]).abs() < 1.0e-9,
                    "{method} coordinate {coordinate}: {finite_difference} != {}",
                    result.gradient[coordinate]
                );
            }
            for axis in 0..3 {
                for other in 0..3 {
                    let mut plus_positions = positions;
                    let mut minus_positions = positions;
                    let mut plus_lattice = lattice;
                    let mut minus_lattice = lattice;
                    for atom in 0..numbers.len() {
                        plus_positions[3 * atom + axis] += step * positions[3 * atom + other];
                        minus_positions[3 * atom + axis] -= step * positions[3 * atom + other];
                    }
                    for column in 0..3 {
                        plus_lattice[axis + 3 * column] += step * lattice[other + 3 * column];
                        minus_lattice[axis + 3 * column] -= step * lattice[other + 3 * column];
                    }
                    let finite_difference = (energy(
                        &numbers,
                        &plus_positions,
                        Some(&plus_lattice),
                        [true; 3],
                        &param,
                        cutoff,
                        WorkPartition::SERIAL,
                    )
                    .unwrap()
                        - energy(
                            &numbers,
                            &minus_positions,
                            Some(&minus_lattice),
                            [true; 3],
                            &param,
                            cutoff,
                            WorkPartition::SERIAL,
                        )
                        .unwrap())
                        / (2.0 * step);
                    assert!(
                        (finite_difference - result.virial[axis + 3 * other]).abs() < 1.0e-9,
                        "{method} strain {axis},{other}: {finite_difference} != {}",
                        result.virial[axis + 3 * other]
                    );
                }
            }
        }
    }

    #[test]
    fn analytical_hessian_matches_gradient_differences() {
        let numbers = [6, 8, 7];
        let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6, 2.1, 5.0, 1.7];
        let lattice = [8.0, 0.0, 0.0, 1.0, 9.0, 0.0, 0.5, 0.7, 10.0];
        let cutoff = GcpCutoff {
            gcp: 12.0,
            srb: 12.0,
        };
        for (method, periodic) in [
            ("pbeh3c", [false; 3]),
            ("pbeh3c", [true; 3]),
            ("hf3c", [false; 3]),
            ("b973c", [false; 3]),
            ("r2scan3c", [false; 3]),
        ] {
            let param = load(&numbers, Some(method), None).unwrap();
            let cell = periodic.iter().any(|value| *value).then_some(&lattice);
            let result = hessian(
                &numbers,
                &positions,
                cell,
                periodic,
                &param,
                cutoff,
                WorkPartition::SERIAL,
            )
            .unwrap();
            let size = positions.len();
            let step = 1.0e-5;
            for coordinate in 0..size {
                let mut plus = positions;
                let mut minus = positions;
                plus[coordinate] += step;
                minus[coordinate] -= step;
                let plus = derivatives(
                    &numbers,
                    &plus,
                    cell,
                    periodic,
                    &param,
                    cutoff,
                    WorkPartition::SERIAL,
                )
                .unwrap();
                let minus = derivatives(
                    &numbers,
                    &minus,
                    cell,
                    periodic,
                    &param,
                    cutoff,
                    WorkPartition::SERIAL,
                )
                .unwrap();
                for row in 0..size {
                    let finite_difference =
                        (plus.gradient[row] - minus.gradient[row]) / (2.0 * step);
                    assert!(
                        (finite_difference - result.hessian[row * size + coordinate]).abs()
                            < 2.0e-8,
                        "{method} periodic={periodic:?} ({row},{coordinate}): {finite_difference} != {}",
                        result.hessian[row * size + coordinate]
                    );
                }
            }
            for row in 0..size {
                for column in 0..size {
                    assert!(
                        (result.hessian[row * size + column] - result.hessian[column * size + row])
                            .abs()
                            < 1.0e-14
                    );
                }
                for axis in 0..3 {
                    let sum = (axis..size)
                        .step_by(3)
                        .map(|column| result.hessian[row * size + column])
                        .sum::<f64>();
                    assert!(sum.abs() < 1.0e-14, "{method} row {row} axis {axis}: {sum}");
                }
            }
        }
    }
}
