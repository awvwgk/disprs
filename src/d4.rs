use std::cell::UnsafeCell;
use std::convert::TryInto;
use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

pub use crate::d3::EwaldConfig;
use crate::geometry::{
    determinant, inverse, lattice_repetitions, lattice_translation, periodic_reciprocal,
    smooth_cutoff, squared_distance, wrap_positions,
};
use crate::parameters::{lookup, parameter};
pub use crate::WorkPartition;

const ELEMENTS: usize = 118;
const REFERENCES: usize = 7;
const FREQUENCIES: usize = 23;
const EEQBC_ELEMENTS: usize = 103;
const NREF: &[u8; 4 * ELEMENTS] = include_bytes!("../assets/d4_nref.bin");
const NGW: &[u8; 4 * REFERENCES * ELEMENTS] = include_bytes!("../assets/d4_ngw.bin");
const ELEMENT_DATA: &[u8; 8 * 5 * ELEMENTS] = include_bytes!("../assets/d4_element.bin");
const EEQ_DATA: &[u8; 8 * 4 * ELEMENTS] = include_bytes!("../assets/d4_eeq.bin");
const REFCN: &[u8; 8 * REFERENCES * ELEMENTS] = include_bytes!("../assets/d4_refcn.bin");
const REFQ: &[u8; 8 * REFERENCES * ELEMENTS] = include_bytes!("../assets/d4_refq.bin");
const REFALPHA: &[u8; 8 * FREQUENCIES * REFERENCES * ELEMENTS] =
    include_bytes!("../assets/d4_refalpha.bin");
const CHARGE_RCOV: &[u8; 8 * ELEMENTS] = include_bytes!("../assets/d4_charge_rcov.bin");
const PAIR_WEIGHTS: &[u8; 8 * ELEMENTS * ELEMENTS] =
    include_bytes!("../assets/d4_pair_weights.bin");
const REFERENCE_C6: &[u8; 8 * (REFERENCES * ELEMENTS).pow(2)] =
    include_bytes!("../assets/d4_reference_c6.bin");
const RAW_ALPHA: &[u8; 8 * FREQUENCIES * REFERENCES * ELEMENTS] =
    include_bytes!("../assets/d4_raw_alpha.bin");
const SECONDARY: &[u8; 8 * FREQUENCIES * REFERENCES * ELEMENTS] =
    include_bytes!("../assets/d4_secondary.bin");
const RAW_PARAMETERS: &[u8; 8 * 5 * REFERENCES * ELEMENTS] =
    include_bytes!("../assets/d4_raw_parameters.bin");
const QUADRATURE: &[u8; 8 * FREQUENCIES] = include_bytes!("../assets/d4_quadrature.bin");
const EEQBC_REFQ: &[u8; 8 * REFERENCES * ELEMENTS] = include_bytes!("../assets/d4_eeqbc_refq.bin");
const EEQBC_REFH: &[u8; 8 * REFERENCES * ELEMENTS] = include_bytes!("../assets/d4_eeqbc_refh.bin");
const EEQBC_DATA: &[u8; 8 * 10 * EEQBC_ELEMENTS] = include_bytes!("../assets/d4_eeqbc.bin");
const EEQBC_VDW: &[u8; 8 * EEQBC_ELEMENTS * EEQBC_ELEMENTS] =
    include_bytes!("../assets/d4_eeqbc_vdw.bin");

#[derive(Clone)]
pub struct Properties {
    pub coordination: Vec<f64>,
    pub charges: Vec<f64>,
    pub c6: Vec<f64>,
    pub polarizabilities: Vec<f64>,
}

#[derive(Clone, Copy, PartialEq)]
pub struct Param {
    pub s6: f64,
    pub s8: f64,
    pub s9: f64,
    pub a1: f64,
    pub a2: f64,
    pub alpha: f64,
}

#[derive(Clone)]
pub struct Dispersion {
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub virial: [f64; 9],
}

#[derive(Clone, Copy, PartialEq)]
pub struct Model<'a> {
    d4s: bool,
    eeqbc: bool,
    ga: f64,
    gc: f64,
    wf: f64,
    cutoff: Cutoff,
    charge_cutoff: f64,
    fixed_charges: Option<&'a [f64]>,
    ghosts: &'a [bool],
    partition: WorkPartition,
    pub(crate) ewald: Option<EwaldConfig>,
}

impl<'a> Model<'a> {
    pub const D4: Self = Self {
        d4s: false,
        eeqbc: false,
        ga: 3.0,
        gc: 2.0,
        wf: 6.0,
        charge_cutoff: 60.0,
        fixed_charges: None,
        ghosts: &[],
        partition: WorkPartition::SERIAL,
        ewald: None,
        cutoff: Cutoff {
            cn: 30.0,
            disp2: 60.0,
            disp3: 40.0,
            width2: 0.0,
            width3: 0.0,
        },
    };
    pub const D4S: Self = Self {
        d4s: true,
        ..Self::D4
    };

    pub fn custom(d4s: bool, ga: f64, gc: f64, wf: f64) -> Result<Self, &'static str> {
        if [ga, gc, wf]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err("D4 model parameters must be finite and positive");
        }
        if d4s && wf != 6.0 {
            return Err(
                "D4S uses pair-specific weighting widths; scalar wf is only supported for D4",
            );
        }
        Ok(Self {
            ga,
            gc,
            wf,
            ..if d4s { Self::D4S } else { Self::D4 }
        })
    }

    pub fn set_charge_cutoff(&mut self, cutoff: f64) -> Result<(), &'static str> {
        if !cutoff.is_finite() || cutoff <= 0.0 {
            return Err("charge summation range must be finite and positive");
        }
        self.charge_cutoff = cutoff;
        Ok(())
    }

    pub fn set_ewald(&mut self, config: Option<EwaldConfig>) -> Result<(), &'static str> {
        let config = config
            .map(|mut config| {
                if !config.tolerance.is_finite() || !config.kcut.is_finite() {
                    return Err("Ewald controls must be finite");
                }
                if config.mesh > 0
                    && !(config.mesh as usize)
                        .checked_next_power_of_two()
                        .and_then(|size| size.checked_pow(3))
                        .and_then(|count| count.checked_mul(std::mem::size_of::<[f64; 2]>()))
                        .is_some_and(|bytes| bytes <= isize::MAX as usize)
                {
                    return Err("Ewald mesh size is not representable");
                }
                config.tolerance = if config.tolerance > 0.0 {
                    config.tolerance
                } else {
                    1e-4
                };
                config.kcut = config.kcut.max(0.0);
                Ok(config)
            })
            .transpose()?;
        self.ewald = config;
        Ok(())
    }

    pub fn set_charge_model(&mut self, charge_model: i32) -> Result<(), &'static str> {
        self.eeqbc = match charge_model {
            0 => false,
            1 => true,
            _ => return Err("D4 charge model must be 0 (EEQ) or 1 (EEQBC)"),
        };
        Ok(())
    }

    pub fn set_cutoff(&mut self, cutoff: Cutoff) -> Result<(), &'static str> {
        validate_cutoff(cutoff)?;
        self.cutoff = cutoff;
        Ok(())
    }

    pub fn with_ghosts<'b>(self, ghosts: &'b [bool]) -> Model<'b>
    where
        'a: 'b,
    {
        Model { ghosts, ..self }
    }

    /// Use geometry-independent atomic charges without renormalization.
    /// None restores the selected EEQ/EEQBC solver; reference tables are unchanged.
    pub fn with_fixed_charges<'b>(
        self,
        charges: Option<&'b [f64]>,
    ) -> Result<Model<'b>, &'static str>
    where
        'a: 'b,
    {
        if charges.is_some_and(|values| values.iter().any(|value| !value.is_finite())) {
            return Err("fixed D4 charges must be finite");
        }
        Ok(Model {
            fixed_charges: charges,
            ..self
        })
    }

    pub fn set_work_partition(&mut self, part: i32, parts: i32) -> Result<(), &'static str> {
        self.partition = WorkPartition::new(part, parts).ok_or("invalid work partition")?;
        Ok(())
    }

    fn validate_selection(self, atoms: usize) -> Result<(), &'static str> {
        if self
            .fixed_charges
            .is_some_and(|charges| charges.len() != atoms)
        {
            return Err("fixed D4 charges must contain one value per atom");
        }
        if !self.ghosts.is_empty() && self.ghosts.len() != atoms {
            return Err("ghost mask must contain one entry per atom");
        }
        Ok(())
    }

    fn active(self, atom: usize) -> bool {
        self.ghosts.get(atom) != Some(&true)
    }

    fn owns_pair(self, first: usize, second: usize) -> bool {
        self.active(first) && self.active(second) && self.partition.owns_pair(first, second)
    }
}

#[derive(Clone, Copy, PartialEq)]
pub struct Cutoff {
    pub cn: f64,
    pub disp2: f64,
    pub disp3: f64,
    pub width2: f64,
    pub width3: f64,
}

impl Default for Cutoff {
    fn default() -> Self {
        Model::D4.cutoff
    }
}

pub fn load_param(method: &str, atm: bool) -> Option<Param> {
    let method = method.split('/').next()?;
    let line = if atm {
        lookup(method, "d4.bj")
    } else {
        lookup(method, "d4.bj-eeq-two").or_else(|| lookup(method, "d4.bj"))
    }?;
    Some(Param {
        s6: parameter(line, "s6", 1.0),
        s8: parameter(line, "s8", 0.0),
        s9: if atm { parameter(line, "s9", 1.0) } else { 0.0 },
        a1: parameter(line, "a1", 0.0),
        a2: parameter(line, "a2", 0.0),
        alpha: parameter(line, "alp", 16.0),
    })
}

pub fn energy(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
) -> Result<f64, &'static str> {
    energy_with_cutoff(numbers, positions, charge, model, param, model.cutoff)
}

pub fn energy_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    cutoff: Cutoff,
) -> Result<f64, &'static str> {
    if model.ewald.is_some() {
        return Err("D4 Fourier dispersion requires full 3D periodicity");
    }
    validate_cutoff(cutoff)?;
    let properties = properties_with_cutoff(numbers, positions, charge, model, cutoff.cn)?;
    let mut result = two_body_energy(numbers, positions, &properties.c6, param, cutoff, model);
    if param.s9.abs() >= f64::EPSILON {
        let c6 = model_coefficients(
            numbers,
            &properties.coordination,
            &vec![0.0; numbers.len()],
            model,
        );
        result += atm_energy(numbers, positions, &c6, param, cutoff, model);
    }
    Ok(result)
}

pub fn dispersion(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
) -> Result<Dispersion, &'static str> {
    dispersion_with_cutoff(numbers, positions, charge, model, param, model.cutoff)
}

pub fn dispersion_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    cutoff: Cutoff,
) -> Result<Dispersion, &'static str> {
    if model.ewald.is_some() {
        return Err("D4 Fourier dispersion requires full 3D periodicity");
    }
    let (value, gradient) =
        differentiated_energy(numbers, positions, charge, model, param, cutoff)?;
    let mut virial = [0.0; 9];
    for atom in 0..numbers.len() {
        for first in 0..3 {
            for second in 0..3 {
                virial[first + 3 * second] +=
                    gradient[3 * atom + first] * positions[3 * atom + second];
            }
        }
    }
    Ok(Dispersion {
        energy: value,
        gradient,
        virial,
    })
}

pub fn periodic_energy(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<f64, &'static str> {
    if let Some(config) = model.ewald {
        return fourier::dispersion(
            numbers, positions, charge, model, param, lattice, periodic, config,
        )
        .map(|result| result.energy);
    }
    if !periodic.iter().any(|&active| active) {
        return energy(numbers, positions, charge, model, param);
    }
    let positions = wrap_positions(positions, lattice, periodic)?;
    periodic_scalar_energy::<false>(numbers, &positions, charge, model, param, lattice, periodic)
        .map(|result| result.0)
}

pub fn periodic_pairwise(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    if model.ewald.is_some() {
        return Err("D4 Fourier pair matrices are not supported");
    }
    if !periodic.iter().any(|&active| active) {
        return pairwise(numbers, positions, charge, model, param);
    }
    let positions = wrap_positions(positions, lattice, periodic)?;
    periodic_scalar_energy::<true>(numbers, &positions, charge, model, param, lattice, periodic)
        .map(|(_, pair2, pair3)| (pair2, pair3))
}

pub fn periodic_dispersion(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<Dispersion, &'static str> {
    if let Some(config) = model.ewald {
        return fourier::dispersion(
            numbers, positions, charge, model, param, lattice, periodic, config,
        );
    }
    if !periodic.iter().any(|&active| active) {
        return dispersion(numbers, positions, charge, model, param);
    }
    let wrapped = wrap_positions(positions, lattice, periodic)?;
    let positions = wrapped.as_slice();
    let result = periodic_differentiated_energy(
        numbers, positions, charge, model, param, lattice, periodic,
    )?;
    let coordinates = positions.len();
    let derivatives = result.gradient();
    let gradient = derivatives[..coordinates].to_vec();
    let mut virial = [0.0; 9];
    for first in 0..3 {
        for second in 0..3 {
            virial[first + 3 * second] = (0..3)
                .map(|column| {
                    derivatives[coordinates + first + 3 * column] * lattice[second + 3 * column]
                })
                .sum::<f64>()
                + (0..numbers.len())
                    .map(|atom| gradient[3 * atom + first] * positions[3 * atom + second])
                    .sum::<f64>();
        }
    }
    Ok(Dispersion {
        energy: result.value,
        gradient,
        virial,
    })
}

pub fn periodic_properties(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<Properties, &'static str> {
    model.validate_selection(numbers.len())?;
    if !periodic.iter().any(|&active| active) {
        return properties(numbers, positions, charge, model);
    }
    let wrapped = wrap_positions(positions, lattice, periodic)?;
    let positions = wrapped.as_slice();
    let (coordination, charges) =
        periodic_environment(numbers, positions, charge, lattice, periodic, model)?;
    let c6 = model_coefficients(numbers, &coordination, &charges, model);
    let weights: Vec<_> = numbers
        .iter()
        .zip(&coordination)
        .zip(&charges)
        .map(|((&number, &cn), &q)| weights(number as usize - 1, number as usize - 1, cn, q, model))
        .collect();
    let polarizabilities = numbers
        .iter()
        .enumerate()
        .map(|(atom, &number)| {
            (0..int(NREF, number as usize - 1) as usize)
                .map(|reference| {
                    weights[atom][reference] * alpha(number as usize - 1, reference, 0, model)
                })
                .sum()
        })
        .collect();
    Ok(Properties {
        coordination,
        charges,
        c6,
        polarizabilities,
    })
}

pub fn pairwise(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    pairwise_with_cutoff(numbers, positions, charge, model, param, model.cutoff)
}

pub fn pairwise_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    cutoff: Cutoff,
) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    if model.ewald.is_some() {
        return Err("D4 Fourier pair matrices are not supported");
    }
    validate_cutoff(cutoff)?;
    let properties = properties_with_cutoff(numbers, positions, charge, model, cutoff.cn)?;
    let atoms = numbers.len();
    let mut pair2 = vec![0.0; atoms * atoms];
    for first in 0..atoms {
        for second in 0..first {
            if !model.owns_pair(first, second) {
                continue;
            }
            if squared_distance(positions, first, second) > cutoff.disp2 * cutoff.disp2 {
                continue;
            }
            let contribution = pair_energy(
                numbers[first] as usize - 1,
                numbers[second] as usize - 1,
                squared_distance(positions, first, second),
                properties.c6[first * atoms + second],
                param,
            );
            let contribution = contribution
                * smooth_cutoff(
                    squared_distance(positions, first, second),
                    cutoff.disp2,
                    cutoff.width2,
                )
                .0;
            pair2[first * atoms + second] = 0.5 * contribution;
            pair2[second * atoms + first] = 0.5 * contribution;
        }
    }
    let mut pair3 = vec![0.0; atoms * atoms];
    if param.s9.abs() >= f64::EPSILON {
        let pairs = atm_pair_data(numbers, positions, param, cutoff);
        let c6 = model_coefficients(
            numbers,
            &properties.coordination,
            &vec![0.0; numbers.len()],
            model,
        );
        for first in 0..atoms {
            for second in 0..first {
                if !model.owns_pair(first, second) {
                    continue;
                }
                for third in 0..second {
                    if !model.active(third) {
                        continue;
                    }
                    let contribution =
                        atm_triplet(atoms, &pairs, &c6, param, [first, second, third]);
                    for (left, right) in [(first, second), (first, third), (second, third)] {
                        pair3[left * atoms + right] += contribution / 6.0;
                        pair3[right * atoms + left] += contribution / 6.0;
                    }
                }
            }
        }
    }
    Ok((pair2, pair3))
}

fn int(table: &[u8], index: usize) -> i32 {
    let start = 4 * index;
    i32::from_le_bytes(table[start..start + 4].try_into().unwrap())
}

fn real(table: &[u8], index: usize) -> f64 {
    let start = 8 * index;
    f64::from_le_bytes(table[start..start + 8].try_into().unwrap())
}

pub(crate) fn element(element: usize, field: usize) -> f64 {
    real(ELEMENT_DATA, field + 5 * element)
}

pub(crate) fn quadrature(frequency: usize) -> f64 {
    real(QUADRATURE, frequency)
}

fn eeq(element: usize, field: usize) -> f64 {
    real(EEQ_DATA, field + 4 * element)
}

fn reference(table: &[u8], element: usize, reference: usize) -> f64 {
    real(table, reference + REFERENCES * element)
}

fn alpha(element: usize, reference: usize, frequency: usize, model: Model) -> f64 {
    if model.eeqbc || model.ga != 3.0 || model.gc != 2.0 {
        return custom_alpha(element, reference, model)[frequency];
    }
    real(
        REFALPHA,
        frequency + FREQUENCIES * (reference + REFERENCES * element),
    )
}

fn custom_alpha(element: usize, reference: usize, model: Model) -> [f64; FREQUENCIES] {
    let index = reference + REFERENCES * element;
    let scale = real(RAW_PARAMETERS, 5 * index);
    let count = real(RAW_PARAMETERS, 5 * index + 1);
    let effective = real(RAW_PARAMETERS, 5 * index + 2);
    let hardness = real(RAW_PARAMETERS, 5 * index + 3);
    let charge = if model.eeqbc {
        real(EEQBC_REFH, index) + effective
    } else {
        real(RAW_PARAMETERS, 5 * index + 4)
    };
    if scale == 0.0 {
        return [0.0; FREQUENCIES];
    }
    let zeta = if charge < 0.0 {
        model.ga.exp()
    } else {
        (model.ga * (1.0 - (model.gc * hardness * (1.0 - effective / charge)).exp())).exp()
    };
    std::array::from_fn(|frequency| {
        let index = frequency + FREQUENCIES * index;
        (scale * (real(RAW_ALPHA, index) - count * (real(SECONDARY, index) * zeta))).max(0.0)
    })
}

fn reference_c6(
    first_element: usize,
    first_reference: usize,
    second_element: usize,
    second_reference: usize,
    model: Model,
) -> f64 {
    if model.eeqbc || model.ga != 3.0 || model.gc != 2.0 {
        return 3.0 / PI
            * (0..FREQUENCIES)
                .map(|frequency| {
                    alpha(first_element, first_reference, frequency, model)
                        * alpha(second_element, second_reference, frequency, model)
                        * quadrature(frequency)
                })
                .sum::<f64>();
    }
    real(
        REFERENCE_C6,
        first_reference
            + REFERENCES
                * (first_element + ELEMENTS * (second_reference + REFERENCES * second_element)),
    )
}

unsafe extern "C" {
    fn erf(value: f64) -> f64;
}

pub fn properties(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
) -> Result<Properties, &'static str> {
    properties_with_cutoff(numbers, positions, charge, model, model.cutoff.cn)
}

fn properties_with_cutoff(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    cutoff: f64,
) -> Result<Properties, &'static str> {
    validate(numbers, positions)?;
    model.validate_selection(numbers.len())?;
    let coordination = coordination_numbers(numbers, positions, cutoff);
    let charges = charges(numbers, positions, charge, model)?;
    let weights: Vec<_> = numbers
        .iter()
        .zip(coordination.iter())
        .zip(charges.iter())
        .map(|((&number, &cn), &charge)| {
            weights(number as usize - 1, number as usize - 1, cn, charge, model)
        })
        .collect();
    let c6 = model_coefficients(numbers, &coordination, &charges, model);
    let mut polarizabilities = vec![0.0; numbers.len()];
    for atom in 0..numbers.len() {
        let element = numbers[atom] as usize - 1;
        for (reference, &weight) in weights[atom]
            .iter()
            .enumerate()
            .take(int(NREF, element) as usize)
        {
            polarizabilities[atom] += weight * alpha(element, reference, 0, model);
        }
    }
    Ok(Properties {
        coordination,
        charges,
        c6,
        polarizabilities,
    })
}

pub(crate) fn validate(numbers: &[i32], positions: &[f64]) -> Result<(), &'static str> {
    if positions.len() != 3 * numbers.len() || numbers.is_empty() {
        return Err("invalid D4 structure shape");
    }
    if numbers
        .iter()
        .any(|&number| !(1..=118).contains(&number) || (104..=111).contains(&number))
    {
        return Err("D4 does not support this atomic number");
    }
    Ok(())
}

fn validate_cutoff(cutoff: Cutoff) -> Result<(), &'static str> {
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
        return Err("D4 cutoffs must be finite and nonnegative");
    }
    Ok(())
}

fn coordination_numbers(numbers: &[i32], positions: &[f64], cutoff: f64) -> Vec<f64> {
    let mut result = vec![0.0; numbers.len()];
    for first in 0..numbers.len() {
        let first_element = numbers[first] as usize - 1;
        for second in 0..first {
            let second_element = numbers[second] as usize - 1;
            let distance = (0..3)
                .map(|axis| (positions[3 * first + axis] - positions[3 * second + axis]).powi(2))
                .sum::<f64>()
                .sqrt();
            if distance > cutoff || distance < 1.0e-12 {
                continue;
            }
            let radius = element(first_element, 0) + element(second_element, 0);
            let count = 0.5 * (1.0 + unsafe { erf(-7.5 * (distance - radius) / radius) });
            let difference = (element(first_element, 3) - element(second_element, 3)).abs();
            let electronegativity =
                4.10451 * (-(difference + 19.08857).powi(2) / (2.0 * 11.28174_f64.powi(2))).exp();
            let contribution = count * electronegativity;
            result[first] += contribution;
            result[second] += contribution;
        }
    }
    result
}

fn charges(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
) -> Result<Vec<f64>, &'static str> {
    model.validate_selection(numbers.len())?;
    if let Some(charges) = model.fixed_charges {
        return Ok(charges.to_vec());
    }
    if model.eeqbc {
        let coordinates: Vec<_> = positions
            .iter()
            .map(|&value| Dual::constant(value, 0))
            .collect();
        return eeqbc_charges(numbers, &coordinates, charge, None)
            .map(|values| positions_values(&values));
    }
    molecular_charges(numbers, positions, charge, false, 1.0)
}

pub(crate) fn molecular_charges(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    clamp_cn: bool,
    radius_scale: f64,
) -> Result<Vec<f64>, &'static str> {
    let charge_coordination = charge_coordination_numbers(numbers, positions, radius_scale);
    let size = numbers.len() + 1;
    let mut matrix = vec![0.0; size * size];
    let mut rhs = vec![0.0; size];
    for first in 0..numbers.len() {
        let first_element = numbers[first] as usize - 1;
        let coordination = if clamp_cn {
            charge_coordination[first].max(f64::EPSILON)
        } else {
            charge_coordination[first] + 1.0e-14
        };
        rhs[first] = -eeq(first_element, 0) + eeq(first_element, 2) * coordination.sqrt();
        matrix[first * size + first] =
            eeq(first_element, 1) + (2.0 / PI).sqrt() / eeq(first_element, 3);
        matrix[first * size + numbers.len()] = 1.0;
        matrix[numbers.len() * size + first] = 1.0;
        for second in 0..first {
            let second_element = numbers[second] as usize - 1;
            let distance2 = (0..3)
                .map(|axis| (positions[3 * first + axis] - positions[3 * second + axis]).powi(2))
                .sum::<f64>();
            let inverse_width =
                1.0 / (eeq(first_element, 3).powi(2) + eeq(second_element, 3).powi(2)).sqrt();
            let interaction = unsafe { erf(inverse_width * distance2.sqrt()) } / distance2.sqrt();
            matrix[first * size + second] = interaction;
            matrix[second * size + first] = interaction;
        }
    }
    rhs[numbers.len()] = charge;
    solve(&mut matrix, &mut rhs)?;
    rhs.pop();
    Ok(rhs)
}

pub struct ChargeResponse {
    pub charges: Vec<f64>,
    pub cartesian: Vec<f64>,
    pub strain: Vec<f64>,
}

pub struct PropertyResponse {
    pub properties: Properties,
    pub coordination_cartesian: Vec<f64>,
    pub coordination_strain: Vec<f64>,
    pub charge_cartesian: Vec<f64>,
    pub charge_strain: Vec<f64>,
    pub c6_cartesian: Vec<f64>,
    pub c6_strain: Vec<f64>,
    pub polarizability_cartesian: Vec<f64>,
    pub polarizability_strain: Vec<f64>,
}

fn response_strain(
    gradient: &[f64],
    positions: &[f64],
    cell: Option<(&[f64; 9], [bool; 3])>,
) -> [f64; 9] {
    std::array::from_fn(|component| {
        let row = component % 3;
        let column = component / 3;
        let mut value: f64 = positions
            .as_chunks::<3>()
            .0
            .iter()
            .enumerate()
            .map(|(atom, xyz)| gradient[3 * atom + row] * xyz[column])
            .sum();
        if let Some((lattice, _)) = cell {
            value += (0..3)
                .map(|vector| {
                    gradient[positions.len() + row + 3 * vector] * lattice[column + 3 * vector]
                })
                .sum::<f64>();
        }
        value
    })
}

/// Property-major derivatives: coordinate = 3*atom+axis, strain = row+3*column.
/// Includes the charge-model response; strain transforms coordinates and lattice.
/// The full Cartesian C6 Jacobian requires 3*N^3 elements.
pub fn property_response(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    cell: Option<(&[f64; 9], [bool; 3])>,
) -> Result<PropertyResponse, &'static str> {
    let response_size = numbers
        .len()
        .checked_pow(3)
        .and_then(|size| size.checked_mul(3))
        .filter(|&size| size <= isize::MAX as usize / std::mem::size_of::<f64>())
        .ok_or("D4 property response dimensions overflow")?;
    model.validate_selection(numbers.len())?;
    let charge_response = charge_response(numbers, positions, charge, model, cell, true)?;
    let cell = cell.filter(|(_, periodic)| periodic.iter().any(|&active| active));
    let wrapped;
    let positions = if let Some((lattice, periodic)) = cell {
        wrapped = wrap_positions(positions, lattice, periodic)?;
        wrapped.as_slice()
    } else {
        positions
    };
    let atoms = numbers.len();
    let coordinates = positions.len();
    let mut cn_cartesian = Vec::with_capacity(atoms * coordinates);
    let mut cn_strain = Vec::with_capacity(atoms * 9);
    let coordination = if let Some((lattice, periodic)) = cell {
        TAPE.with(|tape| unsafe {
            let tape = &mut *tape.get();
            tape.nodes.clear();
            tape.solves.clear();
        });
        let size = coordinates + 9;
        let xyz: Vec<_> = positions
            .iter()
            .enumerate()
            .map(|(index, &value)| Dual::variable(value, size, index))
            .collect();
        let dual_lattice =
            std::array::from_fn(|index| Dual::variable(lattice[index], size, coordinates + index));
        let values = periodic_dual_coordination(
            numbers,
            &xyz,
            &dual_lattice,
            lattice,
            periodic,
            false,
            model.cutoff.cn,
        );
        for value in &values {
            let gradient = value.gradient();
            cn_cartesian.extend_from_slice(&gradient[..coordinates]);
            cn_strain.extend_from_slice(&response_strain(&gradient, positions, cell));
        }
        positions_values(&values)
    } else {
        let (values, jacobian) =
            molecular_coordination_jacobian(numbers, positions, false, model.cutoff.cn, 1.0);
        for gradient in jacobian.chunks_exact(coordinates) {
            cn_strain.extend_from_slice(&response_strain(gradient, positions, None));
        }
        cn_cartesian = jacobian;
        values
    };
    let partials =
        model_coefficient_partials(numbers, &coordination, &charge_response.charges, model);
    let mut result = PropertyResponse {
        properties: Properties {
            coordination,
            charges: charge_response.charges,
            c6: partials.iter().map(|values| values[0]).collect(),
            polarizabilities: vec![0.0; atoms],
        },
        coordination_cartesian: cn_cartesian,
        coordination_strain: cn_strain,
        charge_cartesian: charge_response.cartesian,
        charge_strain: charge_response.strain,
        c6_cartesian: vec![0.0; response_size],
        c6_strain: vec![0.0; atoms * atoms * 9],
        polarizability_cartesian: vec![0.0; atoms * coordinates],
        polarizability_strain: vec![0.0; atoms * 9],
    };
    for (atom, &number) in numbers.iter().enumerate() {
        let element = number as usize - 1;
        let weights = weight_partials(
            element,
            element,
            result.properties.coordination[atom],
            result.properties.charges[atom],
            model,
        );
        let mut polar = [0.0; 3];
        for (reference, weight) in weights.iter().enumerate().take(int(NREF, element) as usize) {
            let alpha = alpha(element, reference, 0, model);
            for component in 0..3 {
                polar[component] += weight[component] * alpha;
            }
        }
        result.properties.polarizabilities[atom] = polar[0];
        for (output, cn, charge, size) in [
            (
                &mut result.polarizability_cartesian,
                &result.coordination_cartesian,
                &result.charge_cartesian,
                coordinates,
            ),
            (
                &mut result.polarizability_strain,
                &result.coordination_strain,
                &result.charge_strain,
                9,
            ),
        ] {
            for component in 0..size {
                let index = atom * size + component;
                output[index] = polar[1] * cn[index] + polar[2] * charge[index];
            }
        }
    }
    for (pair, partial) in partials.iter().enumerate() {
        let first = pair / atoms;
        let second = pair % atoms;
        for (output, cn, charge, size) in [
            (
                &mut result.c6_cartesian,
                &result.coordination_cartesian,
                &result.charge_cartesian,
                coordinates,
            ),
            (
                &mut result.c6_strain,
                &result.coordination_strain,
                &result.charge_strain,
                9,
            ),
        ] {
            for component in 0..size {
                output[pair * size + component] = partial[1] * cn[first * size + component]
                    + partial[2] * cn[second * size + component]
                    + partial[3] * charge[first * size + component]
                    + partial[4] * charge[second * size + component];
            }
        }
    }
    Ok(result)
}

pub fn charge_response(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    cell: Option<(&[f64; 9], [bool; 3])>,
    derivatives: bool,
) -> Result<ChargeResponse, &'static str> {
    validate(numbers, positions)?;
    model.validate_selection(numbers.len())?;
    if !charge.is_finite() || positions.iter().any(|value| !value.is_finite()) {
        return Err("nonfinite charge or coordinates");
    }
    let cell = cell.filter(|(_, periodic)| periodic.iter().any(|&active| active));
    let wrapped;
    let positions = if let Some((lattice, periodic)) = cell {
        wrapped = wrap_positions(positions, lattice, periodic)?;
        wrapped.as_slice()
    } else {
        positions
    };
    let coordinate_count = positions.len();
    if let Some(charges) = model.fixed_charges {
        return Ok(ChargeResponse {
            charges: charges.to_vec(),
            cartesian: vec![
                0.0;
                if derivatives {
                    numbers.len() * coordinate_count
                } else {
                    0
                }
            ],
            strain: vec![0.0; if derivatives { numbers.len() * 9 } else { 0 }],
        });
    }
    let size = if derivatives { coordinate_count + 9 } else { 0 };
    TAPE.with(|tape| unsafe {
        let tape = &mut *tape.get();
        tape.nodes.clear();
        tape.solves.clear();
    });
    let coordinates: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(index, &value)| Dual::variable(value, size, index))
        .collect();
    let dual_charges = if let Some((lattice, periodic)) = cell {
        let dual_lattice = std::array::from_fn(|index| {
            Dual::variable(lattice[index], size, coordinate_count + index)
        });
        Some(periodic_dual_charges::<false>(
            numbers,
            &coordinates,
            charge,
            &dual_lattice,
            lattice,
            periodic,
            model,
        )?)
    } else if model.eeqbc {
        Some(eeqbc_charges(numbers, &coordinates, charge, None)?)
    } else {
        None
    };
    let values = if let Some(values) = &dual_charges {
        positions_values(values)
    } else {
        charges(numbers, positions, charge, model)?
    };
    let mut result = ChargeResponse {
        charges: values,
        cartesian: Vec::new(),
        strain: Vec::new(),
    };
    if derivatives {
        for atom in 0..numbers.len() {
            let gradient = if let Some(values) = &dual_charges {
                values[atom].gradient()
            } else {
                let mut selection = vec![0.0; numbers.len()];
                selection[atom] = 1.0;
                let mut gradient =
                    molecular_charge_adjoint(numbers, positions, &result.charges, &selection)?;
                gradient.resize(size, 0.0);
                gradient
            };
            result
                .cartesian
                .extend_from_slice(&gradient[..coordinate_count]);
            result
                .strain
                .extend_from_slice(&response_strain(&gradient, positions, cell));
        }
    }
    Ok(result)
}

#[test]
fn property_responses_match_differences() {
    let numbers = [6, 8, 1];
    let positions = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, -1.0, 2.0, 0.5];
    let lattice = [9.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.5, 0.2, 11.0];
    let step = 1.0e-5;
    let mut eeqbc = Model::D4;
    eeqbc.set_charge_model(1).unwrap();
    let fixed = [0.2, -0.4, 0.1];
    for model in [
        Model::D4,
        Model::D4S,
        eeqbc,
        Model::D4.with_fixed_charges(Some(&fixed)).unwrap(),
        Model::D4S.with_fixed_charges(Some(&fixed)).unwrap(),
        eeqbc.with_fixed_charges(Some(&fixed)).unwrap(),
    ] {
        for periodic in [[false; 3], [true, true, false], [true; 3]] {
            let response =
                property_response(&numbers, &positions, 0.2, model, Some((&lattice, periodic)))
                    .unwrap();
            let reference =
                periodic_properties(&numbers, &positions, 0.2, model, &lattice, periodic).unwrap();
            for (actual, expected) in [
                (&response.properties.c6, &reference.c6),
                (
                    &response.properties.polarizabilities,
                    &reference.polarizabilities,
                ),
            ] {
                for (actual, expected) in actual.iter().zip(expected) {
                    assert!((actual - expected).abs() < 1.0e-10);
                }
            }
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
                    periodic_properties(&numbers, &xyz, 0.2, model, &cell, periodic).unwrap()
                };
                let plus = evaluate(1.0);
                let minus = evaluate(-1.0);
                for (plus, minus, cartesian, strain) in [
                    (
                        &plus.coordination,
                        &minus.coordination,
                        &response.coordination_cartesian,
                        &response.coordination_strain,
                    ),
                    (
                        &plus.charges,
                        &minus.charges,
                        &response.charge_cartesian,
                        &response.charge_strain,
                    ),
                    (
                        &plus.c6,
                        &minus.c6,
                        &response.c6_cartesian,
                        &response.c6_strain,
                    ),
                    (
                        &plus.polarizabilities,
                        &minus.polarizabilities,
                        &response.polarizability_cartesian,
                        &response.polarizability_strain,
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
                                < 3.0e-5,
                            "property {property}, component {component}, analytical {analytical}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn standalone_charge_response_matches_properties_and_differences() {
    let numbers = [6, 8, 1];
    let positions = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, -1.0, 2.0, 0.5];
    let lattice = [12.0, 0.0, 0.0, 1.0, 13.0, 0.0, 0.5, 0.2, 14.0];
    for charge_model in 0..=1 {
        let mut model = Model::D4;
        model.set_charge_model(charge_model).unwrap();
        for dimensions in 0..=3 {
            let periodic = std::array::from_fn(|axis| axis < dimensions);
            let cell = Some((&lattice, periodic));
            let response = charge_response(&numbers, &positions, 0.25, model, cell, true).unwrap();
            let reference =
                periodic_properties(&numbers, &positions, 0.25, model, &lattice, periodic).unwrap();
            for atom in 0..3 {
                assert!((response.charges[atom] - reference.charges[atom]).abs() < 1e-12);
            }
            for coordinate in 0..9 {
                let mut plus = positions;
                let mut minus = positions;
                plus[coordinate] += 1e-5;
                minus[coordinate] -= 1e-5;
                let plus = charge_response(&numbers, &plus, 0.25, model, cell, false).unwrap();
                let minus = charge_response(&numbers, &minus, 0.25, model, cell, false).unwrap();
                for atom in 0..3 {
                    let difference = (plus.charges[atom] - minus.charges[atom]) / 2e-5;
                    assert!((response.cartesian[9 * atom + coordinate] - difference).abs() < 1e-7);
                }
            }
        }
    }
}

mod ewald;
mod fourier;
mod hessian;
#[cfg(test)]
mod tests;
pub use hessian::hessian;

fn eeqbc(element: usize, field: usize) -> f64 {
    real(EEQBC_DATA, field + 10 * element)
}

fn eeqbc_charges(
    numbers: &[i32],
    coordinates: &[Dual],
    charge: f64,
    cell: Option<(&[Dual; 9], [bool; 3])>,
) -> Result<Vec<Dual>, &'static str> {
    let periodic = cell.map_or([false; 3], |(_, periodic)| periodic);
    let lattice = cell.map(|(lattice, _)| lattice);
    validate(numbers, &positions_values(coordinates))?;
    if numbers
        .iter()
        .any(|&number| number as usize > EEQBC_ELEMENTS)
    {
        return Err("EEQBC supports atomic numbers 1 through 103");
    }
    if !charge.is_finite() {
        return Err("EEQBC total charge must be finite");
    }
    if coordinates.iter().any(|value| !value.value.is_finite()) {
        return Err("EEQBC coordinates must be finite");
    }
    let atoms = numbers.len();
    let size = coordinates[0].size;
    let constant = |value| Dual::constant(value, size);
    let positions = positions_values(coordinates);
    let lattice_values: [f64; 9] =
        std::array::from_fn(|index| lattice.map_or(0.0, |values| values[index].value));
    let distance = |first: usize, second: usize, image: [i32; 3]| {
        let mut distance2 = constant(0.0);
        for axis in 0..3 {
            let mut difference =
                coordinates[3 * first + axis].clone() - coordinates[3 * second + axis].clone();
            if let Some(lattice) = lattice {
                for column in 0..3 {
                    difference =
                        difference - lattice[axis + 3 * column].clone() * image[column] as f64;
                }
            }
            distance2 += difference.clone() * difference;
        }
        distance2.sqrt()
    };
    let images = if lattice.is_some() {
        lattice_indices(&lattice_values, periodic, 25.0)
    } else {
        vec![[0; 3]]
    };
    let mut coordination = vec![constant(0.0); atoms];
    let mut local_charge = vec![constant(charge / atoms as f64); atoms];
    for first in 0..atoms {
        let first_element = numbers[first] as usize - 1;
        for second in 0..=first {
            let second_element = numbers[second] as usize - 1;
            let radius = eeqbc(first_element, 7) + eeqbc(second_element, 7);
            for &image in &images {
                let vector = image_displacement(&positions, &lattice_values, first, second, image);
                let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
                if first != second && distance2 < 1.0e-12 {
                    return Err("coincident atoms in EEQBC structure");
                }
                if distance2 > 625.0 || distance2 < 1.0e-12 {
                    continue;
                }
                let count = (((distance(first, second, image) - constant(radius))
                    * (-2.0 / radius.powf(0.75)))
                .erf()
                    + constant(1.0))
                    * 0.5;
                coordination[first] += count.clone();
                if first != second {
                    coordination[second] += count.clone();
                    let transfer = count * (eeqbc(second_element, 9) - eeqbc(first_element, 9));
                    local_charge[first] += transfer.clone();
                    local_charge[second] = local_charge[second].clone() - transfer;
                }
            }
        }
    }
    let mut widths = Vec::with_capacity(atoms);
    let mut hardness = Vec::with_capacity(atoms);
    let mut electronegativity = Vec::with_capacity(atoms);
    for atom in 0..atoms {
        let element = numbers[atom] as usize - 1;
        let width = (constant(1.0)
            - coordination[atom].clone() * (0.14 / eeqbc(element, 8).powf(0.75)))
            * eeqbc(element, 2);
        if width.value <= 0.0 || !width.value.is_finite() {
            return Err("nonpositive effective EEQBC charge width");
        }
        hardness.push(
            constant(eeqbc(element, 1))
                + local_charge[atom].clone() * eeqbc(element, 5)
                + constant((2.0 / PI).sqrt()) / width.clone(),
        );
        widths.push(width);
        electronegativity.push(
            constant(-eeqbc(element, 0))
                + coordination[atom].clone() * eeqbc(element, 3)
                + local_charge[atom].clone() * eeqbc(element, 4),
        );
    }
    let dimension = atoms + 1;
    let mut matrix = vec![constant(0.0); dimension * dimension];
    let mut rhs = vec![constant(0.0); dimension];
    let direct_images = if lattice.is_some() {
        fixed_indices(2, true)
            .into_iter()
            .filter(|image| (0..3).all(|axis| periodic[axis] || image[axis] == 0))
            .collect()
    } else {
        vec![[0; 3]]
    };
    for first in 0..atoms {
        let first_element = numbers[first] as usize - 1;
        matrix[first * dimension + first] += constant(1.0);
        matrix[first * dimension + atoms] = constant(1.0);
        matrix[atoms * dimension + first] = constant(1.0);
        for second in 0..=first {
            if first == second && lattice.is_none() {
                continue;
            }
            let second_element = numbers[second] as usize - 1;
            let radius = real(EEQBC_VDW, first_element + EEQBC_ELEMENTS * second_element);
            let inverse_width = constant(1.0)
                / (widths[first].clone() * widths[first].clone()
                    + widths[second].clone() * widths[second].clone())
                .sqrt();
            let mut nearest = if lattice.is_some() {
                closest_directional_images(&positions, &lattice_values, first, second, periodic)
            } else {
                vec![[0; 3]]
            };
            if first == second && periodic.iter().all(|&active| active) {
                let mut ordered = fixed_indices(1, true);
                ordered.sort_by_key(|image| {
                    (
                        image[0].abs(),
                        image[1].abs(),
                        image[2].abs(),
                        image[0] < 0,
                        image[1] < 0,
                        image[2] < 0,
                    )
                });
                for image in &mut nearest {
                    let index = ordered
                        .iter()
                        .position(|candidate| candidate == image)
                        .unwrap();
                    *image = ordered[index - 1];
                }
            }
            let mut capacitance = constant(0.0);
            let mut interaction = constant(0.0);
            for base in &nearest {
                for direct in &direct_images {
                    let image = std::array::from_fn(|axis| base[axis] - direct[axis]);
                    let vector =
                        image_displacement(&positions, &lattice_values, first, second, image);
                    if vector.iter().map(|value| value * value).sum::<f64>() < f64::EPSILON {
                        continue;
                    }
                    let separation = distance(first, second, image);
                    let pair_capacitance =
                        (((separation.clone() - constant(radius)) * (-0.60 / radius)).erf()
                            + constant(1.0))
                            * (0.5 * (eeqbc(first_element, 6) * eeqbc(second_element, 6)).sqrt()
                                / nearest.len() as f64);
                    interaction += -pair_capacitance.clone()
                        * (inverse_width.clone() * separation.clone()).erf()
                        / separation;
                    capacitance += pair_capacitance;
                }
            }
            if first == second {
                matrix[first * dimension + first] +=
                    interaction + capacitance * hardness[first].clone();
                continue;
            }
            matrix[first * dimension + second] = interaction.clone();
            matrix[second * dimension + first] = interaction;
            matrix[first * dimension + first] += capacitance.clone() * hardness[first].clone();
            matrix[second * dimension + second] += capacitance.clone() * hardness[second].clone();
            let transfer = capacitance
                * (electronegativity[first].clone() - electronegativity[second].clone());
            rhs[first] += transfer.clone();
            rhs[second] = rhs[second].clone() - transfer;
        }
    }
    rhs[atoms] = constant(charge);
    dual_solve(&mut matrix, &mut rhs)?;
    rhs.pop();
    Ok(rhs)
}

fn charge_coordination_numbers(numbers: &[i32], positions: &[f64], radius_scale: f64) -> Vec<f64> {
    let mut result = vec![0.0; numbers.len()];
    for first in 0..numbers.len() {
        for second in 0..first {
            let distance = (0..3)
                .map(|axis| (positions[3 * first + axis] - positions[3 * second + axis]).powi(2))
                .sum::<f64>()
                .sqrt();
            if distance > 25.0 {
                continue;
            }
            let radius = (real(CHARGE_RCOV, numbers[first] as usize - 1)
                + real(CHARGE_RCOV, numbers[second] as usize - 1))
                * radius_scale;
            let contribution = 0.5 * (1.0 + unsafe { erf(-7.5 * (distance - radius) / radius) });
            result[first] += contribution;
            result[second] += contribution;
        }
    }
    result
        .into_iter()
        .map(|value| (1.0 + 8.0_f64.exp()).ln() - (1.0 + (8.0 - value).exp()).ln())
        .collect()
}

fn molecular_coordination_jacobian(
    numbers: &[i32],
    positions: &[f64],
    charge_model: bool,
    cutoff: f64,
    radius_scale: f64,
) -> (Vec<f64>, Vec<f64>) {
    let atoms = numbers.len();
    let coordinates = positions.len();
    let mut values = vec![0.0; atoms];
    let mut jacobian = vec![0.0; atoms * coordinates];
    for first in 0..atoms {
        for second in 0..first {
            let vector: [f64; 3] = std::array::from_fn(|axis| {
                positions[3 * first + axis] - positions[3 * second + axis]
            });
            let distance = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
            if distance > cutoff || distance < 1.0e-12 {
                continue;
            }
            let first_element = numbers[first] as usize - 1;
            let second_element = numbers[second] as usize - 1;
            let radius = if charge_model {
                (real(CHARGE_RCOV, first_element) + real(CHARGE_RCOV, second_element))
                    * radius_scale
            } else {
                element(first_element, 0) + element(second_element, 0)
            };
            let argument = -7.5 * (distance - radius) / radius;
            let factor = if charge_model {
                1.0
            } else {
                let difference = (element(first_element, 3) - element(second_element, 3)).abs();
                4.10451 * (-(difference + 19.08857).powi(2) / (2.0 * 11.28174_f64.powi(2))).exp()
            };
            let count = 0.5 * (1.0 + unsafe { erf(argument) }) * factor;
            let radial = -7.5 / radius / PI.sqrt() * (-argument * argument).exp() * factor;
            values[first] += count;
            values[second] += count;
            for axis in 0..3 {
                let derivative = radial * vector[axis] / distance;
                for atom in [first, second] {
                    jacobian[atom * coordinates + 3 * first + axis] += derivative;
                    jacobian[atom * coordinates + 3 * second + axis] -= derivative;
                }
            }
        }
    }
    if charge_model {
        for atom in 0..atoms {
            let scale = 8.0_f64.exp() / (8.0_f64.exp() + values[atom].exp());
            for coordinate in 0..coordinates {
                jacobian[atom * coordinates + coordinate] *= scale;
            }
            values[atom] = (1.0 + 8.0_f64.exp()).ln() - (1.0 + (8.0 - values[atom]).exp()).ln();
        }
    }
    (values, jacobian)
}

fn molecular_charge_adjoint(
    numbers: &[i32],
    positions: &[f64],
    charges: &[f64],
    dedq: &[f64],
) -> Result<Vec<f64>, &'static str> {
    molecular_charge_adjoint_with_cn_clamp(numbers, positions, charges, dedq, false, 1.0)
}

pub(crate) fn molecular_charge_adjoint_with_cn_clamp(
    numbers: &[i32],
    positions: &[f64],
    charges: &[f64],
    dedq: &[f64],
    clamp_cn: bool,
    radius_scale: f64,
) -> Result<Vec<f64>, &'static str> {
    let atoms = numbers.len();
    let coordinates = positions.len();
    let dimension = atoms + 1;
    let (coordination, coordination_jacobian) =
        molecular_coordination_jacobian(numbers, positions, true, 25.0, radius_scale);
    let mut matrix = vec![0.0; dimension * dimension];
    for first in 0..atoms {
        let first_element = numbers[first] as usize - 1;
        matrix[first * dimension + first] =
            eeq(first_element, 1) + (2.0 / PI).sqrt() / eeq(first_element, 3);
        matrix[first * dimension + atoms] = 1.0;
        matrix[atoms * dimension + first] = 1.0;
        for second in 0..first {
            let second_element = numbers[second] as usize - 1;
            let distance = squared_distance(positions, first, second).sqrt();
            let gamma =
                1.0 / (eeq(first_element, 3).powi(2) + eeq(second_element, 3).powi(2)).sqrt();
            let interaction = unsafe { erf(gamma * distance) } / distance;
            matrix[first * dimension + second] = interaction;
            matrix[second * dimension + first] = interaction;
        }
    }
    let mut lambda = vec![0.0; dimension];
    lambda[..atoms].copy_from_slice(dedq);
    solve(&mut matrix, &mut lambda)?;
    let mut gradient = vec![0.0; coordinates];
    for atom in 0..atoms {
        if clamp_cn && coordination[atom] < f64::EPSILON {
            continue;
        }
        let count = if clamp_cn {
            coordination[atom].max(f64::EPSILON)
        } else {
            coordination[atom] + 1.0e-14
        };
        let scale = lambda[atom] * 0.5 * eeq(numbers[atom] as usize - 1, 2) / count.sqrt();
        for coordinate in 0..coordinates {
            gradient[coordinate] += scale * coordination_jacobian[atom * coordinates + coordinate];
        }
    }
    for first in 0..atoms {
        for second in 0..first {
            let vector: [f64; 3] = std::array::from_fn(|component| {
                positions[3 * first + component] - positions[3 * second + component]
            });
            let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
            let distance = distance2.sqrt();
            let gamma = 1.0
                / (eeq(numbers[first] as usize - 1, 3).powi(2)
                    + eeq(numbers[second] as usize - 1, 3).powi(2))
                .sqrt();
            let radial = 2.0 * gamma * (-(gamma * distance).powi(2)).exp() / (PI.sqrt() * distance)
                - unsafe { erf(gamma * distance) } / distance2;
            let scale = -(lambda[first] * charges[second] + lambda[second] * charges[first])
                * radial
                / distance;
            for axis in 0..3 {
                let contribution = scale * vector[axis];
                gradient[3 * first + axis] += contribution;
                gradient[3 * second + axis] -= contribution;
            }
        }
    }
    Ok(gradient)
}

fn solve(matrix: &mut [f64], rhs: &mut [f64]) -> Result<(), &'static str> {
    let size = rhs.len();
    for column in 0..size {
        let pivot = (column..size)
            .max_by(|&left, &right| {
                matrix[left * size + column]
                    .abs()
                    .total_cmp(&matrix[right * size + column].abs())
            })
            .unwrap();
        if matrix[pivot * size + column].abs() < 1.0e-14 {
            return Err("singular EEQ system");
        }
        if pivot != column {
            for index in 0..size {
                matrix.swap(column * size + index, pivot * size + index);
            }
            rhs.swap(column, pivot);
        }
        for row in column + 1..size {
            let factor = matrix[row * size + column] / matrix[column * size + column];
            matrix[row * size + column] = 0.0;
            for index in column + 1..size {
                matrix[row * size + index] -= factor * matrix[column * size + index];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    for row in (0..size).rev() {
        rhs[row] = (rhs[row]
            - (row + 1..size)
                .map(|column| matrix[row * size + column] * rhs[column])
                .sum::<f64>())
            / matrix[row * size + row];
    }
    Ok(())
}

fn weights(
    atom_element: usize,
    other_element: usize,
    coordination: f64,
    charge: f64,
    model: Model,
) -> [f64; REFERENCES] {
    let count = int(NREF, atom_element) as usize;
    let mut base = [0.0; REFERENCES];
    for (reference_index, value) in base.iter_mut().enumerate().take(count) {
        let gaussian_count = int(NGW, reference_index + REFERENCES * atom_element) as usize;
        let width = if !model.d4s {
            model.wf
        } else {
            real(PAIR_WEIGHTS, other_element + ELEMENTS * atom_element)
        };
        *value = (1..=gaussian_count)
            .map(|gaussian| {
                (-width
                    * gaussian as f64
                    * (coordination - reference(REFCN, atom_element, reference_index)).powi(2))
                .exp()
            })
            .sum();
    }
    let normalization = base[..count].iter().sum::<f64>();
    let effective = charge + element(atom_element, 2);
    for (reference_index, value) in base.iter_mut().enumerate().take(count) {
        let reference_charge = reference(
            if model.eeqbc { EEQBC_REFQ } else { REFQ },
            atom_element,
            reference_index,
        ) + element(atom_element, 2);
        let zeta = if effective < 0.0 {
            model.ga.exp()
        } else {
            (model.ga
                * (1.0
                    - (model.gc * element(atom_element, 4) * (1.0 - reference_charge / effective))
                        .exp()))
            .exp()
        };
        *value = *value / normalization * zeta;
    }
    base
}

fn weight_partials(
    atom_element: usize,
    other_element: usize,
    coordination: f64,
    charge: f64,
    model: Model,
) -> [[f64; 3]; REFERENCES] {
    let count = int(NREF, atom_element) as usize;
    let width = if !model.d4s {
        model.wf
    } else {
        real(PAIR_WEIGHTS, other_element + ELEMENTS * atom_element)
    };
    let mut base = [0.0; REFERENCES];
    let mut base_derivative = [0.0; REFERENCES];
    for reference_index in 0..count {
        let difference = coordination - reference(REFCN, atom_element, reference_index);
        for gaussian in 1..=int(NGW, reference_index + REFERENCES * atom_element) {
            let exponent = (-width * gaussian as f64 * difference * difference).exp();
            base[reference_index] += exponent;
            base_derivative[reference_index] +=
                -2.0 * width * gaussian as f64 * difference * exponent;
        }
    }
    let normalization = base[..count].iter().sum::<f64>();
    let normalization_derivative = base_derivative[..count].iter().sum::<f64>();
    let effective = charge + element(atom_element, 2);
    let mut result = [[0.0; 3]; REFERENCES];
    for reference_index in 0..count {
        let normalized = base[reference_index] / normalization;
        let normalized_derivative = (base_derivative[reference_index] * normalization
            - base[reference_index] * normalization_derivative)
            / normalization.powi(2);
        let reference_charge = reference(
            if model.eeqbc { EEQBC_REFQ } else { REFQ },
            atom_element,
            reference_index,
        ) + element(atom_element, 2);
        let (zeta, zeta_derivative) = if effective < 0.0 {
            (model.ga.exp(), 0.0)
        } else {
            let steepness = model.gc * element(atom_element, 4);
            let inner = (steepness * (1.0 - reference_charge / effective)).exp();
            let zeta = (model.ga * (1.0 - inner)).exp();
            (
                zeta,
                -model.ga * zeta * inner * steepness * reference_charge / effective.powi(2),
            )
        };
        result[reference_index] = [
            normalized * zeta,
            normalized_derivative * zeta,
            normalized * zeta_derivative,
        ];
    }
    result
}

fn coefficient_partials(
    first_element: usize,
    second_element: usize,
    first: &[[f64; 3]; REFERENCES],
    second: &[[f64; 3]; REFERENCES],
    model: Model,
) -> [f64; 5] {
    let mut result = [0.0; 5];
    for (first_reference, first_weight) in first
        .iter()
        .enumerate()
        .take(int(NREF, first_element) as usize)
    {
        for (second_reference, second_weight) in second
            .iter()
            .enumerate()
            .take(int(NREF, second_element) as usize)
        {
            let reference = reference_c6(
                first_element,
                first_reference,
                second_element,
                second_reference,
                model,
            );
            result[0] += first_weight[0] * second_weight[0] * reference;
            result[1] += first_weight[1] * second_weight[0] * reference;
            result[2] += first_weight[0] * second_weight[1] * reference;
            result[3] += first_weight[2] * second_weight[0] * reference;
            result[4] += first_weight[0] * second_weight[2] * reference;
        }
    }
    result
}

fn model_coefficient_partials(
    numbers: &[i32],
    coordination: &[f64],
    charges: &[f64],
    model: Model,
) -> Vec<[f64; 5]> {
    let mut species = Vec::new();
    let indices: Vec<_> = numbers
        .iter()
        .map(|&number| {
            if !model.d4s {
                return 0;
            }
            species
                .iter()
                .position(|&known| known == number)
                .unwrap_or_else(|| {
                    species.push(number);
                    species.len() - 1
                })
        })
        .collect();
    if !model.d4s {
        species.push(numbers[0]);
    }
    let weights: Vec<Vec<_>> = numbers
        .iter()
        .enumerate()
        .map(|(atom, &number)| {
            species
                .iter()
                .map(|&other| {
                    weight_partials(
                        number as usize - 1,
                        other as usize - 1,
                        coordination[atom],
                        charges[atom],
                        model,
                    )
                })
                .collect()
        })
        .collect();
    let atoms = numbers.len();
    let mut result = vec![[0.0; 5]; atoms * atoms];
    if model.eeqbc || model.ga != 3.0 || model.gc != 2.0 {
        let mut spectra = vec![[[0.0; 3]; FREQUENCIES]; atoms * species.len()];
        let mut alphas = vec![None; ELEMENTS];
        for atom in 0..atoms {
            let element = numbers[atom] as usize - 1;
            let alphas = alphas[element].get_or_insert_with(|| {
                std::array::from_fn::<_, REFERENCES, _>(|reference| {
                    if reference < int(NREF, element) as usize {
                        custom_alpha(element, reference, model)
                    } else {
                        [0.0; FREQUENCIES]
                    }
                })
            });
            for reference in 0..int(NREF, element) as usize {
                for frequency in 0..FREQUENCIES {
                    let alpha = alphas[reference][frequency];
                    for other in 0..species.len() {
                        for component in 0..3 {
                            spectra[atom * species.len() + other][frequency][component] +=
                                weights[atom][other][reference][component] * alpha;
                        }
                    }
                }
            }
        }
        for first in 0..atoms {
            for second in 0..=first {
                let mut partials = [0.0; 5];
                for (frequency, (left, right)) in spectra[first * species.len() + indices[second]]
                    .iter()
                    .zip(&spectra[second * species.len() + indices[first]])
                    .enumerate()
                {
                    let scale = 3.0 / PI * quadrature(frequency);
                    for (total, value) in partials.iter_mut().zip([
                        left[0] * right[0],
                        left[1] * right[0],
                        left[0] * right[1],
                        left[2] * right[0],
                        left[0] * right[2],
                    ]) {
                        *total += scale * value;
                    }
                }
                result[first * atoms + second] = partials;
                result[second * atoms + first] = [
                    partials[0],
                    partials[2],
                    partials[1],
                    partials[4],
                    partials[3],
                ];
            }
        }
        return result;
    }
    for first in 0..atoms {
        for second in 0..=first {
            let partials = coefficient_partials(
                numbers[first] as usize - 1,
                numbers[second] as usize - 1,
                &weights[first][indices[second]],
                &weights[second][indices[first]],
                model,
            );
            result[first * atoms + second] = partials;
            result[second * atoms + first] = [
                partials[0],
                partials[2],
                partials[1],
                partials[4],
                partials[3],
            ];
        }
    }
    result
}

fn coefficients(numbers: &[i32], weights: &[[f64; REFERENCES]], model: Model) -> Vec<f64> {
    let mut c6 = vec![0.0; numbers.len() * numbers.len()];
    for first in 0..numbers.len() {
        let first_element = numbers[first] as usize - 1;
        for second in 0..=first {
            let second_element = numbers[second] as usize - 1;
            let mut coefficient = 0.0;
            for first_reference in 0..int(NREF, first_element) as usize {
                for second_reference in 0..int(NREF, second_element) as usize {
                    coefficient += weights[first][first_reference]
                        * weights[second][second_reference]
                        * reference_c6(
                            first_element,
                            first_reference,
                            second_element,
                            second_reference,
                            model,
                        );
                }
            }
            c6[first * numbers.len() + second] = coefficient;
            c6[second * numbers.len() + first] = coefficient;
        }
    }
    c6
}

fn model_coefficients(
    numbers: &[i32],
    coordination: &[f64],
    charges: &[f64],
    model: Model,
) -> Vec<f64> {
    if model.eeqbc || model.ga != 3.0 || model.gc != 2.0 {
        return model_coefficient_partials(numbers, coordination, charges, model)
            .into_iter()
            .map(|partials| partials[0])
            .collect();
    }
    if !model.d4s {
        let weights: Vec<_> = numbers
            .iter()
            .zip(coordination)
            .zip(charges)
            .map(|((&number, &cn), &charge)| {
                weights(number as usize - 1, number as usize - 1, cn, charge, model)
            })
            .collect();
        return coefficients(numbers, &weights, model);
    }
    let atoms = numbers.len();
    let mut species = Vec::new();
    let indices: Vec<_> = numbers
        .iter()
        .map(|&number| {
            species
                .iter()
                .position(|&known| known == number)
                .unwrap_or_else(|| {
                    species.push(number);
                    species.len() - 1
                })
        })
        .collect();
    let pair_weights: Vec<Vec<_>> = numbers
        .iter()
        .enumerate()
        .map(|(atom, &number)| {
            species
                .iter()
                .map(|&other| {
                    weights(
                        number as usize - 1,
                        other as usize - 1,
                        coordination[atom],
                        charges[atom],
                        model,
                    )
                })
                .collect()
        })
        .collect();
    let mut result = vec![0.0; atoms * atoms];
    for first in 0..atoms {
        for second in 0..=first {
            let first_element = numbers[first] as usize - 1;
            let second_element = numbers[second] as usize - 1;
            let mut coefficient = 0.0;
            for (first_reference, &first_weight) in pair_weights[first][indices[second]]
                .iter()
                .enumerate()
                .take(int(NREF, first_element) as usize)
            {
                for (second_reference, &second_weight) in pair_weights[second][indices[first]]
                    .iter()
                    .enumerate()
                    .take(int(NREF, second_element) as usize)
                {
                    coefficient += first_weight
                        * second_weight
                        * reference_c6(
                            first_element,
                            first_reference,
                            second_element,
                            second_reference,
                            model,
                        );
                }
            }
            result[first * atoms + second] = coefficient;
            result[second * atoms + first] = coefficient;
        }
    }
    result
}

fn two_body_energy(
    numbers: &[i32],
    positions: &[f64],
    c6: &[f64],
    param: Param,
    cutoff: Cutoff,
    model: Model,
) -> f64 {
    let mut result = 0.0;
    for first in 0..numbers.len() {
        for second in 0..first {
            if !model.owns_pair(first, second) {
                continue;
            }
            if squared_distance(positions, first, second) > cutoff.disp2 * cutoff.disp2 {
                continue;
            }
            result += pair_energy(
                numbers[first] as usize - 1,
                numbers[second] as usize - 1,
                squared_distance(positions, first, second),
                c6[first * numbers.len() + second],
                param,
            ) * smooth_cutoff(
                squared_distance(positions, first, second),
                cutoff.disp2,
                cutoff.width2,
            )
            .0;
        }
    }
    result
}

pub(crate) fn pair_energy(
    first: usize,
    second: usize,
    distance2: f64,
    c6: f64,
    param: Param,
) -> f64 {
    let rrij = 3.0 * element(first, 1) * element(second, 1);
    let radius = param.a1 * rrij.sqrt() + param.a2;
    -c6 * (param.s6 / (distance2.powi(3) + radius.powi(6))
        + param.s8 * rrij / (distance2.powi(4) + radius.powi(8)))
}

fn atm_energy(
    numbers: &[i32],
    positions: &[f64],
    c6: &[f64],
    param: Param,
    cutoff: Cutoff,
    model: Model,
) -> f64 {
    let pairs = atm_pair_data(numbers, positions, param, cutoff);
    crate::parallel::map(numbers.len(), 64, |start, stride| {
        let mut result = 0.0;
        for first in (start..numbers.len()).step_by(stride) {
            for second in 0..first {
                if !model.owns_pair(first, second) {
                    continue;
                }
                for third in 0..second {
                    if !model.active(third) {
                        continue;
                    }
                    result += atm_triplet(numbers.len(), &pairs, c6, param, [first, second, third]);
                }
            }
        }
        result
    })
    .into_iter()
    .sum()
}

fn atm_pair_data(
    numbers: &[i32],
    positions: &[f64],
    param: Param,
    cutoff: Cutoff,
) -> Vec<[f64; 3]> {
    let atoms = numbers.len();
    let mut pairs = vec![[0.0; 3]; atoms * atoms];
    for first in 0..atoms {
        for second in 0..first {
            let distance = squared_distance(positions, first, second);
            let rrij = 3.0
                * element(numbers[first] as usize - 1, 1)
                * element(numbers[second] as usize - 1, 1);
            let pair = [
                distance,
                param.a1 * rrij.sqrt() + param.a2,
                if distance > cutoff.disp3 * cutoff.disp3 {
                    0.0
                } else {
                    smooth_cutoff(distance, cutoff.disp3, cutoff.width3).0
                },
            ];
            pairs[first * atoms + second] = pair;
            pairs[second * atoms + first] = pair;
        }
    }
    pairs
}

fn atm_triplet(
    atoms: usize,
    pairs: &[[f64; 3]],
    c6: &[f64],
    param: Param,
    [first, second, third]: [usize; 3],
) -> f64 {
    let [distance_ij, radius_ij, switch_ij] = pairs[first * atoms + second];
    let [distance_ik, radius_ik, switch_ik] = pairs[first * atoms + third];
    let [distance_jk, radius_jk, switch_jk] = pairs[second * atoms + third];
    let switch = switch_ij * switch_ik * switch_jk;
    if switch == 0.0 {
        return 0.0;
    }
    let product2 = distance_ij * distance_ik * distance_jk;
    let product = product2.sqrt();
    let c9 = -param.s9
        * (c6[first * atoms + second] * c6[first * atoms + third] * c6[second * atoms + third])
            .abs()
            .sqrt();
    let radius = radius_ij * radius_ik * radius_jk;
    let damping = 1.0 / (1.0 + 6.0 * (radius / product).powf(param.alpha / 3.0));
    let angular = 0.375
        * (distance_ij + distance_jk - distance_ik)
        * (distance_ij - distance_jk + distance_ik)
        * (-distance_ij + distance_jk + distance_ik)
        / (product2 * product2 * product)
        + 1.0 / (product2 * product);
    -c9 * damping * angular * switch
}

fn smooth_atm(mut local: Local, distances: [f64; 3], cutoff: Cutoff) -> Local {
    let switches = distances.map(|distance| smooth_cutoff(distance, cutoff.disp3, cutoff.width3));
    let product = switches.iter().map(|switch| switch.0).product::<f64>();
    for index in 0..3 {
        local.gradient[index] = local.gradient[index] * product
            + local.value
                * switches[index].1
                * switches[(index + 1) % 3].0
                * switches[(index + 2) % 3].0;
        local.gradient[index + 3] *= product;
    }
    local.value *= product;
    local
}

#[derive(Clone)]
struct Dual {
    value: f64,
    node: Option<usize>,
    size: usize,
}

struct TapeNode {
    parents: [(usize, f64); 2],
    curvature: [[f64; 2]; 2],
    parent_count: usize,
    input: Option<usize>,
    solve: Option<(usize, usize)>,
}

struct SolveNode {
    matrix: Vec<Option<usize>>,
    rhs: Vec<Option<usize>>,
    values: Vec<f64>,
    solution: Vec<f64>,
}

#[derive(Default)]
struct Tape {
    nodes: Vec<TapeNode>,
    solves: Vec<SolveNode>,
}

thread_local! {
    static TAPE: UnsafeCell<Tape> = const { UnsafeCell::new(Tape { nodes: Vec::new(), solves: Vec::new() }) };
}

impl Dual {
    fn constant(value: f64, size: usize) -> Self {
        Self {
            value,
            node: None,
            size,
        }
    }

    fn variable(value: f64, size: usize, index: usize) -> Self {
        if size == 0 {
            return Self::constant(value, size);
        }
        let node = TAPE.with(|tape| {
            let tape = unsafe { &mut *tape.get() };
            let node = tape.nodes.len();
            tape.nodes.push(TapeNode {
                parents: [(0, 0.0); 2],
                curvature: [[0.0; 2]; 2],
                parent_count: 0,
                input: Some(index),
                solve: None,
            });
            node
        });
        Self {
            value,
            node: Some(node),
            size,
        }
    }

    fn unary(self, value: f64, derivative: f64) -> Self {
        self.unary_second(value, derivative, 0.0)
    }

    fn unary_second(self, value: f64, derivative: f64, curvature: f64) -> Self {
        Self::record_second(
            value,
            self.size,
            self.node.map(|node| (node, derivative)),
            None,
            [[curvature, 0.0], [0.0, 0.0]],
        )
    }

    fn record(
        value: f64,
        size: usize,
        first: Option<(usize, f64)>,
        second: Option<(usize, f64)>,
    ) -> Self {
        Self::record_second(value, size, first, second, [[0.0; 2]; 2])
    }

    fn record_second(
        value: f64,
        size: usize,
        first: Option<(usize, f64)>,
        second: Option<(usize, f64)>,
        curvature: [[f64; 2]; 2],
    ) -> Self {
        if size == 0 || (first.is_none() && second.is_none()) {
            return Self::constant(value, size);
        }
        let mut parents = [(0, 0.0); 2];
        let mut indices = [0; 2];
        let mut parent_count = 0;
        for (index, parent) in [first, second].into_iter().enumerate() {
            if let Some(parent) = parent {
                parents[parent_count] = parent;
                indices[parent_count] = index;
                parent_count += 1;
            }
        }
        let node = TAPE.with(|tape| {
            let tape = unsafe { &mut *tape.get() };
            let node = tape.nodes.len();
            tape.nodes.push(TapeNode {
                parents,
                curvature: std::array::from_fn(|row| {
                    std::array::from_fn(|column| curvature[indices[row]][indices[column]])
                }),
                parent_count,
                input: None,
                solve: None,
            });
            node
        });
        Self {
            value,
            node: Some(node),
            size,
        }
    }

    fn gradient(&self) -> Vec<f64> {
        let mut result = vec![0.0; self.size];
        let Some(output) = self.node else {
            return result;
        };
        TAPE.with(|tape| {
            let tape = unsafe { &*tape.get() };
            let mut adjoints = vec![0.0; tape.nodes.len()];
            adjoints[output] = 1.0;
            for node in (0..=output).rev() {
                if let Some((operation, 0)) = tape.nodes[node].solve {
                    let system = &tape.solves[operation];
                    let dimension = system.solution.len();
                    let mut lambda = adjoints[node..node + dimension].to_vec();
                    let mut matrix = system.values.clone();
                    solve(&mut matrix, &mut lambda).expect("previously solved EEQ system");
                    for (row, &weight) in lambda.iter().enumerate() {
                        if let Some(parent) = system.rhs[row] {
                            adjoints[parent] += weight;
                        }
                        for column in 0..dimension {
                            if let Some(parent) = system.matrix[row * dimension + column] {
                                adjoints[parent] -= weight * system.solution[column];
                            }
                        }
                    }
                }
                let adjoint = adjoints[node];
                if adjoint == 0.0 {
                    continue;
                }
                if let Some(input) = tape.nodes[node].input {
                    result[input] += adjoint;
                }
                for &(parent, derivative) in
                    &tape.nodes[node].parents[..tape.nodes[node].parent_count]
                {
                    adjoints[parent] += adjoint * derivative;
                }
            }
        });
        result
    }

    fn hessian(&self) -> Result<Vec<f64>, &'static str> {
        // ponytail: refactor dense solves per column; cache LU factors if large-system Hessians need it.
        let count = self.size;
        let mut result = vec![0.0; count.checked_mul(count).ok_or("D4 Hessian size overflow")?];
        let Some(output) = self.node else {
            return Ok(result);
        };
        TAPE.with(|tape| {
            let tape = unsafe { &*tape.get() };
            for direction in 0..count {
                let mut tangent = vec![0.0; tape.nodes.len()];
                for (index, node) in tape.nodes[..=output].iter().enumerate() {
                    if let Some((operation, component)) = node.solve {
                        if component == 0 {
                            let system = &tape.solves[operation];
                            let dimension = system.solution.len();
                            let mut rhs: Vec<_> = (0..dimension)
                                .map(|row| {
                                    system.rhs[row].map_or(0.0, |parent| tangent[parent])
                                        - (0..dimension)
                                            .map(|column| {
                                                system.matrix[row * dimension + column]
                                                    .map_or(0.0, |parent| tangent[parent])
                                                    * system.solution[column]
                                            })
                                            .sum::<f64>()
                                })
                                .collect();
                            solve(&mut system.values.clone(), &mut rhs)?;
                            tangent[index..index + dimension].copy_from_slice(&rhs);
                        }
                    } else {
                        tangent[index] = f64::from(node.input == Some(direction))
                            + node.parents[..node.parent_count]
                                .iter()
                                .map(|&(parent, derivative)| tangent[parent] * derivative)
                                .sum::<f64>();
                    }
                }
                let mut adjoint = vec![0.0; tape.nodes.len()];
                let mut response = vec![0.0; tape.nodes.len()];
                adjoint[output] = 1.0;
                for index in (0..=output).rev() {
                    let node = &tape.nodes[index];
                    if let Some((operation, 0)) = node.solve {
                        let system = &tape.solves[operation];
                        let dimension = system.solution.len();
                        let mut lambda = adjoint[index..index + dimension].to_vec();
                        let transpose: Vec<_> = (0..dimension * dimension)
                            .map(|entry| {
                                system.values[(entry % dimension) * dimension + entry / dimension]
                            })
                            .collect();
                        solve(&mut transpose.clone(), &mut lambda)?;
                        let mut delta: Vec<_> = (0..dimension)
                            .map(|row| {
                                response[index + row]
                                    - (0..dimension)
                                        .map(|column| {
                                            system.matrix[column * dimension + row]
                                                .map_or(0.0, |parent| tangent[parent])
                                                * lambda[column]
                                        })
                                        .sum::<f64>()
                            })
                            .collect();
                        solve(&mut transpose.clone(), &mut delta)?;
                        for row in 0..dimension {
                            if let Some(parent) = system.rhs[row] {
                                adjoint[parent] += lambda[row];
                                response[parent] += delta[row];
                            }
                            for column in 0..dimension {
                                if let Some(parent) = system.matrix[row * dimension + column] {
                                    adjoint[parent] -= lambda[row] * system.solution[column];
                                    response[parent] -= delta[row] * system.solution[column]
                                        + lambda[row] * tangent[index + column];
                                }
                            }
                        }
                    }
                    if let Some(input) = node.input {
                        result[direction * count + input] += response[index];
                    }
                    for (slot, &(parent, derivative)) in
                        node.parents[..node.parent_count].iter().enumerate()
                    {
                        let delta = node.parents[..node.parent_count]
                            .iter()
                            .enumerate()
                            .map(|(other, &(source, _))| {
                                node.curvature[slot][other] * tangent[source]
                            })
                            .sum::<f64>();
                        response[parent] += response[index] * derivative + adjoint[index] * delta;
                        adjoint[parent] += adjoint[index] * derivative;
                    }
                }
            }
            Ok(result)
        })
    }

    fn exp(self) -> Self {
        let value = self.value.exp();
        self.unary_second(value, value, value)
    }

    fn erf(self) -> Self {
        let input = self.value;
        let derivative = 2.0 / PI.sqrt() * (-input * input).exp();
        self.unary_second(unsafe { erf(input) }, derivative, -2.0 * input * derivative)
    }

    fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        self.unary_second(value, 0.5 / value, -0.25 / value.powi(3))
    }

    fn ln(self) -> Self {
        let input = self.value;
        self.unary_second(input.ln(), input.recip(), -input.recip().powi(2))
    }

    fn powf(self, exponent: f64) -> Self {
        let input = self.value;
        self.unary_second(
            input.powf(exponent),
            exponent * input.powf(exponent - 1.0),
            exponent * (exponent - 1.0) * input.powf(exponent - 2.0),
        )
    }
}

impl Add for Dual {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::record(
            self.value + rhs.value,
            self.size,
            self.node.map(|node| (node, 1.0)),
            rhs.node.map(|node| (node, 1.0)),
        )
    }
}

impl AddAssign for Dual {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl Add<f64> for Dual {
    type Output = Self;
    fn add(mut self, rhs: f64) -> Self {
        self.value += rhs;
        self
    }
}

impl Sub for Dual {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl Neg for Dual {
    type Output = Self;
    fn neg(self) -> Self {
        let value = -self.value;
        self.unary(value, -1.0)
    }
}

impl Mul<f64> for Dual {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let value = self.value * rhs;
        self.unary(value, rhs)
    }
}

impl Mul for Dual {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::record_second(
            self.value * rhs.value,
            self.size,
            self.node.map(|node| (node, rhs.value)),
            rhs.node.map(|node| (node, self.value)),
            [[0.0, 1.0], [1.0, 0.0]],
        )
    }
}

impl Div for Dual {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let numerator = self.value;
        let denominator = rhs.value;
        Self::record_second(
            numerator / denominator,
            self.size,
            self.node.map(|node| (node, denominator.recip())),
            rhs.node
                .map(|node| (node, -numerator / denominator.powi(2))),
            [
                [0.0, -denominator.recip().powi(2)],
                [
                    -denominator.recip().powi(2),
                    2.0 * numerator / denominator.powi(3),
                ],
            ],
        )
    }
}

impl Div<f64> for Dual {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        let value = self.value / rhs;
        self.unary(value, rhs.recip())
    }
}

fn differentiated_energy(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    cutoff: Cutoff,
) -> Result<(f64, Vec<f64>), &'static str> {
    validate_cutoff(cutoff)?;
    validate(numbers, positions)?;
    model.validate_selection(numbers.len())?;
    let atoms = numbers.len();
    let coordinate_count = positions.len();
    let (coordination_values, coordination_jacobian) =
        molecular_coordination_jacobian(numbers, positions, false, cutoff.cn, 1.0);
    let charge_values = charges(numbers, positions, charge, model)?;
    let charged = model_coefficient_partials(numbers, &coordination_values, &charge_values, model);
    let mut energy = 0.0;
    let mut gradient = vec![0.0; coordinate_count];
    let mut dedcn = vec![0.0; atoms];
    let mut dedq = vec![0.0; atoms];
    for first in 0..numbers.len() {
        for second in 0..first {
            if !model.owns_pair(first, second) {
                continue;
            }
            let distance2 = squared_distance(positions, first, second);
            if distance2 > cutoff.disp2 * cutoff.disp2 {
                continue;
            }
            let c6 = charged[first * atoms + second][0];
            let rrij = 3.0
                * element(numbers[first] as usize - 1, 1)
                * element(numbers[second] as usize - 1, 1);
            let radius = param.a1 * rrij.sqrt() + param.a2;
            let sixth = 1.0 / (distance2.powi(3) + radius.powi(6));
            let eighth = 1.0 / (distance2.powi(4) + radius.powi(8));
            let potential = param.s6 * sixth + param.s8 * rrij * eighth;
            let switch = smooth_cutoff(distance2, cutoff.disp2, cutoff.width2);
            energy -= c6 * potential * switch.0;
            add_pair_gradient(
                &mut gradient,
                positions,
                first,
                second,
                c6 * (switch.0
                    * (3.0 * param.s6 * distance2.powi(2) * sixth.powi(2)
                        + 4.0 * param.s8 * rrij * distance2.powi(3) * eighth.powi(2))
                    - potential * switch.1),
            );
            let potential = potential * switch.0;
            let partials = charged[first * atoms + second];
            dedcn[first] -= potential * partials[1];
            dedcn[second] -= potential * partials[2];
            dedq[first] -= potential * partials[3];
            dedq[second] -= potential * partials[4];
        }
    }
    if param.s9.abs() >= f64::EPSILON {
        let zero =
            model_coefficient_partials(numbers, &coordination_values, &vec![0.0; atoms], model);
        let contributions = crate::parallel::map(atoms, 64, |start, stride| {
            let mut energy = 0.0;
            let mut gradient = vec![0.0; coordinate_count];
            let mut coefficient_derivatives = vec![0.0; atoms * atoms];
            for first in (start..atoms).step_by(stride) {
                for second in 0..first {
                    if !model.owns_pair(first, second) {
                        continue;
                    }
                    for third in 0..second {
                        if !model.active(third) {
                            continue;
                        }
                        let distances = [
                            squared_distance(positions, first, second),
                            squared_distance(positions, first, third),
                            squared_distance(positions, second, third),
                        ];
                        if distances
                            .iter()
                            .any(|&distance| distance > cutoff.disp3 * cutoff.disp3)
                        {
                            continue;
                        }
                        let coefficients = [
                            zero[first * atoms + second][0],
                            zero[first * atoms + third][0],
                            zero[second * atoms + third][0],
                        ];
                        let radius = [(first, second), (first, third), (second, third)]
                            .iter()
                            .map(|&(left, right)| {
                                let rrij = 3.0
                                    * element(numbers[left] as usize - 1, 1)
                                    * element(numbers[right] as usize - 1, 1);
                                param.a1 * rrij.sqrt() + param.a2
                            })
                            .product();
                        let local = smooth_atm(
                            local_atm(distances, coefficients, radius, param),
                            distances,
                            cutoff,
                        );
                        energy += local.value;
                        for (&(left, right), &derivative) in
                            [(first, second), (first, third), (second, third)]
                                .iter()
                                .zip(local.gradient[..3].iter())
                        {
                            add_pair_gradient(&mut gradient, positions, left, right, derivative);
                        }
                        coefficient_derivatives[first * atoms + second] += local.gradient[3];
                        coefficient_derivatives[first * atoms + third] += local.gradient[4];
                        coefficient_derivatives[second * atoms + third] += local.gradient[5];
                    }
                }
            }
            (energy, gradient, coefficient_derivatives)
        });
        for (local_energy, local_gradient, coefficient_derivatives) in contributions {
            energy += local_energy;
            for (total, value) in gradient.iter_mut().zip(local_gradient) {
                *total += value;
            }
            for first in 0..atoms {
                for second in 0..first {
                    let partials = zero[first * atoms + second];
                    let derivative = coefficient_derivatives[first * atoms + second];
                    dedcn[first] += derivative * partials[1];
                    dedcn[second] += derivative * partials[2];
                }
            }
        }
    }
    for coordinate in 0..coordinate_count {
        for atom in 0..atoms {
            gradient[coordinate] +=
                dedcn[atom] * coordination_jacobian[atom * coordinate_count + coordinate];
        }
    }
    let charge_gradient = if model.fixed_charges.is_some() {
        vec![0.0; coordinate_count]
    } else if model.eeqbc {
        TAPE.with(|tape| unsafe {
            *tape.get() = Tape::default();
        });
        let coordinates: Vec<_> = positions
            .iter()
            .enumerate()
            .map(|(index, &value)| Dual::variable(value, coordinate_count, index))
            .collect();
        let charges = eeqbc_charges(numbers, &coordinates, charge, None)?;
        charges
            .into_iter()
            .zip(&dedq)
            .fold(
                Dual::constant(0.0, coordinate_count),
                |value, (charge, &partial)| value + charge * partial,
            )
            .gradient()
    } else {
        molecular_charge_adjoint(numbers, positions, &charge_values, &dedq)?
    };
    for (value, contribution) in gradient.iter_mut().zip(charge_gradient) {
        *value += contribution;
    }
    Ok((energy, gradient))
}

fn add_pair_gradient(
    gradient: &mut [f64],
    positions: &[f64],
    first: usize,
    second: usize,
    distance2_derivative: f64,
) {
    for axis in 0..3 {
        let contribution = 2.0
            * distance2_derivative
            * (positions[3 * first + axis] - positions[3 * second + axis]);
        gradient[3 * first + axis] += contribution;
        gradient[3 * second + axis] -= contribution;
    }
}

fn periodic_environment(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    lattice: &[f64; 9],
    periodic: [bool; 3],
    model: Model,
) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    validate(numbers, positions)?;
    periodic_reciprocal(lattice, periodic)?;
    let coordinates: Vec<_> = positions
        .iter()
        .map(|&value| Dual::constant(value, 0))
        .collect();
    let dual_lattice: [Dual; 9] = std::array::from_fn(|index| Dual::constant(lattice[index], 0));
    let coordination = periodic_dual_coordination(
        numbers,
        &coordinates,
        &dual_lattice,
        lattice,
        periodic,
        false,
        model.cutoff.cn,
    );
    let charges = periodic_dual_charges::<false>(
        numbers,
        &coordinates,
        charge,
        &dual_lattice,
        lattice,
        periodic,
        model,
    )?;
    Ok((
        coordination.into_iter().map(|value| value.value).collect(),
        charges.into_iter().map(|value| value.value).collect(),
    ))
}

fn periodic_dual_coordination(
    numbers: &[i32],
    coordinates: &[Dual],
    lattice: &[Dual; 9],
    lattice_values: &[f64; 9],
    periodic: [bool; 3],
    charge_model: bool,
    cutoff: f64,
) -> Vec<Dual> {
    let size = coordinates[0].size;
    let images = lattice_indices(lattice_values, periodic, cutoff);
    let positions = positions_values(coordinates);
    let mut result = vec![Dual::constant(0.0, size); numbers.len()];
    for first in 0..numbers.len() {
        let first_element = numbers[first] as usize - 1;
        for second in 0..=first {
            let second_element = numbers[second] as usize - 1;
            let radius = if charge_model {
                real(CHARGE_RCOV, first_element) + real(CHARGE_RCOV, second_element)
            } else {
                element(first_element, 0) + element(second_element, 0)
            };
            let factor = if charge_model {
                1.0
            } else {
                let difference = (element(first_element, 3) - element(second_element, 3)).abs();
                4.10451 * (-(difference + 19.08857).powi(2) / (2.0 * 11.28174_f64.powi(2))).exp()
            };
            let mut value = 0.0;
            let mut partials = [0.0; 12];
            for &image in &images {
                let vector = image_displacement(&positions, lattice_values, first, second, image);
                let distance2: f64 = vector.iter().map(|value| value * value).sum();
                if distance2 > cutoff * cutoff || distance2 < 1.0e-12 {
                    continue;
                }
                let distance = distance2.sqrt();
                let argument = -7.5 * (distance - radius) / radius;
                value += 0.5 * (1.0 + unsafe { erf(argument) }) * factor;
                if size != 0 {
                    let radial = -7.5 / (radius * PI.sqrt() * distance)
                        * (-argument * argument).exp()
                        * factor;
                    for axis in 0..3 {
                        let derivative = radial * vector[axis];
                        partials[axis] += derivative;
                        for column in 0..3 {
                            partials[3 + axis + 3 * column] -= derivative * image[column] as f64;
                        }
                    }
                }
            }
            let contribution =
                periodic_pair_dual(value, partials, coordinates, lattice, first, second);
            result[first] += contribution.clone();
            if first != second {
                result[second] += contribution;
            }
        }
    }
    if charge_model {
        result
            .into_iter()
            .map(|value| {
                Dual::constant((1.0 + 8.0_f64.exp()).ln(), size)
                    - (Dual::constant(1.0, size) + (Dual::constant(8.0, size) - value).exp()).ln()
            })
            .collect()
    } else {
        result
    }
}

fn periodic_pair_dual(
    value: f64,
    partials: [f64; 12],
    coordinates: &[Dual],
    lattice: &[Dual; 9],
    first: usize,
    second: usize,
) -> Dual {
    let mut result = Dual::constant(value, coordinates[0].size);
    if result.size != 0 {
        for axis in 0..3 {
            if first != second {
                for (atom, sign) in [(first, 1.0), (second, -1.0)] {
                    result = Dual::record(
                        value,
                        result.size,
                        result.node.map(|node| (node, 1.0)),
                        coordinates[3 * atom + axis]
                            .node
                            .map(|node| (node, sign * partials[axis])),
                    );
                }
            }
        }
        for index in 0..9 {
            result = Dual::record(
                value,
                result.size,
                result.node.map(|node| (node, 1.0)),
                lattice[index].node.map(|node| (node, partials[3 + index])),
            );
        }
    }
    result
}

fn lattice_indices(lattice: &[f64; 9], periodic: [bool; 3], cutoff: f64) -> Vec<[i32; 3]> {
    let repetitions = lattice_repetitions(lattice, periodic, cutoff);
    let repetitions: [i32; 3] =
        std::array::from_fn(|axis| repetitions[axis] + i32::from(periodic[axis]));
    let mut result = Vec::new();
    for first in -repetitions[0]..=repetitions[0] {
        for second in -repetitions[1]..=repetitions[1] {
            for third in -repetitions[2]..=repetitions[2] {
                result.push([first, second, third]);
            }
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn partial_ewald_eeq_charges<const SECOND: bool>(
    numbers: &[i32],
    coordinates: &[Dual],
    charge: f64,
    lattice: &[Dual; 9],
    lattice_values: &[f64; 9],
    periodic: [bool; 3],
    cutoff: f64,
) -> Result<Vec<Dual>, &'static str> {
    if !charge.is_finite() {
        return Err("D4 total charge must be finite");
    }
    let atoms = numbers.len();
    let dimension = atoms + 1;
    let size = coordinates[0].size;
    let constant = |value| Dual::constant(value, size);
    let coordination_function = if SECOND {
        hessian::coordination
    } else {
        periodic_dual_coordination
    };
    let coordination = coordination_function(
        numbers,
        coordinates,
        lattice,
        lattice_values,
        periodic,
        true,
        25.0,
    );
    periodic_reciprocal(lattice_values, periodic)?;
    let shortest = (0..3)
        .filter(|&column| periodic[column])
        .map(|column| {
            (0..3)
                .map(|axis| lattice_values[axis + 3 * column].powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .fold(f64::INFINITY, f64::min);
    let alpha = (PI.sqrt() / shortest).max(7.0 / cutoff.max(shortest));
    let mut matrix = vec![constant(0.0); dimension * dimension];
    let mut rhs = vec![constant(0.0); dimension];
    for first in 0..atoms {
        let first_element = numbers[first] as usize - 1;
        rhs[first] = (coordination[first].clone() + 1.0e-14).sqrt() * eeq(first_element, 2)
            + -eeq(first_element, 0);
        matrix[first * dimension + atoms] = constant(1.0);
        matrix[atoms * dimension + first] = constant(1.0);
        for second in 0..=first {
            let second_element = numbers[second] as usize - 1;
            let gamma =
                1.0 / (eeq(first_element, 3).powi(2) + eeq(second_element, 3).powi(2)).sqrt();
            let displacement = std::array::from_fn(|axis| {
                coordinates[3 * first + axis].clone() - coordinates[3 * second + axis].clone()
            });
            let mut interaction =
                ewald::gaussian_pair(&displacement, lattice, periodic, constant(gamma), alpha);
            if first == second {
                interaction = interaction + eeq(first_element, 1);
            }
            matrix[first * dimension + second] = interaction.clone();
            matrix[second * dimension + first] = interaction;
        }
    }
    rhs[atoms] = constant(charge);
    dual_solve(&mut matrix, &mut rhs)?;
    rhs.pop();
    Ok(rhs)
}

fn periodic_dual_charges<const SECOND: bool>(
    numbers: &[i32],
    coordinates: &[Dual],
    charge: f64,
    lattice: &[Dual; 9],
    lattice_values: &[f64; 9],
    periodic: [bool; 3],
    model: Model,
) -> Result<Vec<Dual>, &'static str> {
    model.validate_selection(numbers.len())?;
    if let Some(charges) = model.fixed_charges {
        return Ok(charges
            .iter()
            .map(|&value| Dual::constant(value, coordinates[0].size))
            .collect());
    }
    if model.eeqbc {
        return eeqbc_charges(numbers, coordinates, charge, Some((lattice, periodic)));
    }
    if !periodic.iter().all(|&active| active) {
        return partial_ewald_eeq_charges::<SECOND>(
            numbers,
            coordinates,
            charge,
            lattice,
            lattice_values,
            periodic,
            model.charge_cutoff,
        );
    }
    let derivatives = coordinates[0].size;
    let atoms = numbers.len();
    let size = atoms + 1;
    let coordination_function = if SECOND {
        hessian::coordination
    } else {
        periodic_dual_coordination
    };
    let coordination = coordination_function(
        numbers,
        coordinates,
        lattice,
        lattice_values,
        periodic,
        true,
        25.0,
    );
    let mut matrix = vec![Dual::constant(0.0, derivatives); size * size];
    let mut rhs = vec![Dual::constant(0.0, derivatives); size];
    let alpha = ewald_alpha(lattice_values);
    let volume = determinant(lattice_values);
    let inverse = inverse(lattice_values);
    let reciprocal = std::array::from_fn(|index| inverse[(index % 3) * 3 + index / 3] * 2.0 * PI);
    let direct_indices = fixed_indices(2, true);
    let reciprocal_indices = fixed_indices(2, false);
    let position_values = positions_values(coordinates);
    let direct_vectors: Vec<_> = direct_indices
        .iter()
        .map(|&image| lattice_translation(lattice_values, image))
        .collect();
    let reciprocal_terms: Vec<_> = reciprocal_indices
        .iter()
        .map(|&image| {
            let wave = lattice_translation(&reciprocal, image);
            let wave2: f64 = wave.iter().map(|value| value * value).sum();
            let factor = (-wave2 / (4.0 * alpha * alpha)).exp() / wave2 / volume * (4.0 * PI);
            let partials: [f64; 9] = std::array::from_fn(|index| {
                let axis = index % 3;
                let column = index / 3;
                let transformed = (0..3)
                    .map(|component| inverse[column + 3 * component] * wave[component])
                    .sum::<f64>();
                factor
                    * (2.0 * wave[axis] * transformed * (1.0 / (4.0 * alpha * alpha) + 1.0 / wave2)
                        - inverse[column + 3 * axis])
            });
            (wave, factor, partials)
        })
        .collect();
    for first in 0..atoms {
        let first_element = numbers[first] as usize - 1;
        rhs[first] = (coordination[first].clone() + 1.0e-14).sqrt() * eeq(first_element, 2)
            + -eeq(first_element, 0);
        matrix[first * size + atoms] = Dual::constant(1.0, derivatives);
        matrix[atoms * size + first] = Dual::constant(1.0, derivatives);
        for second in 0..=first {
            let second_element = numbers[second] as usize - 1;
            let gamma =
                1.0 / (eeq(first_element, 3).powi(2) + eeq(second_element, 3).powi(2)).sqrt();
            if SECOND {
                let mut interaction =
                    hessian::eeq_pair(coordinates, lattice_values, first, second, gamma, alpha);
                if first == second {
                    interaction = interaction
                        + eeq(first_element, 1)
                        + (2.0 / PI).sqrt() / eeq(first_element, 3)
                        - Dual::constant(2.0 * alpha / PI.sqrt(), derivatives);
                }
                matrix[first * size + second] = interaction.clone();
                matrix[second * size + first] = interaction;
                continue;
            }
            let images = closest_images(&position_values, lattice_values, first, second);
            let mut value = 0.0;
            let mut partials = [0.0; 12];
            for image in &images {
                let base =
                    image_displacement(&position_values, lattice_values, first, second, *image);
                for (direct, index) in direct_vectors.iter().zip(&direct_indices) {
                    let vector: [f64; 3] = std::array::from_fn(|axis| base[axis] + direct[axis]);
                    let distance2: f64 = vector.iter().map(|value| value * value).sum();
                    let distance = distance2.sqrt();
                    if distance > f64::EPSILON.sqrt() {
                        let contribution =
                            (unsafe { erf(distance * gamma) - erf(distance * alpha) }) / distance;
                        value += contribution;
                        if derivatives != 0 {
                            let radial = (2.0 / PI.sqrt()
                                * (gamma * (-distance2 * gamma * gamma).exp()
                                    - alpha * (-distance2 * alpha * alpha).exp())
                                - contribution)
                                / distance2;
                            for axis in 0..3 {
                                let derivative = radial * vector[axis];
                                partials[axis] += derivative;
                                for column in 0..3 {
                                    partials[3 + axis + 3 * column] +=
                                        derivative * (index[column] - image[column]) as f64;
                                }
                            }
                        }
                    }
                }
                let fractional: [f64; 3] = std::array::from_fn(|column| {
                    image[column] as f64
                        + (0..3)
                            .map(|axis| inverse[column + 3 * axis] * base[axis])
                            .sum::<f64>()
                });
                for (wave, factor, factor_partials) in &reciprocal_terms {
                    let phase: f64 = (0..3).map(|axis| base[axis] * wave[axis]).sum();
                    let cosine = phase.cos();
                    value += cosine * factor;
                    if derivatives != 0 {
                        let sine = phase.sin() * factor;
                        for axis in 0..3 {
                            partials[axis] -= sine * wave[axis];
                            for column in 0..3 {
                                partials[3 + axis + 3 * column] += cosine
                                    * factor_partials[axis + 3 * column]
                                    + sine * wave[axis] * fractional[column];
                            }
                        }
                    }
                }
            }
            value /= images.len() as f64;
            partials
                .iter_mut()
                .for_each(|partial| *partial /= images.len() as f64);
            if first == second {
                value += eeq(first_element, 1) + (2.0 / PI).sqrt() / eeq(first_element, 3)
                    - 2.0 * alpha / PI.sqrt();
            }
            let interaction =
                periodic_pair_dual(value, partials, coordinates, lattice, first, second);
            matrix[first * size + second] = interaction.clone();
            matrix[second * size + first] = interaction;
        }
    }
    rhs[atoms] = Dual::constant(charge, derivatives);
    dual_solve(&mut matrix, &mut rhs)?;
    rhs.pop();
    Ok(rhs)
}

fn positions_values(coordinates: &[Dual]) -> Vec<f64> {
    coordinates.iter().map(|value| value.value).collect()
}

fn fixed_indices(repetitions: i32, origin: bool) -> Vec<[i32; 3]> {
    let mut result = Vec::new();
    for first in -repetitions..=repetitions {
        for second in -repetitions..=repetitions {
            for third in -repetitions..=repetitions {
                if origin || first != 0 || second != 0 || third != 0 {
                    result.push([first, second, third]);
                }
            }
        }
    }
    result
}

fn closest_images(
    positions: &[f64],
    lattice: &[f64; 9],
    first: usize,
    second: usize,
) -> Vec<[i32; 3]> {
    closest_directional_images(positions, lattice, first, second, [true; 3])
}

fn closest_directional_images(
    positions: &[f64],
    lattice: &[f64; 9],
    first: usize,
    second: usize,
    periodic: [bool; 3],
) -> Vec<[i32; 3]> {
    let candidates = fixed_indices(1, true);
    let mut distances: Vec<_> = candidates
        .into_iter()
        .filter(|image| (0..3).all(|axis| periodic[axis] || image[axis] == 0))
        .filter_map(|image| {
            let distance2 = (0..3)
                .map(|axis| {
                    positions[3 * first + axis]
                        - positions[3 * second + axis]
                        - (0..3)
                            .map(|column| lattice[axis + 3 * column] * image[column] as f64)
                            .sum::<f64>()
                })
                .map(|value| value * value)
                .sum::<f64>();
            (distance2 >= f64::EPSILON.sqrt()).then_some((image, distance2))
        })
        .collect();
    distances.sort_by(|left, right| left.1.total_cmp(&right.1));
    let minimum = distances[0].1;
    distances
        .into_iter()
        .take_while(|(_, distance)| (distance - minimum).abs() <= 0.01)
        .map(|(image, _)| image)
        .collect()
}

fn ewald_alpha(lattice: &[f64; 9]) -> f64 {
    let volume = determinant(lattice).abs();
    let inverse = inverse(lattice);
    let reciprocal: [f64; 9] =
        std::array::from_fn(|index| inverse[(index % 3) * 3 + index / 3] * 2.0 * PI);
    let shortest = |matrix: &[f64; 9]| {
        (0..3)
            .map(|column| {
                (0..3)
                    .map(|axis| matrix[axis + 3 * column].powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .fold(f64::INFINITY, f64::min)
    };
    let direct = shortest(lattice);
    let wave = shortest(&reciprocal);
    let difference = |alpha: f64| {
        let reciprocal_term = |length: f64| {
            4.0 * PI * (-0.25 * length * length / alpha.powi(2)).exp() / (volume * length * length)
        };
        let direct_term = |length: f64| 1.0 - unsafe { erf(alpha * length) };
        reciprocal_term(4.0 * wave)
            - reciprocal_term(5.0 * wave)
            - direct_term(2.0 * direct) / (2.0 * direct)
            + direct_term(3.0 * direct) / (3.0 * direct)
    };
    let tolerance = f64::EPSILON.sqrt();
    let mut alpha = 1.0e-8;
    let mut diff = difference(alpha);
    while diff < -tolerance {
        alpha *= 2.0;
        diff = difference(alpha);
    }
    let mut left = 0.5 * alpha;
    while diff < tolerance {
        alpha *= 2.0;
        diff = difference(alpha);
    }
    let mut right = alpha;
    alpha = 0.5 * (left + right);
    diff = difference(alpha);
    let mut iterations = 0;
    while diff.abs() > tolerance && iterations <= 30 {
        if diff < 0.0 {
            left = alpha;
        } else {
            right = alpha;
        }
        alpha = 0.5 * (left + right);
        diff = difference(alpha);
        iterations += 1;
    }
    alpha
}

fn periodic_differentiated_energy(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    lattice_values: &[f64; 9],
    periodic: [bool; 3],
) -> Result<Dual, &'static str> {
    TAPE.with(|tape| unsafe {
        let tape = &mut *tape.get();
        tape.nodes.clear();
        tape.solves.clear();
    });
    validate(numbers, positions)?;
    model.validate_selection(numbers.len())?;
    periodic_reciprocal(lattice_values, periodic)?;
    let coordinate_count = positions.len();
    let size = coordinate_count + 9;
    let coordinates: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(index, &value)| Dual::variable(value, size, index))
        .collect();
    let lattice: [Dual; 9] = std::array::from_fn(|index| {
        Dual::variable(lattice_values[index], size, coordinate_count + index)
    });
    periodic_dual_energy(
        numbers,
        &coordinates,
        charge,
        model,
        param,
        &lattice,
        lattice_values,
        periodic,
    )
}

#[allow(clippy::type_complexity)]
fn periodic_scalar_energy<const PAIRWISE: bool>(
    numbers: &[i32],
    positions: &[f64],
    charge: f64,
    model: Model,
    param: Param,
    lattice_values: &[f64; 9],
    periodic: [bool; 3],
) -> Result<(f64, Vec<f64>, Vec<f64>), &'static str> {
    let cutoff = model.cutoff;
    validate(numbers, positions)?;
    model.validate_selection(numbers.len())?;
    periodic_reciprocal(lattice_values, periodic)?;
    let (coordination, charges) =
        periodic_environment(numbers, positions, charge, lattice_values, periodic, model)?;
    let c6 = model_coefficients(numbers, &coordination, &charges, model);
    let mut energy = 0.0;
    let mut pair2 = vec![
        0.0;
        if PAIRWISE {
            numbers.len() * numbers.len()
        } else {
            0
        }
    ];
    let mut pair3 = pair2.clone();
    let pair_images = lattice_indices(lattice_values, periodic, cutoff.disp2);
    for first in 0..numbers.len() {
        for second in 0..=first {
            if !model.owns_pair(first, second) {
                continue;
            }
            for &image in &pair_images {
                let vector = image_displacement(positions, lattice_values, first, second, image);
                let distance2 = vector.iter().map(|value| value * value).sum::<f64>();
                if distance2 > cutoff.disp2.powi(2) || distance2 < f64::EPSILON {
                    continue;
                }
                let contribution = pair_energy(
                    numbers[first] as usize - 1,
                    numbers[second] as usize - 1,
                    distance2,
                    c6[first * numbers.len() + second],
                    param,
                ) * smooth_cutoff(distance2, cutoff.disp2, cutoff.width2).0
                    * if first == second { 0.5 } else { 1.0 };
                energy += contribution;
                if PAIRWISE {
                    pair2[first * numbers.len() + second] += contribution * 0.5;
                    pair2[second * numbers.len() + first] += contribution * 0.5;
                }
            }
        }
    }
    if param.s9.abs() >= f64::EPSILON {
        let c6 = model_coefficients(numbers, &coordination, &vec![0.0; numbers.len()], model);
        let images = lattice_indices(lattice_values, periodic, cutoff.disp3);
        let contributions = crate::parallel::map(numbers.len(), 8, |start, stride| {
            let mut energy = 0.0;
            let mut pair3 = vec![
                0.0;
                if PAIRWISE {
                    numbers.len() * numbers.len()
                } else {
                    0
                }
            ];
            for first in (start..numbers.len()).step_by(stride) {
                if !model.active(first) {
                    continue;
                }
                let vectors: Vec<Vec<([f64; 3], f64)>> = (0..=first)
                    .map(|atom| {
                        images
                            .iter()
                            .filter_map(|&image| {
                                let vector = image_displacement(
                                    positions,
                                    lattice_values,
                                    atom,
                                    first,
                                    image.map(|value| -value),
                                );
                                let distance =
                                    vector.iter().map(|value| value * value).sum::<f64>();
                                (distance <= cutoff.disp3.powi(2) && distance >= f64::EPSILON)
                                    .then_some((vector, distance))
                            })
                            .collect()
                    })
                    .collect();
                for second in 0..=first {
                    if !model.owns_pair(first, second) {
                        continue;
                    }
                    for (third, third_vectors) in vectors.iter().enumerate().take(second + 1) {
                        if !model.active(third) {
                            continue;
                        }
                        let radius = [(first, second), (first, third), (second, third)]
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
                        let triple = if first == second {
                            if first == third {
                                1.0 / 6.0
                            } else {
                                0.5
                            }
                        } else if first != third && second != third {
                            1.0
                        } else {
                            0.5
                        };
                        for (first_second, distance_ij) in &vectors[second] {
                            for (first_third, distance_ik) in third_vectors {
                                let distance_jk = (0..3)
                                    .map(|axis| first_third[axis] - first_second[axis])
                                    .map(|value| value * value)
                                    .sum::<f64>();
                                if distance_jk > cutoff.disp3.powi(2) || distance_jk < f64::EPSILON
                                {
                                    continue;
                                }
                                let contribution = scalar_atm(
                                    [*distance_ij, *distance_ik, distance_jk],
                                    [
                                        c6[first * numbers.len() + second],
                                        c6[first * numbers.len() + third],
                                        c6[second * numbers.len() + third],
                                    ],
                                    radius,
                                    param,
                                ) * triple
                                    * [*distance_ij, *distance_ik, distance_jk]
                                        .map(|distance| {
                                            smooth_cutoff(distance, cutoff.disp3, cutoff.width3).0
                                        })
                                        .iter()
                                        .product::<f64>();
                                energy += contribution;
                                if PAIRWISE {
                                    for (left, right) in
                                        [(first, second), (first, third), (second, third)]
                                    {
                                        pair3[left * numbers.len() + right] += contribution / 6.0;
                                        pair3[right * numbers.len() + left] += contribution / 6.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            (energy, pair3)
        });
        for (contribution, pairs) in contributions {
            energy += contribution;
            for (total, value) in pair3.iter_mut().zip(pairs) {
                *total += value;
            }
        }
    }
    Ok((energy, pair2, pair3))
}

fn scalar_atm(distance: [f64; 3], coefficient: [f64; 3], radius: f64, param: Param) -> f64 {
    let product2 = distance[0] * distance[1] * distance[2];
    let product = product2.sqrt();
    let c9 = -param.s9 * (coefficient[0] * coefficient[1] * coefficient[2]).sqrt();
    let damping = 1.0 / (1.0 + 6.0 * (radius / product).powf(param.alpha / 3.0));
    let angular = 0.375
        * (distance[0] + distance[2] - distance[1])
        * (distance[0] - distance[2] + distance[1])
        * (-distance[0] + distance[2] + distance[1])
        / (product2 * product2 * product)
        + 1.0 / (product2 * product);
    -c9 * damping * angular
}

fn image_displacement(
    positions: &[f64],
    lattice: &[f64; 9],
    first: usize,
    second: usize,
    image: [i32; 3],
) -> [f64; 3] {
    let translation = lattice_translation(lattice, image);
    std::array::from_fn(|axis| {
        positions[3 * first + axis] - positions[3 * second + axis] - translation[axis]
    })
}

#[allow(clippy::too_many_arguments)]
fn periodic_dual_energy(
    numbers: &[i32],
    coordinates: &[Dual],
    charge: f64,
    model: Model,
    param: Param,
    lattice: &[Dual; 9],
    lattice_values: &[f64; 9],
    periodic: [bool; 3],
) -> Result<Dual, &'static str> {
    let cutoff = model.cutoff;
    let size = coordinates[0].size;
    let coordination = periodic_dual_coordination(
        numbers,
        coordinates,
        lattice,
        lattice_values,
        periodic,
        false,
        cutoff.cn,
    );
    let charges = periodic_dual_charges::<false>(
        numbers,
        coordinates,
        charge,
        lattice,
        lattice_values,
        periodic,
        model,
    )?;
    let c6 = dual_coefficients(numbers, &coordination, &charges, model);
    let pair_images = lattice_indices(lattice_values, periodic, cutoff.disp2);
    let mut total = Dual::constant(0.0, size);
    let positions = positions_values(coordinates);
    let atoms = numbers.len();
    let mut gradient = vec![0.0; size];
    for first in 0..numbers.len() {
        for second in 0..=first {
            if !model.owns_pair(first, second) {
                continue;
            }
            let rrij = 3.0
                * element(numbers[first] as usize - 1, 1)
                * element(numbers[second] as usize - 1, 1);
            let radius = param.a1 * rrij.sqrt() + param.a2;
            let weight = if first == second { 0.5 } else { 1.0 };
            let mut potential = 0.0;
            for &image in &pair_images {
                let vector = image_displacement(&positions, lattice_values, first, second, image);
                let distance2: f64 = vector.iter().map(|value| value * value).sum();
                if distance2 > cutoff.disp2.powi(2) || distance2 < f64::EPSILON {
                    continue;
                }
                let sixth = 1.0 / (distance2.powi(3) + radius.powi(6));
                let eighth = 1.0 / (distance2.powi(4) + radius.powi(8));
                let switch = smooth_cutoff(distance2, cutoff.disp2, cutoff.width2);
                potential -= weight * (param.s6 * sixth + param.s8 * rrij * eighth) * switch.0;
                add_image_gradient(
                    &mut gradient,
                    first,
                    second,
                    vector,
                    image,
                    weight
                        * c6[first * atoms + second].value
                        * (switch.0
                            * (3.0 * param.s6 * distance2.powi(2) * sixth.powi(2)
                                + 4.0 * param.s8 * rrij * distance2.powi(3) * eighth.powi(2))
                            - (param.s6 * sixth + param.s8 * rrij * eighth) * switch.1),
                );
            }
            total += c6[first * atoms + second].clone() * potential;
        }
    }
    if param.s9.abs() >= f64::EPSILON {
        let zero = vec![Dual::constant(0.0, size); numbers.len()];
        let c6 = dual_coefficients(numbers, &coordination, &zero, model);
        let images = lattice_indices(lattice_values, periodic, cutoff.disp3);
        let contributions = crate::parallel::map(atoms, 8, |start, stride| {
            let mut gradient = vec![0.0; size];
            let mut coefficient_derivatives = vec![0.0; atoms * atoms];
            let mut energy = 0.0;
            for first in (start..atoms).step_by(stride) {
                if !model.active(first) {
                    continue;
                }
                let vectors: Vec<Vec<_>> = (0..=first)
                    .map(|atom| {
                        images
                            .iter()
                            .filter_map(|&image| {
                                let vector = image_displacement(
                                    &positions,
                                    lattice_values,
                                    first,
                                    atom,
                                    image,
                                );
                                let distance =
                                    vector.iter().map(|value| value * value).sum::<f64>();
                                (distance <= cutoff.disp3.powi(2) && distance >= f64::EPSILON)
                                    .then_some((vector, distance, image))
                            })
                            .collect()
                    })
                    .collect();
                for second in 0..=first {
                    if !model.owns_pair(first, second) {
                        continue;
                    }
                    for (third, third_vectors) in vectors.iter().enumerate().take(second + 1) {
                        if !model.active(third) {
                            continue;
                        }
                        let pairs = [(first, second), (first, third), (second, third)];
                        let radius = pairs
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
                        let coefficients =
                            pairs.map(|(left, right)| c6[left * atoms + right].value);
                        let weight = if first == third {
                            1.0 / 6.0
                        } else if first == second || second == third {
                            0.5
                        } else {
                            1.0
                        };
                        for (first_second, distance_ij, image_ij) in &vectors[second] {
                            for (first_third, distance_ik, image_ik) in third_vectors {
                                let second_third: [f64; 3] = std::array::from_fn(|axis| {
                                    first_third[axis] - first_second[axis]
                                });
                                let distance_jk: f64 =
                                    second_third.iter().map(|value| value * value).sum();
                                if distance_jk > cutoff.disp3.powi(2) || distance_jk < f64::EPSILON
                                {
                                    continue;
                                }
                                let distances = [*distance_ij, *distance_ik, distance_jk];
                                let local = smooth_atm(
                                    local_atm(
                                        [*distance_ij, *distance_ik, distance_jk],
                                        coefficients,
                                        radius,
                                        param,
                                    ),
                                    distances,
                                    cutoff,
                                )
                                .scale(weight);
                                energy += local.value;
                                let vectors = [*first_second, *first_third, second_third];
                                let images = [
                                    *image_ij,
                                    *image_ik,
                                    std::array::from_fn(|axis| image_ik[axis] - image_ij[axis]),
                                ];
                                for (index, &(left, right)) in pairs.iter().enumerate() {
                                    add_image_gradient(
                                        &mut gradient,
                                        left,
                                        right,
                                        vectors[index],
                                        images[index],
                                        local.gradient[index],
                                    );
                                    coefficient_derivatives[left * atoms + right] +=
                                        local.gradient[3 + index];
                                }
                            }
                        }
                    }
                }
            }
            (energy, gradient, coefficient_derivatives)
        });
        for (energy, local_gradient, coefficient_derivatives) in contributions {
            for (total, value) in gradient.iter_mut().zip(local_gradient) {
                *total += value;
            }
            total.value += energy;
            for (coefficient, derivative) in c6.iter().zip(coefficient_derivatives) {
                total = Dual::record(
                    total.value,
                    size,
                    total.node.map(|node| (node, 1.0)),
                    coefficient.node.map(|node| (node, derivative)),
                );
            }
        }
    }
    for (input, derivative) in coordinates.iter().chain(lattice).zip(gradient) {
        total = Dual::record(
            total.value,
            size,
            total.node.map(|node| (node, 1.0)),
            input.node.map(|node| (node, derivative)),
        );
    }
    Ok(total)
}

fn add_image_gradient(
    gradient: &mut [f64],
    first: usize,
    second: usize,
    vector: [f64; 3],
    image: [i32; 3],
    derivative: f64,
) {
    let coordinates = gradient.len() - 9;
    for axis in 0..3 {
        let contribution = 2.0 * derivative * vector[axis];
        gradient[3 * first + axis] += contribution;
        gradient[3 * second + axis] -= contribution;
        for column in 0..3 {
            gradient[coordinates + axis + 3 * column] -= contribution * image[column] as f64;
        }
    }
}

fn dual_solve(matrix: &mut [Dual], rhs: &mut [Dual]) -> Result<(), &'static str> {
    let derivatives = rhs[0].size;
    let mut matrix_values: Vec<_> = matrix.iter().map(|value| value.value).collect();
    if derivatives == 0 {
        let mut values: Vec<_> = rhs.iter().map(|value| value.value).collect();
        solve(&mut matrix_values, &mut values)?;
        for (output, value) in rhs.iter_mut().zip(values) {
            *output = Dual::constant(value, 0);
        }
        return Ok(());
    }
    let values = matrix_values.clone();
    let mut solution: Vec<_> = rhs.iter().map(|value| value.value).collect();
    solve(&mut matrix_values, &mut solution)?;
    let matrix_nodes = matrix.iter().map(|value| value.node).collect();
    let rhs_nodes = rhs.iter().map(|value| value.node).collect();
    TAPE.with(|tape| {
        let tape = unsafe { &mut *tape.get() };
        let operation = tape.solves.len();
        tape.solves.push(SolveNode {
            matrix: matrix_nodes,
            rhs: rhs_nodes,
            values,
            solution: solution.clone(),
        });
        for (output_index, (output, value)) in rhs.iter_mut().zip(solution).enumerate() {
            let node = tape.nodes.len();
            tape.nodes.push(TapeNode {
                parents: [(0, 0.0); 2],
                curvature: [[0.0; 2]; 2],
                parent_count: 0,
                input: None,
                solve: Some((operation, output_index)),
            });
            *output = Dual {
                value,
                node: Some(node),
                size: derivatives,
            };
        }
    });
    Ok(())
}

fn dual_coefficients(
    numbers: &[i32],
    coordination: &[Dual],
    charges: &[Dual],
    model: Model,
) -> Vec<Dual> {
    let size = coordination[0].size;
    let atoms = numbers.len();
    let mut result = vec![Dual::constant(0.0, size); atoms * atoms];
    let partials = model_coefficient_partials(
        numbers,
        &positions_values(coordination),
        &positions_values(charges),
        model,
    );
    for first in 0..atoms {
        for second in 0..=first {
            let partial = partials[first * atoms + second];
            let cn = Dual::record(
                partial[0],
                size,
                coordination[first].node.map(|node| (node, partial[1])),
                coordination[second].node.map(|node| (node, partial[2])),
            );
            let charge = Dual::record(
                0.0,
                size,
                charges[first].node.map(|node| (node, partial[3])),
                charges[second].node.map(|node| (node, partial[4])),
            );
            let coefficient = cn + charge;
            result[first * atoms + second] = coefficient.clone();
            result[second * atoms + first] = coefficient;
        }
    }
    result
}

#[derive(Clone, Copy)]
pub(crate) struct Local {
    pub(crate) value: f64,
    pub(crate) gradient: [f64; 6],
}

impl Local {
    fn scale(self, factor: f64) -> Self {
        Self {
            value: self.value * factor,
            gradient: self.gradient.map(|value| value * factor),
        }
    }
}

pub(crate) fn local_atm(
    distances: [f64; 3],
    coefficients: [f64; 3],
    radius: f64,
    param: Param,
) -> Local {
    let product2 = distances.iter().product::<f64>();
    let product = product2.sqrt();
    let inverse3 = 1.0 / (product2 * product);
    let inverse5 = inverse3 / product2;
    let factors = [
        distances[0] + distances[2] - distances[1],
        distances[0] - distances[2] + distances[1],
        -distances[0] + distances[2] + distances[1],
    ];
    let angular3 = 0.375 * factors.iter().product::<f64>() * inverse5;
    let angular = angular3 + inverse3;
    let c9 = param.s9 * coefficients.iter().product::<f64>().sqrt();
    let exponent = param.alpha / 3.0;
    let damping = 1.0 / (1.0 + 6.0 * (radius / product).powf(exponent));
    let value = c9 * damping * angular;
    let mut gradient = [0.0; 6];
    let factor_derivatives = [[1.0, 1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, -1.0, 1.0]];
    for index in 0..3 {
        let signs = factor_derivatives[index];
        let angular_derivative = 0.375
            * inverse5
            * (signs[0] * factors[1] * factors[2]
                + signs[1] * factors[0] * factors[2]
                + signs[2] * factors[0] * factors[1])
            - (2.5 * angular3 + 1.5 * inverse3) / distances[index];
        gradient[index] = c9
            * damping
            * (angular_derivative + angular * 0.5 * exponent * (1.0 - damping) / distances[index]);
        gradient[3 + index] = 0.5 * value / coefficients[index];
    }
    Local { value, gradient }
}
