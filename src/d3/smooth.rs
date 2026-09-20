use super::*;

const ELEMENTS: usize = 94;

/// C6 interpolation model. D3S uses the pair widths published with
/// Tkachenko and Head-Gordon, DOI: 10.1021/acs.jctc.4c01105 (SI 004).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Model {
    #[default]
    D3,
    D3S,
}

impl Model {
    /// Coordination numbers and the row-major C6 matrix.
    pub fn properties(
        self,
        numbers: &[i32],
        positions: &[f64],
        cell: Option<(&[f64; 9], [bool; 3])>,
        cn_cutoff: f64,
    ) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
        let (lattice, periodic) = cell.map_or((None, [false; 3]), |(lattice, periodic)| {
            (Some(lattice), periodic)
        });
        properties_with_model(numbers, positions, lattice, periodic, cn_cutoff, self)
    }

    /// Property-major Cartesian and strain responses, as in `property_response`.
    pub fn property_response(
        self,
        numbers: &[i32],
        positions: &[f64],
        cell: Option<(&[f64; 9], [bool; 3])>,
        cn_cutoff: f64,
    ) -> Result<PropertyResponse, &'static str> {
        let (lattice, periodic) = cell.map_or((None, [false; 3]), |(lattice, periodic)| {
            (Some(lattice), periodic)
        });
        property_response_with_model(numbers, positions, lattice, periodic, cn_cutoff, self)
    }

    pub(super) fn validate(self, numbers: &[i32]) -> Result<(), &'static str> {
        if self == Self::D3S && numbers.iter().any(|number| !(1..=94).contains(number)) {
            return Err("D3S supports atomic numbers 1 through 94");
        }
        Ok(())
    }

    /// Energy, Cartesian gradient, virial, and two-/three-body pair matrices.
    /// Coordinates and cutoffs are in bohr. None selects a nonperiodic molecule.
    /// D3S retains the chosen D3 damping; ATM is an optional model extension.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn dispersion(
        self,
        numbers: &[i32],
        positions: &[f64],
        damping: Damping,
        atm: Option<Atm>,
        cutoff: RealspaceCutoff,
        cell: Option<(&[f64; 9], [bool; 3])>,
        ghosts: &[bool],
        partition: WorkPartition,
    ) -> Result<(PeriodicResult, Vec<f64>, Vec<f64>), &'static str> {
        let (lattice, periodic) = cell.unwrap_or((&[0.0; 9], [false; 3]));
        realspace::evaluate(
            numbers, positions, lattice, periodic, damping, atm, cutoff, ghosts, partition, None,
            self,
        )
    }

    /// Energy and Cartesian Hessian at fixed lattice, in row-major order.
    #[allow(clippy::too_many_arguments)]
    pub fn hessian(
        self,
        numbers: &[i32],
        positions: &[f64],
        damping: Damping,
        atm: Option<Atm>,
        cutoff: RealspaceCutoff,
        cell: Option<(&[f64; 9], [bool; 3])>,
        ghosts: &[bool],
        partition: WorkPartition,
    ) -> Result<(f64, Vec<f64>), &'static str> {
        let (lattice, periodic) = cell.unwrap_or((&[0.0; 9], [false; 3]));
        let mut hessian = Vec::new();
        let result = realspace::evaluate(
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
            self,
        )?;
        Ok((result.0.energy, hessian))
    }
}

fn widths() -> &'static [f64] {
    static WIDTHS: std::sync::OnceLock<Vec<f64>> = std::sync::OnceLock::new();
    WIDTHS.get_or_init(|| {
        include_str!("../../assets/d3s_weights.txt")
            .split_whitespace()
            .map(|value| value.parse().expect("invalid embedded D3S width"))
            .collect()
    })
}

pub(super) fn pair_coefficients(
    first: usize,
    second: usize,
    cn_first: f64,
    cn_second: f64,
) -> [f64; 6] {
    coefficients(
        first,
        second,
        cn_first,
        cn_second,
        widths()[first * ELEMENTS + second],
    )
}

fn interpolation(element: usize, coordination: f64, width: f64) -> [[f64; REFERENCES]; 3] {
    let count = reference_count(element);
    let offsets: [f64; REFERENCES] = std::array::from_fn(|reference| {
        value(REFERENCE_CN, element * REFERENCES + reference) - coordination
    });
    let nearest = offsets[..count]
        .iter()
        .map(|offset| offset * offset)
        .fold(f64::INFINITY, f64::min);
    let mut weights = [0.0; REFERENCES];
    for reference in 0..count {
        weights[reference] = (-width * (offsets[reference].powi(2) - nearest)).exp();
    }
    let norm = weights.iter().sum::<f64>();
    weights.iter_mut().for_each(|weight| *weight /= norm);
    let slopes = offsets.map(|offset| 2.0 * width * offset);
    let mean = weights
        .iter()
        .zip(slopes)
        .map(|(weight, slope)| weight * slope)
        .sum::<f64>();
    let variance = weights
        .iter()
        .zip(slopes)
        .map(|(weight, slope)| weight * (slope - mean).powi(2))
        .sum::<f64>();
    let first = std::array::from_fn(|reference| weights[reference] * (slopes[reference] - mean));
    let second = std::array::from_fn(|reference| {
        weights[reference] * ((slopes[reference] - mean).powi(2) - variance)
    });
    [weights, first, second]
}

fn coefficients(
    first: usize,
    second: usize,
    cn_first: f64,
    cn_second: f64,
    width: f64,
) -> [f64; 6] {
    let left = interpolation(first, cn_first, width);
    let right = interpolation(second, cn_second, width);
    atomic_c6_second_derivatives(
        first, second, &left[0], &right[0], &left[1], &right[1], &left[2], &right[2],
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn virial_ghosts_and_partitions() {
        let numbers = [1, 28, 1];
        let positions = [0.0, 0.0, 0.0, 3.0, 0.2, 0.0, 0.3, 1.4, 0.1];
        let lattice = [7.0, 0.1, 0.0, 0.2, 8.0, 0.0, 0.0, 0.0, 9.0];
        let damping = load_named("blyp", 1, false).unwrap();
        let cutoff = RealspaceCutoff {
            cn: 8.0,
            disp2: 9.0,
            disp3: 5.0,
            width2: 1.0,
            width3: 1.0,
        };
        for periodic in [[false; 3], [true; 3]] {
            for ghosts in [[false; 3], [true, false, false]] {
                let evaluate = |positions: &[f64], lattice: &[f64; 9], partition| {
                    Model::D3S
                        .dispersion(
                            &numbers,
                            positions,
                            damping,
                            Some(Atm {
                                s9: 1.0,
                                alpha: 16.0,
                            }),
                            cutoff,
                            Some((lattice, periodic)),
                            &ghosts,
                            partition,
                        )
                        .unwrap()
                };
                let (result, pair2, pair3) = evaluate(&positions, &lattice, WorkPartition::SERIAL);
                if ghosts[0] {
                    assert!(pair2[..3]
                        .iter()
                        .chain(&pair3[..3])
                        .all(|value| *value == 0.0));
                }
                let left = evaluate(&positions, &lattice, WorkPartition::new(0, 2).unwrap()).0;
                let right = evaluate(&positions, &lattice, WorkPartition::new(1, 2).unwrap()).0;
                assert!((result.energy - left.energy - right.energy).abs() < 1e-13);
                for ((total, left), right) in result
                    .gradient
                    .iter()
                    .chain(&result.virial)
                    .zip(left.gradient.iter().chain(&left.virial))
                    .zip(right.gradient.iter().chain(&right.virial))
                {
                    assert!((total - left - right).abs() < 1e-12);
                }
                let step = 1e-5;
                for row in 0..3 {
                    for column in 0..3 {
                        let deform = |values: &[f64], step: f64| {
                            let mut result = values.to_vec();
                            for (target, source) in result
                                .as_chunks_mut::<3>()
                                .0
                                .iter_mut()
                                .zip(values.as_chunks::<3>().0)
                            {
                                target[row] += step * source[column];
                            }
                            result
                        };
                        let plus = evaluate(
                            &deform(&positions, step),
                            &deform(&lattice, step).try_into().unwrap(),
                            WorkPartition::SERIAL,
                        )
                        .0
                        .energy;
                        let minus = evaluate(
                            &deform(&positions, -step),
                            &deform(&lattice, -step).try_into().unwrap(),
                            WorkPartition::SERIAL,
                        )
                        .0
                        .energy;
                        assert!(
                            ((plus - minus) / (2.0 * step) - result.virial[row + 3 * column]).abs()
                                < 1e-8
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn published_blyp_scan() {
        let mut lines = include_str!("../../tests/reference/d3s_scan.xyz")
            .lines()
            .filter(|line| !line.trim().is_empty() && line.trim() != ">");
        let damping = load_named("blyp", 1, false).unwrap();
        let mut count = 0;
        for expected in
            include_str!("../../tests/reference/d3s_scan_energies.txt").split_whitespace()
        {
            let atoms: usize = lines.next().unwrap().trim().parse().unwrap();
            lines.next().unwrap();
            let mut numbers = Vec::new();
            let mut positions = Vec::new();
            for _ in 0..atoms {
                let mut fields = lines.next().unwrap().split_whitespace();
                numbers.push(match fields.next().unwrap() {
                    "H" => 1,
                    "C" => 6,
                    "N" => 7,
                    "Cl" => 17,
                    "Ni" => 28,
                    "Zn" => 30,
                    symbol => panic!("unexpected reference element {symbol}"),
                });
                positions.extend(fields.map(|value| value.parse::<f64>().unwrap() / 0.52917726));
            }
            let energy = Model::D3S
                .dispersion(
                    &numbers,
                    &positions,
                    damping,
                    None,
                    RealspaceCutoff::default(),
                    None,
                    &[],
                    WorkPartition::SERIAL,
                )
                .unwrap()
                .0
                .energy;
            let expected: f64 = expected.parse().unwrap();
            assert!(
                (energy - expected).abs() < 1e-7,
                "frame {count}: native={energy:.14}, reference={expected:.14}"
            );
            count += 1;
        }
        assert_eq!(count, 41);
        assert!(lines.next().is_none());
    }

    #[test]
    fn molecular_and_periodic_derivatives() {
        let numbers = [1, 28, 1];
        let positions = [0.0, 0.0, 0.0, 3.0, 0.2, 0.0, 0.3, 1.4, 0.1];
        let lattice = [7.0, 0.1, 0.0, 0.2, 8.0, 0.0, 0.0, 0.0, 9.0];
        let cutoff = RealspaceCutoff {
            cn: 8.0,
            disp2: 9.0,
            disp3: 5.0,
            width2: 1.0,
            width3: 1.0,
        };
        for periodic in [[false; 3], [true, false, false], [true; 3]] {
            let cell = Some((&lattice, periodic));
            for kind in 0..=6 {
                let damping = if kind == 6 {
                    Damping::Z {
                        s6: 1.0,
                        s8: 1.0,
                        a1: 200770.0,
                    }
                } else {
                    load_named("pbe", kind, false).unwrap()
                };
                let atm = Some(Atm {
                    s9: 1.0,
                    alpha: 16.0,
                });
                let evaluate = |positions: &[f64]| {
                    Model::D3S
                        .dispersion(
                            &numbers,
                            positions,
                            damping,
                            atm,
                            cutoff,
                            cell,
                            &[],
                            WorkPartition::SERIAL,
                        )
                        .unwrap()
                };
                let (result, pair2, pair3) = evaluate(&positions);
                assert!((result.energy - pair2.iter().chain(&pair3).sum::<f64>()).abs() < 1e-13);
                let (energy, hessian) = Model::D3S
                    .hessian(
                        &numbers,
                        &positions,
                        damping,
                        atm,
                        cutoff,
                        cell,
                        &[],
                        WorkPartition::SERIAL,
                    )
                    .unwrap();
                assert!((energy - result.energy).abs() < 1e-13);
                let step = 1e-5;
                for coordinate in 0..positions.len() {
                    let mut plus = positions;
                    let mut minus = positions;
                    plus[coordinate] += step;
                    minus[coordinate] -= step;
                    let plus = evaluate(&plus).0;
                    let minus = evaluate(&minus).0;
                    assert!(
                        ((plus.energy - minus.energy) / (2.0 * step) - result.gradient[coordinate])
                            .abs()
                            < 1e-8
                    );
                    for row in 0..positions.len() {
                        assert!(
                            ((plus.gradient[row] - minus.gradient[row]) / (2.0 * step)
                                - hessian[row * positions.len() + coordinate])
                                .abs()
                                < 1e-7
                        );
                    }
                }
            }
            let response = Model::D3S
                .property_response(&numbers, &positions, cell, cutoff.cn)
                .unwrap();
            let step = 1e-5;
            for coordinate in 0..positions.len() {
                let mut plus = positions;
                let mut minus = positions;
                plus[coordinate] += step;
                minus[coordinate] -= step;
                let plus = Model::D3S
                    .properties(&numbers, &plus, cell, cutoff.cn)
                    .unwrap()
                    .1;
                let minus = Model::D3S
                    .properties(&numbers, &minus, cell, cutoff.cn)
                    .unwrap()
                    .1;
                for pair in 0..numbers.len().pow(2) {
                    assert!(
                        ((plus[pair] - minus[pair]) / (2.0 * step)
                            - response.c6_cartesian[pair * positions.len() + coordinate])
                            .abs()
                            < 1e-6
                    );
                }
            }
        }
    }

    #[test]
    fn published_widths_and_model_selection() {
        assert_eq!(widths().len(), ELEMENTS * ELEMENTS);
        for first in 0..ELEMENTS {
            for second in 0..ELEMENTS {
                let width = widths()[first * ELEMENTS + second];
                assert!(width > 0.0 && width <= 4.0);
                assert_eq!(width, widths()[second * ELEMENTS + first]);
            }
        }
        assert_eq!(widths()[27], 1.3221476510);
        let numbers = [1, 28, 1];
        let positions = [0.0, 0.0, 0.0, 3.0, 0.2, 0.0, 0.3, 1.4, 0.1];
        let damping = load_named("blyp", 1, false).unwrap();
        let evaluate = |model: Model, numbers: &[i32]| {
            model.dispersion(
                numbers,
                &positions,
                damping,
                None,
                RealspaceCutoff::default(),
                None,
                &[],
                WorkPartition::SERIAL,
            )
        };
        let d3 = evaluate(Model::D3, &numbers).unwrap().0.energy;
        assert!(
            (d3 - energy_partitioned(&numbers, &positions, damping, WorkPartition::SERIAL)
                .unwrap())
            .abs()
                < 1e-14
        );
        let d3s = evaluate(Model::D3S, &numbers).unwrap().0.energy;
        assert!((d3s - d3).abs() > 1e-6);
        assert!(evaluate(Model::D3S, &[1, 95, 1]).is_err());
        assert!(evaluate(Model::D3S, &[1, 0, 1]).is_err());
    }

    #[test]
    fn pair_width_interpolation_derivatives() {
        let first = 5;
        let second = 7;
        let cn_first = 1.3;
        let cn_second = 0.8;
        let ordinary = atomic_c6(
            first,
            second,
            &weights(first, cn_first),
            &weights(second, cn_second),
        );
        assert!(
            (coefficients(first, second, cn_first, cn_second, 4.0)[0] - ordinary).abs() < 1e-12
        );
        for width in [0.5, 1.7, 4.0] {
            let result = coefficients(first, second, cn_first, cn_second, width);
            let step = 1e-5;
            let plus = coefficients(first, second, cn_first + step, cn_second, width);
            let minus = coefficients(first, second, cn_first - step, cn_second, width);
            for (value, derivative) in [(0, 1), (1, 3), (2, 4)] {
                assert!(
                    ((plus[value] - minus[value]) / (2.0 * step) - result[derivative]).abs() < 1e-7
                );
            }
            let plus = coefficients(first, second, cn_first, cn_second + step, width);
            let minus = coefficients(first, second, cn_first, cn_second - step, width);
            for (value, derivative) in [(0, 2), (1, 4), (2, 5)] {
                assert!(
                    ((plus[value] - minus[value]) / (2.0 * step) - result[derivative]).abs() < 1e-7
                );
            }
        }
        assert!(coefficients(first, second, 100.0, 100.0, 4.0)
            .iter()
            .all(|value| value.is_finite()));
    }
}
