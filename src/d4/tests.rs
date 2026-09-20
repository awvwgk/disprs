use super::*;

impl Local {
    fn variable(value: f64, index: usize) -> Self {
        let mut gradient = [0.0; 6];
        gradient[index] = 1.0;
        Self { value, gradient }
    }

    fn constant(value: f64) -> Self {
        Self {
            value,
            gradient: [0.0; 6],
        }
    }

    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
            gradient: std::array::from_fn(|index| self.gradient[index] + rhs.gradient[index]),
        }
    }

    fn sub(self, rhs: Self) -> Self {
        self.add(rhs.scale(-1.0))
    }

    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
            gradient: std::array::from_fn(|index| {
                self.gradient[index] * rhs.value + self.value * rhs.gradient[index]
            }),
        }
    }

    fn div(self, rhs: Self) -> Self {
        Self {
            value: self.value / rhs.value,
            gradient: std::array::from_fn(|index| {
                (self.gradient[index] * rhs.value - self.value * rhs.gradient[index])
                    / rhs.value.powi(2)
            }),
        }
    }

    fn powf(self, exponent: f64) -> Self {
        let value = self.value.powf(exponent);
        let derivative = exponent * self.value.powf(exponent - 1.0);
        Self {
            value,
            gradient: self.gradient.map(|item| item * derivative),
        }
    }

    fn sqrt(self) -> Self {
        self.powf(0.5)
    }
}

fn local_atm_dual(distances: [f64; 3], coefficients: [f64; 3], radius: f64, param: Param) -> Local {
    let distance: [Local; 3] =
        std::array::from_fn(|index| Local::variable(distances[index], index));
    let coefficient: [Local; 3] =
        std::array::from_fn(|index| Local::variable(coefficients[index], index + 3));
    let product2 = distance[0].mul(distance[1]).mul(distance[2]);
    let product = product2.sqrt();
    let c9 = coefficient[0]
        .mul(coefficient[1])
        .mul(coefficient[2])
        .sqrt()
        .scale(-param.s9);
    let damping = Local::constant(1.0).div(
        Local::constant(radius)
            .div(product)
            .powf(param.alpha / 3.0)
            .scale(6.0)
            .add(Local::constant(1.0)),
    );
    let angular = distance[0]
        .add(distance[2])
        .sub(distance[1])
        .mul(distance[0].sub(distance[2]).add(distance[1]))
        .mul(distance[0].scale(-1.0).add(distance[2]).add(distance[1]))
        .scale(0.375)
        .div(product2.powf(2.5))
        .add(Local::constant(1.0).div(product2.powf(1.5)));
    c9.scale(-1.0).mul(damping).mul(angular)
}

#[test]
fn second_order_tape_and_charge_solve() {
    use super::*;
    TAPE.with(|tape| unsafe { *tape.get() = Tape::default() });
    let coordinate = Dual::variable(0.0, 1, 0);
    let square = coordinate.clone() * coordinate;
    assert_eq!(square.hessian().unwrap(), [2.0]);
    TAPE.with(|tape| unsafe { *tape.get() = Tape::default() });
    let coordinate = Dual::variable(2.0, 1, 0);
    let mut matrix = [coordinate.clone()];
    let mut rhs = [Dual::constant(1.0, 1)];
    dual_solve(&mut matrix, &mut rhs).unwrap();
    let energy = rhs[0].clone() * rhs[0].clone();
    assert!((energy.hessian().unwrap()[0] - 0.375).abs() < 1e-14);
    let composed = (coordinate.clone().ln().exp() / coordinate.sqrt()).erf();
    let expected = -2.5 * (-2.0_f64).exp() / (2.0 * PI).sqrt() / 2.0;
    assert!((composed.hessian().unwrap()[0] - expected).abs() < 1e-14);
}

#[test]
fn directional_periodicity_preserves_inactive_axes_and_derivatives() {
    let numbers = [6, 8];
    let positions = [0.4, 0.8, 1.2, 3.1, 2.2, 1.7];
    for periodic in [[true, false, false], [true, false, true], [true; 3]] {
        for eeqbc in [false, true] {
            for d4s in [false, true] {
                let mut model = Model::custom(d4s, 2.0, 1.0, 6.0).unwrap();
                model.set_charge_model(i32::from(eeqbc)).unwrap();
                model.set_charge_cutoff(35.0).unwrap();
                model
                    .set_cutoff(Cutoff {
                        cn: 10.0,
                        disp2: 12.0,
                        disp3: 9.0,
                        width2: 2.0,
                        width3: 2.0,
                    })
                    .unwrap();
                let mut lattice = [7.3, 0.2, 0.0, 1.0, 8.3, 0.3, 0.4, 0.7, 9.3];
                let param = load_param("pbe", true).unwrap();
                let result = periodic_dispersion(
                    &numbers, &positions, 0.5, model, param, &lattice, periodic,
                )
                .unwrap();
                let (pair2, pair3) =
                    periodic_pairwise(&numbers, &positions, 0.5, model, param, &lattice, periodic)
                        .unwrap();
                assert!((pair2.iter().chain(&pair3).sum::<f64>() - result.energy).abs() < 1.0e-13);
                for column in 0..3 {
                    if !periodic[column] {
                        lattice[3 * column..3 * column + 3].fill(0.0);
                    }
                }
                let zero = periodic_dispersion(
                    &numbers, &positions, 0.5, model, param, &lattice, periodic,
                )
                .unwrap();
                assert!((result.energy - zero.energy).abs() < 1.0e-14);
                let properties =
                    periodic_properties(&numbers, &positions, 0.5, model, &lattice, periodic)
                        .unwrap();
                assert!((properties.charges.iter().sum::<f64>() - 0.5).abs() < 1.0e-13);
                for coordinate in 0..positions.len() {
                    let evaluate = |delta| {
                        let mut shifted = positions;
                        shifted[coordinate] += delta;
                        periodic_energy(&numbers, &shifted, 0.5, model, param, &lattice, periodic)
                            .unwrap()
                    };
                    assert!(
                        (result.gradient[coordinate]
                            - (evaluate(1.0e-5) - evaluate(-1.0e-5)) / 2.0e-5)
                            .abs()
                            < 1.0e-8
                    );
                }
                for component in 0..9 {
                    let evaluate = |delta| {
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
                        periodic_energy(&numbers, &shifted, 0.5, model, param, &cell, periodic)
                            .unwrap()
                    };
                    assert!(
                        (result.virial[component]
                            - (evaluate(1.0e-5) - evaluate(-1.0e-5)) / 2.0e-5)
                            .abs()
                            < 1.0e-8,
                        "eeqbc={eeqbc} component={component}"
                    );
                }
                let mut translated = positions;
                for axis in 0..3 {
                    translated[axis] += 3.0 * lattice[axis];
                }
                assert!(
                    (periodic_energy(&numbers, &translated, 0.5, model, param, &lattice, periodic)
                        .unwrap()
                        - result.energy)
                        .abs()
                        < 1.0e-13
                );
            }
        }
    }
}

#[test]
fn eeqbc_charge_response_conserves_charge() {
    let numbers = [8, 1, 1];
    let positions = [0.0, 0.0, 0.0, 0.0, 1.4, 1.0, 0.0, -1.4, 1.0];
    let origin = vec![Dual::constant(0.0, 0); 3];
    assert!(eeqbc_charges(&[112], &origin, 0.0, None).is_err());
    assert!(eeqbc_charges(&[8], &origin, f64::NAN, None).is_err());
    assert!(eeqbc_charges(&[8], &vec![Dual::constant(f64::NAN, 0); 3], 0.0, None).is_err());
    assert!(eeqbc_charges(&[8, 1], &vec![Dual::constant(0.0, 0); 6], 0.0, None).is_err());
    assert_eq!(
        eeqbc_charges(&[103], &origin, 1.0, None).unwrap()[0].value,
        1.0
    );
    for charge in [-1.0, 0.0, 1.0] {
        TAPE.with(|tape| unsafe {
            *tape.get() = Tape::default();
        });
        let coordinates: Vec<_> = positions
            .iter()
            .enumerate()
            .map(|(index, &value)| Dual::variable(value, positions.len(), index))
            .collect();
        let charges = eeqbc_charges(&numbers, &coordinates, charge, None).unwrap();
        assert!((charges.iter().map(|value| value.value).sum::<f64>() - charge).abs() < 1.0e-13);
        let response = charges[0].gradient();
        for coordinate in 0..positions.len() {
            let evaluate = |delta| {
                let mut shifted = positions;
                shifted[coordinate] += delta;
                let coordinates: Vec<_> = shifted
                    .iter()
                    .map(|&value| Dual::constant(value, 0))
                    .collect();
                eeqbc_charges(&numbers, &coordinates, charge, None).unwrap()[0].value
            };
            let expected = (evaluate(1.0e-5) - evaluate(-1.0e-5)) / 2.0e-5;
            assert!((response[coordinate] - expected).abs() < 1.0e-9);
        }
    }
}

#[test]
fn custom_model_and_cutoff_derivatives() {
    let numbers = [6, 8, 7];
    let positions = [0.0, 0.0, 0.0, 5.5, 0.0, 0.0, 1.0, 4.0, 0.0];
    let param = load_param("pbe", true).unwrap();
    for (d4s, eeqbc) in [(false, false), (true, false), (false, true), (true, true)] {
        for width in [0.0, 2.0, 8.0] {
            let mut model = Model::custom(d4s, 2.0, 1.0, if d4s { 6.0 } else { 5.0 }).unwrap();
            model.set_charge_model(i32::from(eeqbc)).unwrap();
            model
                .set_cutoff(Cutoff {
                    cn: 8.0,
                    disp2: 6.0,
                    disp3: 6.0,
                    width2: width,
                    width3: width,
                })
                .unwrap();
            let result = dispersion(&numbers, &positions, 0.0, model, param).unwrap();
            let scalar = energy(&numbers, &positions, 0.0, model, param).unwrap();
            let (pair2, pair3) = pairwise(&numbers, &positions, 0.0, model, param).unwrap();
            assert!((scalar - result.energy).abs() < 1.0e-14);
            assert!((pair2.iter().chain(&pair3).sum::<f64>() - scalar).abs() < 1.0e-14);
            for coordinate in 0..positions.len() {
                let mut shifted = positions;
                shifted[coordinate] += 1.0e-5;
                let plus = energy(&numbers, &shifted, 0.0, model, param).unwrap();
                shifted[coordinate] -= 2.0e-5;
                let minus = energy(&numbers, &shifted, 0.0, model, param).unwrap();
                assert!((result.gradient[coordinate] - (plus - minus) / 2.0e-5).abs() < 1.0e-9);
            }
        }
    }
}

#[test]
fn ghosts_and_partitions_preserve_environment_and_derivatives() {
    let numbers = [6, 8, 7, 1];
    let positions = [0.2, 0.3, 0.4, 2.6, 0.7, 0.8, 1.2, 2.8, 0.6, 0.7, 1.1, 2.3];
    let lattice = [5.3, 0.0, 0.0, 0.6, 5.7, 0.0, 0.3, 0.4, 6.2];
    let ghosts = [false, false, false, true];
    for mut base in [Model::D4, Model::D4S] {
        base.set_cutoff(Cutoff {
            cn: 6.0,
            disp2: 7.0,
            disp3: 6.0,
            width2: 1.0,
            width3: 1.0,
        })
        .unwrap();
        base.set_charge_cutoff(10.0).unwrap();
        for charge_model in [0, 1] {
            base.set_charge_model(charge_model).unwrap();
            for periodic in [
                [false; 3],
                [true, false, false],
                [true, true, false],
                [true; 3],
            ] {
                let selected = base.with_ghosts(&ghosts);
                let original =
                    periodic_properties(&numbers, &positions, 0.5, base, &lattice, periodic)
                        .unwrap();
                let environment =
                    periodic_properties(&numbers, &positions, 0.5, selected, &lattice, periodic)
                        .unwrap();
                assert_eq!(original.coordination, environment.coordination);
                assert_eq!(original.charges, environment.charges);
                assert_eq!(original.c6, environment.c6);
                assert_eq!(original.polarizabilities, environment.polarizabilities);
                for atm in [false, true] {
                    let param = load_param("pbe", atm).unwrap();
                    let serial = periodic_dispersion(
                        &numbers, &positions, 0.5, selected, param, &lattice, periodic,
                    )
                    .unwrap();
                    let pairs = periodic_pairwise(
                        &numbers, &positions, 0.5, selected, param, &lattice, periodic,
                    )
                    .unwrap();
                    assert!(
                        (serial.energy - pairs.0.iter().chain(&pairs.1).sum::<f64>()).abs()
                            < 1.0e-13
                    );
                    for matrix in [&pairs.0, &pairs.1] {
                        for atom in 0..4 {
                            assert_eq!(matrix[3 * 4 + atom], 0.0);
                            assert_eq!(matrix[atom * 4 + 3], 0.0);
                        }
                    }
                    assert!(serial.gradient[9..]
                        .iter()
                        .any(|value| value.abs() > 1.0e-12));
                    let mut summed_energy = 0.0;
                    let mut summed_gradient = [0.0; 12];
                    let mut summed_virial = [0.0; 9];
                    let mut summed_pairs = [0.0; 32];
                    for part in 0..3 {
                        let mut model = selected;
                        model.set_work_partition(part, 3).unwrap();
                        let result = periodic_dispersion(
                            &numbers, &positions, 0.5, model, param, &lattice, periodic,
                        )
                        .unwrap();
                        let value = periodic_energy(
                            &numbers, &positions, 0.5, model, param, &lattice, periodic,
                        )
                        .unwrap();
                        let pairs = periodic_pairwise(
                            &numbers, &positions, 0.5, model, param, &lattice, periodic,
                        )
                        .unwrap();
                        assert!((value - result.energy).abs() < 1.0e-13);
                        assert!(
                            (value - pairs.0.iter().chain(&pairs.1).sum::<f64>()).abs() < 1.0e-13
                        );
                        summed_energy += value;
                        for (total, value) in summed_gradient.iter_mut().zip(&result.gradient) {
                            *total += value;
                        }
                        for (total, value) in summed_virial.iter_mut().zip(result.virial) {
                            *total += value;
                        }
                        for (total, value) in
                            summed_pairs.iter_mut().zip(pairs.0.iter().chain(&pairs.1))
                        {
                            *total += value;
                        }
                        if part == 1 {
                            let step = 1.0e-5;
                            for coordinate in [0, 9] {
                                let mut plus = positions;
                                let mut minus = positions;
                                plus[coordinate] += step;
                                minus[coordinate] -= step;
                                let plus = periodic_energy(
                                    &numbers, &plus, 0.5, model, param, &lattice, periodic,
                                )
                                .unwrap();
                                let minus = periodic_energy(
                                    &numbers, &minus, 0.5, model, param, &lattice, periodic,
                                )
                                .unwrap();
                                assert!(
                                    (result.gradient[coordinate] - (plus - minus) / (2.0 * step))
                                        .abs()
                                        < 1.0e-8
                                );
                            }
                            let strained = |strain: f64| {
                                let mut shifted = positions;
                                let mut cell = lattice;
                                for coordinate in shifted.iter_mut().step_by(3) {
                                    *coordinate *= 1.0 + strain;
                                }
                                for coordinate in cell.iter_mut().step_by(3) {
                                    *coordinate *= 1.0 + strain;
                                }
                                periodic_energy(
                                    &numbers, &shifted, 0.5, model, param, &cell, periodic,
                                )
                                .unwrap()
                            };
                            let numerical = (strained(step) - strained(-step)) / (2.0 * step);
                            assert!((result.virial[0] - numerical).abs() < 1.0e-8,
                                "D4S={} charge_model={charge_model} periodic={periodic:?} ATM={atm}: {} != {numerical}",
                                base.d4s, result.virial[0]);
                        }
                    }
                    assert!((serial.energy - summed_energy).abs() < 1.0e-13);
                    for (expected, actual) in serial
                        .gradient
                        .iter()
                        .chain(&serial.virial)
                        .zip(summed_gradient.iter().chain(&summed_virial))
                    {
                        assert!((expected - actual).abs() < 1.0e-12);
                    }
                    for (expected, actual) in pairs.0.iter().chain(&pairs.1).zip(summed_pairs) {
                        assert!((expected - actual).abs() < 1.0e-13);
                    }
                    let mut empty = selected;
                    empty.set_work_partition(15, 16).unwrap();
                    for model in [empty, base.with_ghosts(&[true; 4])] {
                        let result = periodic_dispersion(
                            &numbers, &positions, 0.5, model, param, &lattice, periodic,
                        )
                        .unwrap();
                        assert_eq!(result.energy, 0.0);
                        assert!(result
                            .gradient
                            .iter()
                            .chain(&result.virial)
                            .all(|&value| value == 0.0));
                        let pairs = periodic_pairwise(
                            &numbers, &positions, 0.5, model, param, &lattice, periodic,
                        )
                        .unwrap();
                        assert!(pairs.0.iter().chain(&pairs.1).all(|&value| value == 0.0));
                    }
                }
            }
        }
    }
    let mut model = Model::D4;
    for (part, parts) in [(-1, 2), (2, 2), (0, 0)] {
        assert!(model.set_work_partition(part, parts).is_err());
        assert_eq!(model.partition, WorkPartition::SERIAL);
    }
    let invalid = model.with_ghosts(&[true]);
    let param = load_param("pbe", true).unwrap();
    for periodic in [[false; 3], [true; 3]] {
        assert!(
            periodic_energy(&numbers, &positions, 0.5, invalid, param, &lattice, periodic).is_err()
        );
        assert!(
            periodic_dispersion(&numbers, &positions, 0.5, invalid, param, &lattice, periodic)
                .is_err()
        );
        assert!(
            periodic_pairwise(&numbers, &positions, 0.5, invalid, param, &lattice, periodic)
                .is_err()
        );
        assert!(
            periodic_properties(&numbers, &positions, 0.5, invalid, &lattice, periodic).is_err()
        );
    }
}

#[test]
fn custom_reference_reproduces_default_data() {
    assert!(Model::custom(true, 3.0, 2.0, 5.0).is_err());
    for element in 0..ELEMENTS {
        for reference in 0..int(NREF, element) as usize {
            for frequency in 0..FREQUENCIES {
                let expected = alpha(element, reference, frequency, Model::D4);
                let actual = custom_alpha(element, reference, Model::D4)[frequency];
                assert!(
                    (actual - expected).abs() < 1.0e-12,
                    "{element} {reference} {frequency}: {actual} != {expected}"
                );
            }
        }
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(Model::custom(false, bad, 2.0, 6.0).is_err());
        assert!(Model::custom(false, 3.0, bad, 6.0).is_err());
        assert!(Model::custom(false, 3.0, 2.0, bad).is_err());
    }
}

#[test]
fn threaded_molecular_atm_matches_pairwise_and_differences() {
    let numbers: Vec<_> = (0..64).map(|index| [6, 8, 7, 1][index % 4]).collect();
    let ghosts: Vec<_> = (0..64).map(|index| index % 4 == 0).collect();
    let positions: Vec<_> = (0..64)
        .flat_map(|index| {
            [
                3.2 * (index % 4) as f64 + 0.11 * (index % 3) as f64,
                3.2 * ((index / 4) % 4) as f64 + 0.07 * (index % 5) as f64,
                3.2 * (index / 16) as f64 + 0.03 * (index % 7) as f64,
            ]
        })
        .collect();
    for model in [Model::D4, Model::D4S].into_iter().flat_map(|model| {
        let mut first = model.with_ghosts(&ghosts);
        first.set_work_partition(0, 2).unwrap();
        let mut second = first;
        second.set_work_partition(1, 2).unwrap();
        [model, first, second]
    }) {
        let param = load_param("pbe", true).unwrap();
        let actual = dispersion(&numbers, &positions, 0.5, model, param).unwrap();
        let (pair2, pair3) = pairwise(&numbers, &positions, 0.5, model, param).unwrap();
        assert!((actual.energy - pair2.iter().chain(&pair3).sum::<f64>()).abs() < 1.0e-12);
        for coordinate in [0, 73, 191] {
            let mut shifted = positions.clone();
            shifted[coordinate] += 1.0e-5;
            let plus = energy(&numbers, &shifted, 0.5, model, param).unwrap();
            shifted[coordinate] -= 2.0e-5;
            let minus = energy(&numbers, &shifted, 0.5, model, param).unwrap();
            assert!((actual.gradient[coordinate] - (plus - minus) / 2.0e-5).abs() < 1.0e-8);
        }
    }
}

#[test]
fn analytical_atm_matches_dual() {
    for distances in [
        [4.0, 5.0, 7.0],
        [1.0, 1.0, 2.0],
        [1.0, 1.0, 1.0 + (11.0_f64 / 3.0).sqrt()],
        [100.0, 300.0, 400.0],
    ] {
        for alpha in [12.0, 16.0, 18.0] {
            let param = Param {
                alpha,
                ..load_param("pbe", true).unwrap()
            };
            let actual = local_atm(distances, [10.0, 20.0, 30.0], 80.0, param);
            let expected = local_atm_dual(distances, [10.0, 20.0, 30.0], 80.0, param);
            for (actual, expected) in std::iter::once(actual.value)
                .chain(actual.gradient)
                .zip(std::iter::once(expected.value).chain(expected.gradient))
            {
                assert!(
                    (actual - expected).abs() < 1.0e-12 * expected.abs().max(1.0e-8),
                    "{actual} != {expected}"
                );
            }
        }
    }
}

#[test]
fn merged_parameters_preserve_aliases() {
    for (name, alias) in [
        ("blyp", "b-lyp"),
        ("r2scan", "r²scan"),
        ("pr2scan50", "pr²scan50"),
    ] {
        let expected = load_param(name, true).unwrap();
        let actual = load_param(alias, true).unwrap();
        assert!(actual == expected);
        let disabled =
            load_param(&format!("{}/def2-svp", alias.to_ascii_uppercase()), false).unwrap();
        assert_eq!(disabled.s9, 0.0);
        assert!(
            Param {
                s9: expected.s9,
                ..disabled
            } == expected
        );
    }
    assert_eq!(load_param("dsdblyp", true).unwrap().s6, 0.54);
    assert!(load_param("unknown", true).is_none());
}

#[test]
fn water_properties_match_dftd4() {
    let result = properties(
        &[8, 1, 1],
        &[0.0, 0.0, 0.0, 1.8, 0.0, 0.0, -0.45, 1.74, 0.0],
        0.0,
        Model::D4,
    )
    .unwrap();
    let expected_cn = [1.612225497817756, 0.8060518376720625, 0.8061736601456939];
    let expected_q = [-0.6039414350352019, 0.3017316160011496, 0.302209819034053];
    let expected_alpha = [6.774236258024333, 1.32860433606882, 1.3272508171202293];
    let expected_c6 = [
        25.092156987958635,
        4.160326733267335,
        4.156125142231804,
        4.160326733267335,
        0.715563313414557,
        0.7148398821073546,
        4.156125142231804,
        0.7148398821073546,
        0.7141171822166075,
    ];
    for (actual, expected) in result.coordination.iter().zip(expected_cn) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "CN {actual} != {expected}"
        );
    }
    for (actual, expected) in result.charges.iter().zip(expected_q) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "q {actual} != {expected}"
        );
    }
    for (actual, expected) in result.polarizabilities.iter().zip(expected_alpha) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "alpha {actual} != {expected}"
        );
    }
    for (actual, expected) in result.c6.iter().zip(expected_c6) {
        assert!(
            (actual - expected).abs() < 1.0e-11,
            "C6 {actual} != {expected}"
        );
    }
}

#[test]
fn water_dispersion_matches_dftd4() {
    let numbers = [8, 1, 1];
    let positions = [0.0, 0.0, 0.0, 1.8, 0.0, 0.0, -0.45, 1.74, 0.0];
    let param = load_param("pbe", true).unwrap();
    let result = dispersion(&numbers, &positions, 0.0, Model::D4, param).unwrap();
    assert!(
        (result.energy - -0.0001944600681561215).abs() < 1.0e-15,
        "energy {}",
        result.energy
    );
    let expected_gradient = [
        3.647740835467619e-5,
        4.7208592365321494e-5,
        0.0,
        -3.7211815935106094e-5,
        -8.873776610865173e-6,
        0.0,
        7.344075804299021e-7,
        -3.8334815754456324e-5,
        0.0,
    ];
    for (actual, expected) in result.gradient.iter().zip(expected_gradient) {
        assert!(
            (actual - expected).abs() < 2.0e-11,
            "gradient {actual} != {expected}"
        );
    }
    let expected_virial = [
        -6.73117520943844e-5,
        1.277869189948032e-6,
        0.0,
        1.277869189948032e-6,
        -6.670257941275398e-5,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    for (actual, expected) in result.virial.iter().zip(expected_virial) {
        assert!(
            (actual - expected).abs() < 3.0e-11,
            "virial {actual} != {expected}"
        );
    }
    let (pair2, pair3) = pairwise(&numbers, &positions, 0.0, Model::D4, param).unwrap();
    assert!(
        (pair2.iter().sum::<f64>() + pair3.iter().sum::<f64>() - result.energy).abs() < 1.0e-15
    );
}

#[test]
fn disabling_atm_preserves_two_body_parameters() {
    let enabled = load_param("pbe", true).unwrap();
    let disabled = load_param("pbe", false).unwrap();
    assert_eq!(disabled.s9, 0.0);
    for (method, expected) in [
        ("dftb(3ob)", [0.4727337, 0.5467502, 4.4955068]),
        ("dftb(mio)", [1.1948145, 0.6074567, 4.9336133]),
        ("dftb(ob2)", [2.7611320, 0.6037249, 5.3900004]),
        ("lc-dftb", [2.7611320, 0.6037249, 5.3900004]),
        ("dftb(matsci)", [2.7711819, 0.4681712, 5.2918629]),
        ("dftb(pbc)", [1.7303734, 0.5546548, 4.7973454]),
    ] {
        let param = load_param(method, false).unwrap();
        assert_eq!([param.s8, param.a1, param.a2], expected);
        assert_eq!(param.s9, 0.0);
        assert_ne!(param.s8, load_param(method, true).unwrap().s8);
    }
    assert_eq!(
        (
            disabled.s6,
            disabled.s8,
            disabled.a1,
            disabled.a2,
            disabled.alpha
        ),
        (
            enabled.s6,
            enabled.s8,
            enabled.a1,
            enabled.a2,
            enabled.alpha
        )
    );
}

#[test]
fn water_d4s_matches_dftd4() {
    let numbers = [8, 1, 1];
    let positions = [0.0, 0.0, 0.0, 1.8, 0.0, 0.0, -0.45, 1.74, 0.0];
    let properties = properties(&numbers, &positions, 0.0, Model::D4S).unwrap();
    let expected_c6 = [
        25.092156987958635,
        4.25161156746484,
        4.247256524860727,
        4.25161156746484,
        0.7756668537960016,
        0.774863383112792,
        4.247256524860727,
        0.774863383112792,
        0.774060744905505,
    ];
    for (actual, expected) in properties.c6.iter().zip(expected_c6) {
        assert!(
            (actual - expected).abs() < 1.0e-11,
            "C6 {actual} != {expected}"
        );
    }
    let result = dispersion(
        &numbers,
        &positions,
        0.0,
        Model::D4S,
        load_param("pbe", true).unwrap(),
    )
    .unwrap();
    assert!(
        (result.energy - -0.00019978999257333973).abs() < 1.0e-15,
        "energy {}",
        result.energy
    );
    let expected_gradient = [
        3.8436679444272216e-5,
        4.971354097523977e-5,
        0.0,
        -3.9359231111480783e-5,
        -9.229268239073332e-6,
        0.0,
        9.225516672085593e-7,
        -4.0484272736166436e-5,
        0.0,
    ];
    for (actual, expected) in result.gradient.iter().zip(expected_gradient) {
        assert!((actual - expected).abs() < 2.0e-11);
    }
    let expected_virial = [
        -7.126176425090927e-5,
        1.605239900942897e-6,
        0.0,
        1.6052399009428953e-6,
        -7.044263456092959e-5,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    for (actual, expected) in result.virial.iter().zip(expected_virial) {
        assert!((actual - expected).abs() < 3.0e-11);
    }
}

#[test]
fn analytical_gradient_matches_energy_differences() {
    let numbers = [6, 8, 1, 1];
    let positions = [
        0.1, -0.2, 0.3, 2.4, 0.4, -0.1, -0.7, 1.8, 0.5, 2.9, -1.2, 0.8,
    ];
    let param = load_param("pbe", true).unwrap();
    for model in [Model::D4, Model::D4S] {
        let result = dispersion(&numbers, &positions, 1.0, model, param).unwrap();
        let step = 1.0e-5;
        for coordinate in 0..positions.len() {
            let mut displaced = positions;
            displaced[coordinate] += step;
            let plus = energy(&numbers, &displaced, 1.0, model, param).unwrap();
            displaced[coordinate] -= 2.0 * step;
            let minus = energy(&numbers, &displaced, 1.0, model, param).unwrap();
            let numerical = (plus - minus) / (2.0 * step);
            assert!(
                (result.gradient[coordinate] - numerical).abs() < 2.0e-10,
                "gradient {} != {numerical}",
                result.gradient[coordinate]
            );
        }
    }
}

#[test]
fn periodic_gradient_matches_energy_differences() {
    let numbers = [6, 8, 1];
    let positions = [0.4, 0.8, 1.2, 2.9, 1.8, 1.6, 0.7, 2.6, 1.5];
    let lattice = [14.0, 0.2, 0.1, 1.0, 15.0, 0.3, 0.5, 0.7, 16.0];
    let step = 1.0e-5;
    let custom = [false, true].map(|d4s| {
        let mut model = Model::custom(d4s, 2.0, 1.0, 6.0).unwrap();
        model
            .set_cutoff(Cutoff {
                cn: 9.0,
                disp2: 16.0,
                disp3: 14.0,
                width2: 8.0,
                width3: 12.0,
            })
            .unwrap();
        model
    });
    let eeqbc = custom.map(|mut model| {
        model.set_charge_model(1).unwrap();
        model
    });
    for model in [
        Model::D4,
        Model::D4S,
        custom[0],
        custom[1],
        eeqbc[0],
        eeqbc[1],
    ] {
        for atm in [false, true] {
            let param = load_param("pbe", atm).unwrap();
            let result = periodic_differentiated_energy(
                &numbers, &positions, 0.5, model, param, &lattice, [true; 3],
            )
            .unwrap();
            let derivatives = result.gradient();
            for (index, actual) in derivatives.iter().enumerate() {
                let evaluate = |delta| {
                    let mut displaced = positions;
                    let mut strained = lattice;
                    if index < positions.len() {
                        displaced[index] += delta;
                    } else {
                        strained[index - positions.len()] += delta;
                    }
                    periodic_energy(
                        &numbers, &displaced, 0.5, model, param, &strained, [true; 3],
                    )
                    .unwrap()
                };
                let expected = (evaluate(step) - evaluate(-step)) / (2.0 * step);
                assert!(
                    (actual - expected).abs() < 2.0e-8,
                    "ATM={atm}, derivative {index}: {actual} != {expected}"
                );
            }
        }
    }
}

#[test]
#[allow(clippy::excessive_precision)]
fn periodic_eeqbc_matches_dftd4() {
    let numbers = [6, 8];
    let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
    let lattice = [8.0, 0.0, 0.0, 1.0, 9.0, 0.0, 0.5, 0.7, 10.0];
    let expected = [
        [
            -1.4513839361165394e-3,
            3.9099122196875746e-1,
            1.0900877803124251e-1,
            3.5856512090856342e1,
            2.1700820304864489e1,
            2.1700820304864489e1,
            1.4151693450490892e1,
            1.0524306369365418e1,
            4.9547770946788772,
            -6.8975180376508378e-6,
            1.4473815241690331e-5,
            -1.3378566454957086e-5,
            6.8975180376508353e-6,
            -1.4473815241690339e-5,
            1.3378566454957071e-5,
            2.0634963958677728e-3,
            -5.0604663350893989e-6,
            -1.1000792544108611e-7,
            -5.0604663350894192e-6,
            1.6609195776080693e-3,
            -8.2115561226596726e-5,
            -1.1000792544108844e-7,
            -8.2115561226596686e-5,
            1.3950722001261798e-3,
        ],
        [
            -1.4273910211090738e-3,
            3.9099122196875746e-1,
            1.0900877803124251e-1,
            3.4059976750485454e1,
            2.1542926489004987e1,
            2.1542926489004987e1,
            1.4151693450490892e1,
            1.0171326466955001e1,
            4.9547770946788772,
            -6.8508802980062081e-6,
            1.4304667769478148e-5,
            -1.3368093522766357e-5,
            6.8508802980062047e-6,
            -1.4304667769478152e-5,
            1.3368093522766422e-5,
            2.0137323221758067e-3,
            -5.8408578087453065e-6,
            -2.9277111447313049e-7,
            -5.8408578087452862e-6,
            1.6253572692413352e-3,
            -8.1750903910576876e-5,
            -2.9277111447313049e-7,
            -8.1750903910576876e-5,
            1.3697633754320409e-3,
        ],
    ];
    for (mut model, expected) in [Model::D4, Model::D4S].iter().copied().zip(expected.iter()) {
        model.set_charge_model(1).unwrap();
        let properties =
            periodic_properties(&numbers, &positions, 0.5, model, &lattice, [true; 3]).unwrap();
        let result = periodic_dispersion(
            &numbers,
            &positions,
            0.5,
            model,
            load_param("pbe", true).unwrap(),
            &lattice,
            [true; 3],
        )
        .unwrap();
        let actual: Vec<_> = [result.energy]
            .iter()
            .chain(&properties.charges)
            .chain(&properties.c6)
            .chain(&properties.polarizabilities)
            .chain(&result.gradient)
            .chain(&result.virial)
            .copied()
            .collect();
        for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (actual - expected).abs() < 1.0e-11,
                "component {index}: {actual} != {expected}"
            );
        }
    }
}

#[test]
fn periodic_skew_cell_matches_dftd4() {
    let numbers = [6, 8];
    let positions = [0.4, 0.8, 1.2, 4.4, 2.8, 3.6];
    let lattice = [8.0, 0.0, 0.0, 1.0, 9.0, 0.0, 0.5, 0.7, 10.0];
    let properties =
        periodic_properties(&numbers, &positions, 0.0, Model::D4, &lattice, [true; 3]).unwrap();
    for (actual, expected) in properties.coordination.iter().zip([9.28210730562645e-7; 2]) {
        assert!(
            (actual - expected).abs() < 1.0e-13,
            "CN {actual} != {expected}"
        );
    }
    for (actual, expected) in properties
        .charges
        .iter()
        .zip([0.22315075937800366, -0.22315075937800366])
    {
        assert!(
            (actual - expected).abs() < 1.0e-11,
            "q {actual} != {expected}"
        );
    }
    let result = periodic_dispersion(
        &numbers,
        &positions,
        0.0,
        Model::D4,
        load_param("pbe", true).unwrap(),
        &lattice,
        [true; 3],
    )
    .unwrap();
    assert!((result.energy - -0.0017844552381591887).abs() < 1.0e-11);
    let expected_gradient = [
        -8.521759140145904e-6,
        1.8232600190767648e-5,
        -1.6009068929161304e-5,
        8.521759140145902e-6,
        -1.8232600190767665e-5,
        1.600906892916125e-5,
    ];
    for (actual, expected) in result.gradient.iter().zip(expected_gradient) {
        assert!(
            (actual - expected).abs() < 1.0e-10,
            "gradient {actual} != {expected}"
        );
    }
    let expected_virial = [
        0.0025063560142267196,
        -8.07404988915801e-6,
        -6.738213250197521e-7,
        -8.074049889157988e-6,
        0.0020290128351334548,
        -0.00010268081679581214,
        -6.738213250197706e-7,
        -0.00010268081679581214,
        0.0017130518874741672,
    ];
    for (actual, expected) in result.virial.iter().zip(expected_virial) {
        assert!(
            (actual - expected).abs() < 1.0e-9,
            "virial {actual} != {expected}"
        );
    }
    let d4s = periodic_dispersion(
        &numbers,
        &positions,
        0.0,
        Model::D4S,
        load_param("pbe", true).unwrap(),
        &lattice,
        [true; 3],
    )
    .unwrap();
    assert!((d4s.energy - -0.0017580448717701273).abs() < 1.0e-11);
    let expected_d4s_gradient = [
        -8.46625649705616e-6,
        1.7981993098558662e-5,
        -1.6060825994597856e-5,
        8.466256497056154e-6,
        -1.798199309855868e-5,
        1.606082599459792e-5,
    ];
    for (actual, expected) in d4s.gradient.iter().zip(expected_d4s_gradient) {
        assert!((actual - expected).abs() < 1.0e-10);
    }
    let expected_d4s_virial = [
        0.0024540668369781695,
        -8.856063491659248e-6,
        -8.533819862013643e-7,
        -8.856063491659256e-6,
        0.001991313446331876,
        -0.00010205532772350778,
        -8.533819862013641e-7,
        -0.00010205532772350781,
        0.001686094254830103,
    ];
    for (actual, expected) in d4s.virial.iter().zip(expected_d4s_virial) {
        assert!((actual - expected).abs() < 1.0e-9);
    }
    let partial = periodic_energy(
        &numbers,
        &positions,
        0.0,
        Model::D4,
        load_param("pbe", true).unwrap(),
        &lattice,
        [true, false, false],
    )
    .unwrap();
    assert!((partial - result.energy).abs() > 1.0e-8);
}
