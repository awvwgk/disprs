use super::{Dual, PI};

unsafe extern "C" {
    fn erfc(value: f64) -> f64;
}

fn dot(left: &[Dual; 3], right: &[Dual; 3]) -> Dual {
    (0..3).fold(Dual::constant(0.0, left[0].size), |sum, axis| {
        sum + left[axis].clone() * right[axis].clone()
    })
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "callers exclude the zero mode, so wave is positive and the count is finite"
)]
fn incomplete_bessel(wave: f64, transverse: f64) -> [f64; 6] {
    let nodes = [
        0.09501250983763744,
        0.2816035507792589,
        0.4580167776572274,
        0.6178762444026438,
        0.755404408355003,
        0.8656312023878318,
        0.9445750230732326,
        0.9894009349916499,
    ];
    let weights = [
        0.1894506104550685,
        0.1826034150449236,
        0.16915651939500254,
        0.14959598881657673,
        0.12462897125553388,
        0.09515851168249278,
        0.06225352393864789,
        0.027152459411754096,
    ];
    let upper = (50.0 / wave).ln().max(0.0);
    let intervals = (upper / 0.5).ceil().max(1.0) as usize;
    let half = upper / (2.0 * intervals as f64);
    let mut result = [0.0; 6];
    for interval in 0..intervals {
        let center = (2 * interval + 1) as f64 * half;
        for (node, weight) in nodes.iter().zip(weights) {
            for sign in [-1.0, 1.0] {
                let scale = (center + sign * half * node).exp();
                let value = half * weight * (-wave * scale - transverse / scale).exp();
                result[0] += value;
                result[1] -= value * scale;
                result[2] -= value / scale;
                result[3] += value * scale * scale;
                result[4] += value;
                result[5] += value / (scale * scale);
            }
        }
    }
    result
}

fn wire_mode(wave: Dual, transverse: Dual) -> Dual {
    let [value, wave_partial, transverse_partial, wave_second, cross, transverse_second] =
        incomplete_bessel(wave.value, transverse.value);
    Dual::record_second(
        value,
        wave.size,
        wave.node.map(|node| (node, wave_partial)),
        transverse.node.map(|node| (node, transverse_partial)),
        [[wave_second, cross], [cross, transverse_second]],
    )
}

fn entire_exponential_integral(argument: Dual) -> Dual {
    let input = argument.value;
    let value = if input <= 1.0 {
        let mut term = input;
        let mut sum = term;
        for order in 2..=32 {
            term *= -input / order as f64;
            sum += term / order as f64;
        }
        sum
    } else {
        0.577_215_664_901_532_9 + input.ln() + incomplete_bessel(input, 0.0)[0]
    };
    argument.unary_second(
        value,
        if input == 0.0 {
            1.0
        } else {
            -(-input).exp_m1() / input
        },
        if input.abs() < 1e-4 {
            -0.5 + input * (1.0 / 3.0 + input * (-0.125 + input * (1.0 / 30.0 - input / 144.0)))
        } else {
            ((input + 1.0) * (-input).exp() - 1.0) / input.powi(2)
        },
    )
}

fn slab_mode(wave: Dual, height: Dual, alpha: f64) -> Dual {
    let wave_value = wave.value;
    let height_value = height.value;
    let mut value = 0.0;
    let mut wave_partial = 0.0;
    let mut height_partial = 0.0;
    let mut wave_second = 0.0;
    let mut cross = 0.0;
    let mut height_second = 0.0;
    let gaussian = (-(wave_value / (2.0 * alpha)).powi(2) - (alpha * height_value).powi(2)).exp();
    for sign in [-1.0, 1.0] {
        let argument = wave_value / (2.0 * alpha) + sign * alpha * height_value;
        let term = if argument > 26.0 {
            let mut series = 1.0;
            let mut term = 1.0;
            for order in 1..=12 {
                term *= -((2 * order - 1) as f64) / (2.0 * argument * argument);
                series += term;
            }
            gaussian * series / (PI.sqrt() * argument)
        } else {
            (sign * wave_value * height_value).exp() * unsafe { erfc(argument) }
        };
        value += term / wave_value;
        wave_partial += (sign * height_value * term - gaussian / (alpha * PI.sqrt())) / wave_value
            - term / wave_value.powi(2);
        height_partial += sign * term;
        let term_wave = sign * height_value * term - gaussian / (alpha * PI.sqrt());
        wave_second += (height_value.powi(2) * term
            - sign * height_value * gaussian / (alpha * PI.sqrt())
            + wave_value * gaussian / (2.0 * alpha.powi(3) * PI.sqrt()))
            / wave_value
            - 2.0 * term_wave / wave_value.powi(2)
            + 2.0 * term / wave_value.powi(3);
        cross += sign * term_wave;
        height_second += wave_value * term - 2.0 * alpha * gaussian / PI.sqrt();
    }
    Dual::record_second(
        value,
        wave.size,
        wave.node.map(|node| (node, wave_partial)),
        height.node.map(|node| (node, height_partial)),
        [[wave_second, cross], [cross, height_second]],
    )
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "periodic_reciprocal rejects degenerate cells, so the image bounds stay finite"
)]
pub(super) fn gaussian_pair(
    displacement: &[Dual; 3],
    lattice: &[Dual; 9],
    periodic: [bool; 3],
    gamma: Dual,
    alpha: f64,
) -> Dual {
    let size = gamma.size;
    let constant = |value| Dual::constant(value, size);
    let active: Vec<_> = (0..3).filter(|&column| periodic[column]).collect();
    assert!((1..=2).contains(&active.len()));
    let vectors: Vec<[Dual; 3]> = active
        .iter()
        .map(|&column| std::array::from_fn(|axis| lattice[axis + 3 * column].clone()))
        .collect();
    let norm = dot(&vectors[0], &vectors[0]);
    let (measure, reciprocal, transverse): (Dual, Vec<[Dual; 3]>, Dual) = if active.len() == 1 {
        let projection = dot(displacement, &vectors[0]) / norm.clone();
        let normal = std::array::from_fn(|axis| {
            displacement[axis].clone() - vectors[0][axis].clone() * projection.clone()
        });
        (
            norm.clone().sqrt(),
            vec![std::array::from_fn(|axis| {
                vectors[0][axis].clone() / norm.clone()
            })],
            dot(&normal, &normal),
        )
    } else {
        let normal: [Dual; 3] = std::array::from_fn(|axis| {
            vectors[0][(axis + 1) % 3].clone() * vectors[1][(axis + 2) % 3].clone()
                - vectors[0][(axis + 2) % 3].clone() * vectors[1][(axis + 1) % 3].clone()
        });
        let gram = dot(&normal, &normal);
        let area = gram.clone().sqrt();
        let other_norm = dot(&vectors[1], &vectors[1]);
        let overlap = dot(&vectors[0], &vectors[1]);
        let reciprocal = vec![
            std::array::from_fn(|axis| {
                (vectors[0][axis].clone() * other_norm.clone()
                    - vectors[1][axis].clone() * overlap.clone())
                    / gram.clone()
            }),
            std::array::from_fn(|axis| {
                (vectors[1][axis].clone() * norm.clone()
                    - vectors[0][axis].clone() * overlap.clone())
                    / gram.clone()
            }),
        ];
        let height = dot(displacement, &normal) / area.clone();
        (area, reciprocal, height)
    };
    let center: Vec<_> = reciprocal
        .iter()
        .map(|vector| dot(displacement, vector).value.round())
        .collect();
    let base: [Dual; 3] = std::array::from_fn(|axis| {
        vectors
            .iter()
            .zip(&center)
            .fold(displacement[axis].clone(), |value, (vector, image)| {
                value - vector[axis].clone() * *image
            })
    });
    let cutoff = 7.0 / alpha.min(gamma.value);
    let bounds: Vec<_> = reciprocal
        .iter()
        .map(|vector| (cutoff * dot(vector, vector).value.sqrt() + 0.5).ceil() as i32)
        .collect();
    let second_bound = if active.len() == 2 { bounds[1] } else { 0 };
    let mut result = constant(0.0);
    for first in -bounds[0]..=bounds[0] {
        for second in -second_bound..=second_bound {
            let image = [first, second];
            let vector = std::array::from_fn(|axis| {
                vectors
                    .iter()
                    .enumerate()
                    .fold(base[axis].clone(), |value, (column, vector)| {
                        value - vector[axis].clone() * image[column] as f64
                    })
            });
            let distance2 = dot(&vector, &vector);
            if distance2.value > cutoff * cutoff {
                continue;
            }
            if distance2.value < 1.0e-24 {
                result += (gamma.clone() - constant(alpha)) * (2.0 / PI.sqrt());
                result = result
                    - distance2
                        * (gamma.clone().powf(3.0) - constant(alpha.powi(3)))
                        * (2.0 / (3.0 * PI.sqrt()));
            } else {
                let distance = distance2.sqrt();
                result += ((gamma.clone() * distance.clone()).erf()
                    - (distance.clone() * alpha).erf())
                    / distance;
            }
        }
    }
    let wave_cutoff = 14.0 * alpha;
    let wave_bounds: Vec<_> = vectors
        .iter()
        .map(|vector| (wave_cutoff * dot(vector, vector).value.sqrt() / (2.0 * PI)).ceil() as i32)
        .collect();
    let second_bound = if active.len() == 2 { wave_bounds[1] } else { 0 };
    for first in -wave_bounds[0]..=wave_bounds[0] {
        for second in -second_bound..=second_bound {
            if first == 0 && second == 0 {
                continue;
            }
            let image = [first, second];
            let wave = std::array::from_fn(|axis| {
                reciprocal
                    .iter()
                    .enumerate()
                    .fold(constant(0.0), |value, (column, vector)| {
                        value + vector[axis].clone() * (2.0 * PI * image[column] as f64)
                    })
            });
            let wave2 = dot(&wave, &wave);
            if wave2.value > wave_cutoff * wave_cutoff {
                continue;
            }
            let phase = dot(&wave, &base);
            let cosine = phase.clone().unary_second(
                phase.value.cos(),
                -phase.value.sin(),
                -phase.value.cos(),
            );
            let mode = if active.len() == 1 {
                wire_mode(
                    wave2 * (0.25 / alpha.powi(2)),
                    transverse.clone() * alpha.powi(2),
                )
            } else {
                slab_mode(wave2.sqrt(), transverse.clone(), alpha) * PI
            };
            result += cosine * mode / measure.clone();
        }
    }
    if active.len() == 1 {
        result += (constant(2.0 * alpha.ln() + 0.577_215_664_901_532_9)
            - entire_exponential_integral(transverse * alpha.powi(2)))
            / measure;
    } else {
        result = result
            - ((transverse.clone() * alpha).erf() * transverse.clone()
                + (transverse.clone() * transverse * (-alpha * alpha)).exp()
                    * (1.0 / (alpha * PI.sqrt())))
                * (2.0 * PI)
                / measure;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::super::erf;
    use super::*;

    #[test]
    fn gaussian_ewald_split_and_neutral_image_sum() {
        let lattice = [6.1, 0.3, 0.2, 1.1, 7.2, 0.5, 0.0, 0.0, 0.0];
        let displacement = [1.2, 0.8, 1.7];
        for periodic in [[true, false, false], [true, true, false]] {
            let evaluate = |position: [f64; 3], alpha| {
                gaussian_pair(
                    &position.map(|value| Dual::constant(value, 0)),
                    &lattice.map(|value| Dual::constant(value, 0)),
                    periodic,
                    Dual::constant(0.7, 0),
                    alpha,
                )
                .value
            };
            let expected = evaluate(displacement, 0.4);
            for alpha in [0.25, 0.6, 0.9] {
                assert!((evaluate(displacement, alpha) - expected).abs() < 2.0e-12);
            }
            let neutral = evaluate([0.0; 3], 0.4) - expected;
            let direct = |extent: i32| {
                let mut sum = 0.0;
                let second_extent = if periodic[1] { extent } else { 0 };
                for first in -extent..=extent {
                    for second in -second_extent..=second_extent {
                        let shift: [f64; 3] = std::array::from_fn(|axis| {
                            first as f64 * lattice[axis] + second as f64 * lattice[axis + 3]
                        });
                        let potential = |offset: [f64; 3]| {
                            let distance = (0..3)
                                .map(|axis| (shift[axis] + offset[axis]).powi(2))
                                .sum::<f64>()
                                .sqrt();
                            if distance == 0.0 {
                                1.4 / PI.sqrt()
                            } else {
                                (unsafe { erf(0.7 * distance) }) / distance
                            }
                        };
                        sum += potential([0.0; 3])
                            - 0.5
                                * (potential(displacement)
                                    + potential(displacement.map(|value| -value)));
                    }
                }
                sum
            };
            let small = direct(80);
            let large = direct(160);
            assert!((large - neutral).abs() < (small - neutral).abs());
            let exponent = if periodic[1] { 1 } else { 2 };
            let small_radius = 80.5_f64.powi(exponent);
            let large_radius = 160.5_f64.powi(exponent);
            let extrapolated =
                (large_radius * large - small_radius * small) / (large_radius - small_radius);
            assert!(
                (extrapolated - neutral).abs() < 1.0e-8,
                "direct={extrapolated}, Ewald={neutral}"
            );
        }
    }

    #[test]
    fn gaussian_ewald_special_functions_and_responses() {
        for (wave, transverse, expected) in [
            (1.0, 0.0, 0.219_383_934_395_520_3),
            (1.0, 1.0, 0.113_893_872_749_533_4),
            (0.5, 0.5, 0.421_024_438_240_708_3),
        ] {
            assert!((incomplete_bessel(wave, transverse)[0] - expected).abs() < 2.0e-15);
        }
        for periodic in [[true, false, false], [true, false, true]] {
            for displacement in [[0.0; 3], [1.2, 0.8, 1.7], [0.2, 150.0, -180.0]] {
                let lattice = [6.1, 0.3, 0.2, 0.0, 0.0, 0.0, 1.1, 7.2, 0.5];
                let mut inputs = [0.0; 13];
                inputs[..3].copy_from_slice(&displacement);
                inputs[3..12].copy_from_slice(&lattice);
                inputs[12] = 0.7;
                let evaluate = |values: [f64; 13], derivatives, alpha| {
                    super::super::TAPE
                        .with(|tape| unsafe { *tape.get() = super::super::Tape::default() });
                    let variables: [Dual; 13] = std::array::from_fn(|index| {
                        Dual::variable(values[index], derivatives, index)
                    });
                    gaussian_pair(
                        &std::array::from_fn(|axis| variables[axis].clone()),
                        &std::array::from_fn(|index| variables[3 + index].clone()),
                        periodic,
                        variables[12].clone(),
                        alpha,
                    )
                };
                let result = evaluate(inputs, 13, 0.4);
                let gradient = result.gradient();
                let hessian = result.hessian().unwrap();
                assert!(result.value.is_finite());
                assert!((result.value - evaluate(inputs, 0, 0.7).value).abs() < 2.0e-12);
                for index in 0..13 {
                    let mut plus = inputs;
                    let mut minus = inputs;
                    plus[index] += 1.0e-5;
                    minus[index] -= 1.0e-5;
                    let numerical =
                        (evaluate(plus, 0, 0.4).value - evaluate(minus, 0, 0.4).value) / 2.0e-5;
                    assert!(
                        (gradient[index] - numerical).abs() < 2.0e-8,
                        "{periodic:?} {displacement:?} derivative {index}: {} != {numerical}",
                        gradient[index]
                    );
                    plus[index] = inputs[index] + 1e-4;
                    minus[index] = inputs[index] - 1e-4;
                    let plus_gradient = evaluate(plus, 13, 0.4).gradient();
                    let minus_gradient = evaluate(minus, 13, 0.4).gradient();
                    for row in 0..13 {
                        let expected = (plus_gradient[row] - minus_gradient[row]) / 2e-4;
                        assert!(
                            (hessian[index * 13 + row] - expected).abs() < 2e-7,
                            "{periodic:?} {displacement:?} Hessian {row},{index}: {} != {expected}",
                            hessian[index * 13 + row]
                        );
                    }
                }
                for column in 0..3 {
                    if periodic[column] {
                        let mut shifted = inputs;
                        for axis in 0..3 {
                            shifted[axis] += 11.0 * lattice[axis + 3 * column];
                        }
                        assert!((evaluate(shifted, 0, 0.4).value - result.value).abs() < 2.0e-12);
                    }
                }
            }
        }
    }
}
