pub(crate) fn smooth_cutoff(distance2: f64, cutoff: f64, width: f64) -> (f64, f64) {
    let distance = distance2.sqrt();
    let width = width.min(cutoff);
    if width <= 0.0 || distance <= cutoff - width {
        (1.0, 0.0)
    } else if distance >= cutoff {
        (0.0, 0.0)
    } else {
        let fraction = (cutoff - distance) / width;
        (
            fraction.powi(3) * (10.0 + fraction * (-15.0 + 6.0 * fraction)),
            -15.0 * fraction.powi(2) * (1.0 - fraction).powi(2) / (width * distance),
        )
    }
}

pub(crate) fn displacement(positions: &[f64], first: usize, second: usize) -> [f64; 3] {
    std::array::from_fn(|axis| positions[3 * first + axis] - positions[3 * second + axis])
}

pub(crate) fn squared_distance(positions: &[f64], first: usize, second: usize) -> f64 {
    displacement(positions, first, second)
        .iter()
        .map(|value| value * value)
        .sum()
}

pub(crate) fn determinant(matrix: &[f64; 9]) -> f64 {
    matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7])
        - matrix[3] * (matrix[1] * matrix[8] - matrix[2] * matrix[7])
        + matrix[6] * (matrix[1] * matrix[5] - matrix[2] * matrix[4])
}

pub(crate) fn inverse(matrix: &[f64; 9]) -> [f64; 9] {
    let determinant = determinant(matrix);
    [
        matrix[4] * matrix[8] - matrix[5] * matrix[7],
        matrix[2] * matrix[7] - matrix[1] * matrix[8],
        matrix[1] * matrix[5] - matrix[2] * matrix[4],
        matrix[5] * matrix[6] - matrix[3] * matrix[8],
        matrix[0] * matrix[8] - matrix[2] * matrix[6],
        matrix[2] * matrix[3] - matrix[0] * matrix[5],
        matrix[3] * matrix[7] - matrix[4] * matrix[6],
        matrix[1] * matrix[6] - matrix[0] * matrix[7],
        matrix[0] * matrix[4] - matrix[1] * matrix[3],
    ]
    .map(|value| value / determinant)
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "a saturated repetition count is already an infeasible image sum"
)]
pub(crate) fn lattice_repetitions(
    lattice: &[f64; 9],
    periodic: [bool; 3],
    cutoff: f64,
) -> [i32; 3] {
    let reciprocal = periodic_reciprocal(lattice, periodic).unwrap_or([0.0; 9]);
    std::array::from_fn(|row| {
        if periodic[row] {
            (cutoff
                * (0..3)
                    .map(|axis| reciprocal[axis + 3 * row].powi(2))
                    .sum::<f64>()
                    .sqrt())
            .ceil() as i32
        } else {
            0
        }
    })
}

pub(crate) fn periodic_reciprocal(
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<[f64; 9], &'static str> {
    if lattice.iter().any(|value| !value.is_finite()) {
        return Err("lattice must be finite");
    }
    let active: Vec<_> = (0..3).filter(|&axis| periodic[axis]).collect();
    let mut result = [0.0; 9];
    let dot = |first: usize, second: usize| {
        (0..3)
            .map(|axis| lattice[axis + 3 * first] * lattice[axis + 3 * second])
            .sum::<f64>()
    };
    match active.as_slice() {
        [] => (),
        &[first] => {
            let norm = dot(first, first);
            if !norm.is_finite() || norm < 1.0e-20 {
                return Err("active lattice vectors are linearly dependent");
            }
            for axis in 0..3 {
                result[axis + 3 * first] = lattice[axis + 3 * first] / norm;
            }
        }
        &[first, second] => {
            let first_norm = dot(first, first);
            let second_norm = dot(second, second);
            let overlap = dot(first, second);
            let gram = first_norm * second_norm - overlap * overlap;
            if !gram.is_finite() || gram <= 1.0e-12 * first_norm * second_norm || gram < 1.0e-20 {
                return Err("active lattice vectors are linearly dependent");
            }
            for axis in 0..3 {
                result[axis + 3 * first] = (second_norm * lattice[axis + 3 * first]
                    - overlap * lattice[axis + 3 * second])
                    / gram;
                result[axis + 3 * second] = (first_norm * lattice[axis + 3 * second]
                    - overlap * lattice[axis + 3 * first])
                    / gram;
            }
        }
        _ => {
            let volume = determinant(lattice);
            if !volume.is_finite() || volume.abs() < 1.0e-12 {
                return Err("active lattice vectors are linearly dependent");
            }
            let inverse = inverse(lattice);
            result = std::array::from_fn(|index| inverse[index / 3 + 3 * (index % 3)]);
        }
    }
    if result.iter().any(|value| !value.is_finite()) {
        return Err("reciprocal lattice overflow");
    }
    Ok(result)
}

pub(crate) fn lattice_translation(lattice: &[f64; 9], image: [i32; 3]) -> [f64; 3] {
    std::array::from_fn(|axis| {
        (0..3)
            .map(|column| lattice[axis + 3 * column] * image[column] as f64)
            .sum()
    })
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "saturating only affects positions already astronomically far outside the cell"
)]
pub(crate) fn wrap_positions(
    positions: &[f64],
    lattice: &[f64; 9],
    periodic: [bool; 3],
) -> Result<Vec<f64>, &'static str> {
    let reciprocal = periodic_reciprocal(lattice, periodic)?;
    if !positions.len().is_multiple_of(3) || positions.iter().any(|value| !value.is_finite()) {
        return Err("positions must contain finite Cartesian triples");
    }
    let mut result = positions.to_vec();
    for (input, output) in positions
        .as_chunks::<3>()
        .0
        .iter()
        .zip(result.as_chunks_mut::<3>().0)
    {
        let image = std::array::from_fn(|column| {
            (0..3)
                .map(|axis| input[axis] * reciprocal[axis + 3 * column])
                .sum::<f64>()
                .floor() as i32
        });
        let shift = lattice_translation(lattice, image);
        for axis in 0..3 {
            output[axis] -= shift[axis];
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cutoff_switch_boundaries_and_squared_distance_derivative() {
        assert_eq!(smooth_cutoff(100.0, 5.0, 0.0), (1.0, 0.0));
        assert_eq!(smooth_cutoff(9.0, 5.0, 2.0), (1.0, 0.0));
        assert_eq!(smooth_cutoff(25.0, 5.0, 2.0), (0.0, 0.0));
        assert_eq!(smooth_cutoff(36.0, 5.0, 2.0), (0.0, 0.0));
        assert_eq!(smooth_cutoff(0.0, 5.0, 5.0), (1.0, 0.0));
        assert_eq!(smooth_cutoff(16.0, 5.0, 8.0), smooth_cutoff(16.0, 5.0, 5.0));
        for width in [2.0, 5.0, 8.0] {
            let actual = smooth_cutoff(16.0, 5.0, width);
            let step = 1.0e-5;
            let derivative = (smooth_cutoff(16.0 + step, 5.0, width).0
                - smooth_cutoff(16.0 - step, 5.0, width).0)
                / (2.0 * step);
            assert!((actual.1 - derivative).abs() < 1.0e-10);
        }
    }

    #[test]
    fn directional_image_ranges_use_active_reciprocal_vectors() {
        for periodic in [[true, false, false], [true, false, true], [true; 3]] {
            let mut lattice = [3.0, 1.0, 0.2, 2.5, 2.0, 0.1, 1.0, 0.5, 4.0];
            for column in 0..3 {
                if !periodic[column] {
                    lattice[3 * column..3 * column + 3].fill(0.0);
                }
            }
            let reciprocal = periodic_reciprocal(&lattice, periodic).unwrap();
            for first in 0..3 {
                for second in 0..3 {
                    let product: f64 = (0..3)
                        .map(|axis| reciprocal[axis + 3 * first] * lattice[axis + 3 * second])
                        .sum();
                    assert!(
                        (product
                            - if first == second && periodic[first] {
                                1.0
                            } else {
                                0.0
                            })
                        .abs()
                            < 1.0e-12
                    );
                }
            }
            let repetitions = lattice_repetitions(&lattice, periodic, 5.0);
            for column in 0..3 {
                assert_eq!(repetitions[column] == 0, !periodic[column]);
            }
        }
        assert!(periodic_reciprocal(&[0.0; 9], [true, false, false]).is_err());
    }
}
