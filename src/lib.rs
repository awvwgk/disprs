//! Native D3/D3S, gCP, and D4/D4S dispersion corrections.
//!
//! Use the [D3 module](crate::d3) for D3-family models and gCP, or the
//! [D4 module](crate::d4) for charge-dependent D4.
//! Coordinates and lattice vectors are in bohr; energies are in hartree.
//! Cartesian gradients are energy derivatives, not forces.
//! Geometry arrays store consecutive `(x, y, z)` triples for each atom.
//! See each module for supported elements, parameter choices, and periodic limits.

use std::ffi::c_char;

pub mod d3;
mod d3_ext;
pub mod d4;
mod d4_ext;
pub(crate) mod dual;
mod ffi;
mod geometry;
mod parallel;
mod parameters;

/// Cyclic ownership of interaction terms for externally partitioned calculations.
///
/// Use [`WorkPartition::SERIAL`] for the full calculation, or [`WorkPartition::new`] for worker
/// `part` of `parts`. Sum partial energies, gradients, Hessians, virials, and
/// pair matrices over every worker to recover the full result. All workers must
/// use identical inputs and settings except their worker index.
/// This does not create threads or perform communication.
///
/// Also available as [`d3::WorkPartition`] and [`d4::WorkPartition`].
/// The worker count is always positive and the index is less than the count.
///
/// ```
/// use disprs::WorkPartition;
/// let worker = WorkPartition::new(1, 4).unwrap();
/// assert_eq!(worker.part(), 1);
/// assert_eq!(worker.parts(), 4);
/// assert_eq!(WorkPartition::new(0, 1), Some(WorkPartition::SERIAL));
/// for (part, parts) in [(-1, 4), (0, 0), (0, -1), (4, 4)] {
///     assert!(WorkPartition::new(part, parts).is_none());
/// }
/// let _: disprs::d3::WorkPartition = worker;
/// let _: disprs::d4::WorkPartition = worker;
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WorkPartition {
    /// Invalid states cannot be constructed directly:
    ///
    /// ```compile_fail
    /// let invalid = disprs::WorkPartition { part: 0, parts: 0 };
    /// ```
    part: usize,
    parts: usize,
}

impl WorkPartition {
    /// Own every interaction term: worker zero of one.
    pub const SERIAL: Self = Self { part: 0, parts: 1 };

    /// Create worker `part` of `parts`; return `None` unless `0 <= part < parts`.
    pub fn new(part: i32, parts: i32) -> Option<Self> {
        (parts > 0 && part >= 0 && part < parts).then_some(Self {
            part: part as usize,
            parts: parts as usize,
        })
    }

    /// Zero-based worker index, strictly less than [`WorkPartition::parts`].
    pub const fn part(self) -> usize {
        self.part
    }

    /// Positive number of workers whose partial results must be summed.
    pub const fn parts(self) -> usize {
        self.parts
    }

    pub(crate) fn owns_pair(self, first: usize, second: usize) -> bool {
        self.owns_index(first * (first + 1) / 2 + second)
    }

    pub(crate) fn owns_index(self, index: usize) -> bool {
        index % self.parts == self.part
    }
}

static VERSION: &[u8] = b"0.1.0\0";

#[no_mangle]
/// Return the library version as a static, NUL-terminated C string.
///
/// The pointer remains valid for the process lifetime; do not modify or free it.
///
/// ```
/// let version = unsafe { std::ffi::CStr::from_ptr(disprs::disprs_get_version()) };
/// assert_eq!(version.to_str().unwrap(), "0.1.0");
/// ```
pub extern "C" fn disprs_get_version() -> *const c_char {
    VERSION.as_ptr().cast()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reports_package_version() {
        assert_eq!(
            &VERSION[..VERSION.len() - 1],
            env!("CARGO_PKG_VERSION").as_bytes()
        );
    }
}
