use std::ffi::{c_char, c_int, c_void};
use std::ptr::null_mut;

pub(crate) type Handle = *mut c_void;

pub(crate) struct Structure {
    pub numbers: Vec<i32>,
    pub positions: Vec<f64>,
    pub charge: f64,
    pub lattice: Option<[f64; 9]>,
    pub periodic: [bool; 3],
}

impl Structure {
    #[allow(clippy::type_complexity)]
    pub fn cell(&self) -> Result<Option<(&[f64; 9], [bool; 3])>, &'static str> {
        if self.periodic.iter().any(|&active| active) {
            Ok(Some((
                self.lattice
                    .as_ref()
                    .ok_or("periodic structure requires a lattice")?,
                self.periodic,
            )))
        } else {
            Ok(None)
        }
    }
}

pub(crate) unsafe fn as_structure<'a>(handle: Handle) -> &'a mut Structure {
    unsafe { &mut *handle.cast() }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_new_structure(
    error: Handle,
    natoms: c_int,
    numbers: *const c_int,
    positions: *const f64,
    charge: *const f64,
    lattice: *const f64,
    periodic: *const bool,
) -> Handle {
    unsafe {
        guard(error, || {
            if numbers.is_null() || positions.is_null() {
                set_error(error, "atomic numbers and positions are required");
                return null_mut();
            }
            let count = match atom_count(natoms) {
                Ok(count) => count,
                Err(message) => {
                    set_error(error, message);
                    return null_mut();
                }
            };
            let structure = Structure {
                numbers: std::slice::from_raw_parts(numbers, count).to_vec(),
                positions: std::slice::from_raw_parts(positions, 3 * count).to_vec(),
                charge: charge.as_ref().copied().unwrap_or(0.0),
                lattice: (!lattice.is_null())
                    .then(|| std::slice::from_raw_parts(lattice, 9).try_into().unwrap()),
                periodic: if periodic.is_null() {
                    [!lattice.is_null(); 3]
                } else {
                    std::slice::from_raw_parts(periodic, 3).try_into().unwrap()
                },
            };
            if !structure.charge.is_finite() {
                set_error(error, "charge must be finite");
                return null_mut();
            }
            if let Err(message) = validate_structure(
                &structure.numbers,
                &structure.positions,
                structure.lattice.as_ref(),
                structure.periodic,
                118,
            ) {
                set_error(error, message);
                return null_mut();
            }
            Box::into_raw(Box::new(structure)).cast()
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_update_structure(
    error: Handle,
    structure: Handle,
    positions: *const f64,
    lattice: *const f64,
) {
    unsafe {
        guard(error, || {
            if !require_handles(error, &[structure]) || positions.is_null() {
                set_error(error, "structure and positions are required");
                return;
            }
            let structure = as_structure(structure);
            let positions = std::slice::from_raw_parts(positions, structure.positions.len());
            let cell = if lattice.is_null() {
                structure.lattice
            } else {
                Some(std::slice::from_raw_parts(lattice, 9).try_into().unwrap())
            };
            if let Err(message) = validate_structure(
                &structure.numbers,
                positions,
                cell.as_ref(),
                structure.periodic,
                118,
            ) {
                set_error(error, message);
                return;
            }
            structure.positions.copy_from_slice(positions);
            structure.lattice = cell;
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_delete_structure(handle: *mut Handle) {
    unsafe {
        if let Some(handle) = handle.as_mut() {
            if !handle.is_null() {
                drop(Box::from_raw(handle.cast::<Structure>()));
                *handle = null_mut();
            }
        }
    }
}

pub(crate) unsafe fn guard<T: Default>(error: Handle, operation: impl FnOnce() -> T) -> T {
    unsafe {
        evaluate(error, || {
            #[cfg(test)]
            PANIC_ON_SETUP.with(|pending| {
                assert!(!pending.replace(false), "injected setup failure");
            });
            Ok(operation())
        })
        .unwrap_or_default()
    }
}

#[cfg(test)]
std::thread_local! {
    pub(crate) static PANIC_ON_SETUP: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

#[cfg(test)]
pub(crate) unsafe fn assert_setup_panic<T: Default + PartialEq>(
    error: Handle,
    operation: impl FnOnce() -> T,
) {
    PANIC_ON_SETUP.with(|pending| pending.set(true));
    assert!(operation() == T::default());
    assert!(!PANIC_ON_SETUP.with(|pending| pending.get()));
    assert_eq!(
        unsafe { as_error(error) }.message.as_deref(),
        Some("native calculation panicked")
    );
}

pub(crate) unsafe fn evaluate<T>(
    error: Handle,
    operation: impl FnOnce() -> Result<T, &'static str>,
) -> Option<T> {
    unsafe {
        clear_error(error);
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(operation)) {
            Ok(Ok(value)) => Some(value),
            Ok(Err(message)) => {
                set_error(error, message);
                None
            }
            Err(_) => {
                set_error(error, "native calculation panicked");
                None
            }
        }
    }
}

pub(crate) fn atom_count(natoms: c_int) -> Result<usize, &'static str> {
    let count = usize::try_from(natoms).ok().filter(|&count| count > 0);
    count
        .filter(|&count| {
            count
                .checked_mul(3)
                .and_then(|coordinates| coordinates.checked_pow(2))
                .and_then(|entries| entries.checked_mul(std::mem::size_of::<f64>()))
                .is_some_and(|bytes| bytes <= isize::MAX as usize)
        })
        .ok_or("invalid or unrepresentable atom count")
}

pub(crate) fn validate_structure(
    numbers: &[i32],
    positions: &[f64],
    lattice: Option<&[f64; 9]>,
    periodic: [bool; 3],
    max_element: i32,
) -> Result<(), &'static str> {
    if numbers.is_empty() || numbers.len().checked_mul(3) != Some(positions.len()) {
        return Err("positions must contain one Cartesian triple per atom");
    }
    if numbers
        .iter()
        .any(|&number| number < 1 || number > max_element)
        || positions.iter().any(|value| !value.is_finite())
        || lattice.is_some_and(|cell| cell.iter().any(|value| !value.is_finite()))
    {
        return Err("invalid atomic number or nonfinite geometry");
    }
    if periodic.iter().any(|&active| active) {
        crate::geometry::periodic_reciprocal(
            lattice.ok_or("periodic structure requires lattice vectors")?,
            periodic,
        )?;
    }
    for first in 0..numbers.len() {
        for second in 0..first {
            let distance2: f64 = (0..3)
                .map(|axis| (positions[3 * first + axis] - positions[3 * second + axis]).powi(2))
                .sum();
            if !distance2.is_finite() {
                return Err("interatomic distance overflow");
            }
            if distance2 < 1.0e-12 {
                return Err("coincident atoms in molecular structure");
            }
        }
    }
    Ok(())
}

pub(crate) unsafe fn require_handles(error: Handle, handles: &[Handle]) -> bool {
    unsafe {
        clear_error(error);
        if handles.iter().any(|handle| handle.is_null()) {
            set_error(error, "required API handle is missing");
            return false;
        }
        true
    }
}

pub(crate) unsafe fn finite_parameters(error: Handle, values: &[f64]) -> bool {
    unsafe {
        clear_error(error);
        if values.iter().any(|value| !value.is_finite()) {
            set_error(error, "parameters must be finite");
            return false;
        }
        true
    }
}

pub(crate) fn validate_output(arrays: &[&[f64]]) -> Result<(), &'static str> {
    if arrays
        .iter()
        .any(|values| values.iter().any(|value| !value.is_finite()))
    {
        return Err("native calculation returned nonfinite results");
    }
    Ok(())
}

pub(crate) struct Error {
    pub message: Option<String>,
}

/// # Safety
/// `handle` must be a live pointer from a common or D3/D4 error constructor, and
/// the returned reference must not alias another reference to the same object.
pub(crate) unsafe fn as_error<'a>(handle: Handle) -> &'a mut Error {
    unsafe { &mut *handle.cast() }
}

pub(crate) unsafe fn set_error(handle: Handle, message: impl Into<String>) {
    unsafe {
        if !handle.is_null() {
            as_error(handle).message = Some(message.into());
        }
    }
}

pub(crate) unsafe fn clear_error(handle: Handle) {
    unsafe {
        if !handle.is_null() {
            as_error(handle).message = None;
        }
    }
}

pub(crate) fn new_error() -> Handle {
    Box::into_raw(Box::new(Error { message: None })).cast()
}

#[no_mangle]
pub extern "C" fn disprs_new_error() -> Handle {
    new_error()
}

#[no_mangle]
pub unsafe extern "C" fn disprs_check_error(error: Handle) -> c_int {
    unsafe { check_error(error) }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_get_error(error: Handle, buffer: *mut c_char, size: *const c_int) {
    unsafe { get_error(error, buffer, size) }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_delete_error(error: *mut Handle) {
    unsafe { delete_error(error) }
}

pub(crate) unsafe fn check_error(handle: Handle) -> c_int {
    unsafe { (handle.is_null() || as_error(handle).message.is_some()) as c_int }
}

pub(crate) unsafe fn get_error(handle: Handle, buffer: *mut c_char, size: *const c_int) {
    unsafe {
        if handle.is_null() || buffer.is_null() {
            return;
        }
        let bytes = as_error(handle).message.as_deref().unwrap_or("").as_bytes();
        let capacity = size.as_ref().copied().unwrap_or(512).max(0) as usize;
        if capacity > 0 {
            let length = bytes.len().min(capacity - 1);
            std::ptr::copy_nonoverlapping(bytes.as_ptr(), buffer.cast(), length);
            *buffer.add(length) = 0;
        }
    }
}

pub(crate) unsafe fn delete_error(handle: *mut Handle) {
    unsafe {
        if let Some(handle) = handle.as_mut() {
            if !handle.is_null() {
                drop(Box::from_raw(handle.cast::<Error>()));
                *handle = null_mut();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_handles() {
        use crate::{d3_ext as d3, d4_ext as d4};
        unsafe {
            let mut error = d3::disprs_d3_new_error();
            let numbers = [6, 8];
            let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
            let displaced = [0.0, 0.0, 0.0, 6.0, 0.0, 0.0];
            let mut neutral = d3::disprs_d3_new_structure(
                error,
                2,
                numbers.as_ptr(),
                positions.as_ptr(),
                null_mut(),
                null_mut(),
            );
            assert!(!neutral.is_null());
            assert_eq!(as_structure(neutral).charge, 0.0);
            d4::disprs_d4_update_structure(error, neutral, displaced.as_ptr(), null_mut());
            assert_eq!(d4::disprs_d4_check_error(error), 0);
            assert_eq!(as_structure(neutral).positions, displaced);
            d4::disprs_d4_delete_structure(&mut neutral);
            assert!(neutral.is_null());
            let mut charged = d4::disprs_d4_new_structure(
                error,
                2,
                numbers.as_ptr(),
                positions.as_ptr(),
                &0.5,
                null_mut(),
                null_mut(),
            );
            assert!(!charged.is_null());
            d3::disprs_d3_update_structure(error, charged, displaced.as_ptr(), null_mut());
            assert_eq!(d3::disprs_d3_check_error(error), 0);
            assert_eq!(as_structure(charged).charge, 0.5);
            assert_eq!(as_structure(charged).positions, displaced);
            d3::disprs_d3_update_structure(error, charged, null_mut(), null_mut());
            assert_eq!(d4::disprs_d4_check_error(error), 1);
            assert_eq!(as_structure(charged).positions, displaced);
            d3::disprs_d3_delete_structure(&mut charged);
            assert!(charged.is_null());
            d4::disprs_d4_delete_error(&mut error);
            assert!(error.is_null());
        }
    }

    #[test]
    fn setup_guard_preserves_errors_and_recovers() {
        unsafe {
            let mut error = new_error();
            let handle: Handle = guard(error, || panic!("constructor failure"));
            assert!(handle.is_null());
            assert_eq!(check_error(error), 1);
            guard(error, || set_error(error, "validation failure"));
            assert_eq!(
                as_error(error).message.as_deref(),
                Some("validation failure")
            );
            assert_eq!(guard(error, || 42), 42);
            assert_eq!(check_error(error), 0);
            assert_eq!(guard::<usize>(null_mut(), || panic!("no error handle")), 0);
            delete_error(&mut error);
        }
    }

    #[test]
    fn evaluation_errors_panics_and_recovery() {
        unsafe {
            let mut error = new_error();
            assert_eq!(evaluate::<()>(error, || Err("native failure")), None);
            assert_eq!(as_error(error).message.as_deref(), Some("native failure"));
            assert_eq!(evaluate::<()>(error, || panic!("injected failure")), None);
            assert_eq!(
                as_error(error).message.as_deref(),
                Some("native calculation panicked")
            );
            assert_eq!(evaluate(error, || Ok(42)), Some(42));
            assert_eq!(check_error(error), 0);
            assert_eq!(evaluate::<()>(null_mut(), || Err("no error handle")), None);
            delete_error(&mut error);
        }
    }

    #[test]
    fn error_lifecycle_and_buffer_bounds() {
        unsafe {
            assert_eq!(check_error(null_mut()), 1);
            set_error(null_mut(), "ignored");
            delete_error(null_mut());
            let mut handle = new_error();
            assert_eq!(check_error(handle), 0);
            let mut buffer = [b'x'; 8];
            get_error(handle, buffer.as_mut_ptr().cast(), &8);
            assert_eq!(buffer[0], 0);
            set_error(handle, "failure");
            assert_eq!(check_error(handle), 1);
            for capacity in [-1, 0] {
                buffer.fill(b'x');
                get_error(handle, buffer.as_mut_ptr().cast(), &capacity);
                assert_eq!(buffer, [b'x'; 8]);
            }
            get_error(handle, buffer.as_mut_ptr().cast(), &1);
            assert_eq!(buffer, [0, b'x', b'x', b'x', b'x', b'x', b'x', b'x']);
            get_error(handle, buffer.as_mut_ptr().cast(), &4);
            assert_eq!(&buffer, b"fai\0xxxx");
            get_error(handle, buffer.as_mut_ptr().cast(), &8);
            assert_eq!(&buffer, b"failure\0");
            get_error(null_mut(), buffer.as_mut_ptr().cast(), &8);
            get_error(handle, null_mut(), &8);
            let mut default_buffer = [b'x'; 512];
            get_error(handle, default_buffer.as_mut_ptr().cast(), std::ptr::null());
            assert_eq!(&default_buffer[..8], b"failure\0");
            assert_eq!(&buffer, b"failure\0");
            delete_error(&mut handle);
            assert!(handle.is_null());
            delete_error(&mut handle);
        }
    }
}
