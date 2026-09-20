//! C ABI for the D4 and D4S corrections.
//!
//! # Safety
//!
//! Every exported function in this module is `unsafe` and shares one contract:
//!
//! - Error and structure handles are shared with the common and D3 APIs.
//!   Other handles must come from the matching `disprs_d4_new_*`/`disprs_d4_load_*`
//!   constructor in this module, must still be live, and must not have been
//!   passed to its `disprs_d4_delete_*` destructor. A null handle is reported
//!   through the error handle instead of being dereferenced.
//! - A handle must not be used concurrently from two threads, and no two
//!   arguments of a single call may alias the same object.
//! - Input pointers must point to initialized arrays of the documented length;
//!   lengths are derived from the atom count the structure was built with.
//! - Output pointers must be writable for the documented length. Outputs are
//!   written only after the calculation succeeds and its values are finite, so
//!   a failed call leaves caller buffers unmodified.
//! - Destructors take a pointer to the handle, set it to null, and tolerate a
//!   null or already-nulled slot.

use crate::d4;
#[cfg(test)]
use crate::ffi::as_error as error;
use crate::ffi::{self, as_structure as structure, set_error, Handle, Structure};
use std::ffi::{c_char, c_int, CStr};
use std::ptr::null_mut;

struct Model {
    native: d4::Model<'static>,
    ghosts: Vec<bool>,
    fixed_charges: Option<Vec<f64>>,
}
struct Param {
    native: d4::Param,
}
/// # Safety
/// `handle` must be a live pointer from a `disprs_d4_new_*_model` constructor,
/// and the returned reference must not alias another reference to the same object.
unsafe fn model<'a>(handle: Handle) -> &'a mut Model {
    unsafe { &mut *handle.cast() }
}
/// # Safety
/// `handle` must be a live pointer from a `disprs_d4_*_param` constructor, and
/// the returned reference must not alias another reference to the same object.
unsafe fn param<'a>(handle: Handle) -> &'a mut Param {
    unsafe { &mut *handle.cast() }
}
#[no_mangle]
pub extern "C" fn disprs_d4_get_version() -> c_int {
    40200
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_new_error() -> Handle {
    ffi::new_error()
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_check_error(handle: Handle) -> c_int {
    unsafe { ffi::check_error(handle) }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_error(
    handle: Handle,
    buffer: *mut c_char,
    size: *const c_int,
) {
    unsafe {
        ffi::get_error(handle, buffer, size);
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_delete_error(handle: *mut Handle) {
    unsafe {
        ffi::delete_error(handle);
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_new_structure(
    error_handle: Handle,
    natoms: c_int,
    numbers: *const c_int,
    positions: *const f64,
    charge: *const f64,
    lattice: *const f64,
    periodic: *const bool,
) -> Handle {
    unsafe {
        let mut handle = ffi::disprs_new_structure(
            error_handle,
            natoms,
            numbers,
            positions,
            charge,
            lattice,
            periodic,
        );
        if !handle.is_null() {
            let molecule = structure(handle);
            if let Err(message) = d4::validate(&molecule.numbers, &molecule.positions) {
                set_error(error_handle, message);
                ffi::disprs_delete_structure(&mut handle);
            }
        }
        handle
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_update_structure(
    error_handle: Handle,
    structure_handle: Handle,
    positions: *const f64,
    lattice: *const f64,
) {
    unsafe {
        ffi::disprs_update_structure(error_handle, structure_handle, positions, lattice);
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_delete_structure(handle: *mut Handle) {
    unsafe { ffi::disprs_delete_structure(handle) }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_new_model(
    error_handle: Handle,
    structure_handle: Handle,
    kind: c_int,
) -> Handle {
    unsafe { disprs_d4_new_custom_model(error_handle, structure_handle, kind, 3.0, 2.0, 6.0) }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_new_custom_model(
    error_handle: Handle,
    structure_handle: Handle,
    kind: c_int,
    ga: f64,
    gc: f64,
    wf: f64,
) -> Handle {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if structure_handle.is_null() || !matches!(kind, 0 | 1) {
                set_error(error_handle, "invalid D4 model");
                return null_mut();
            }
            let molecule = structure(structure_handle);
            if let Err(message) = d4::validate(&molecule.numbers, &molecule.positions) {
                set_error(error_handle, message);
                return null_mut();
            }
            match d4::Model::custom(kind == 1, ga, gc, wf) {
                Ok(native) => Box::into_raw(Box::new(Model {
                    native,
                    ghosts: vec![false; structure(structure_handle).numbers.len()],
                    fixed_charges: None,
                }))
                .cast(),
                Err(message) => {
                    set_error(error_handle, message);
                    null_mut()
                }
            }
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_set_model_ghost_index(
    error_handle: Handle,
    model_handle: Handle,
    indices: *const c_int,
    count: c_int,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if model_handle.is_null() || count < 0 || (count > 0 && indices.is_null()) {
                set_error(error_handle, "invalid D4 ghost selection");
                return;
            }
            if count == 0 {
                return;
            }
            let model = model(model_handle);
            let indices = std::slice::from_raw_parts(indices, count as usize);
            if indices
                .iter()
                .any(|&index| index < 0 || index as usize >= model.ghosts.len())
            {
                set_error(error_handle, "D4 ghost index is out of range");
                return;
            }
            for &index in indices {
                model.ghosts[index as usize] = true;
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_set_model_work_partition(
    error_handle: Handle,
    model_handle: Handle,
    part: c_int,
    parts: c_int,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if model_handle.is_null() {
                set_error(error_handle, "D4 model is missing");
                return;
            }
            if let Err(message) = model(model_handle).native.set_work_partition(part, parts) {
                set_error(error_handle, message);
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_set_model_ewald(
    error_handle: Handle,
    model_handle: Handle,
    rank: c_int,
    tolerance: f64,
    kcut: f64,
    mesh: c_int,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if model_handle.is_null() {
                set_error(error_handle, "D4 model is missing");
                return;
            }
            if let Err(message) = model(model_handle).native.set_ewald(Some(d4::EwaldConfig {
                rank: rank.max(0) as usize,
                tolerance,
                kcut,
                mesh,
            })) {
                set_error(error_handle, message);
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_set_charge_cutoff(
    error_handle: Handle,
    model_handle: Handle,
    cutoff: f64,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if model_handle.is_null() {
                set_error(error_handle, "D4 model is missing");
                return;
            }
            if let Err(message) = model(model_handle).native.set_charge_cutoff(cutoff) {
                set_error(error_handle, message);
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_set_fixed_charges(
    error_handle: Handle,
    model_handle: Handle,
    charges: *const f64,
    count: c_int,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if model_handle.is_null() {
                set_error(error_handle, "D4 model is missing");
                return;
            }
            let model = model(model_handle);
            if charges.is_null() && count == 0 {
                model.fixed_charges = None;
                return;
            }
            if count < 0 || count as usize != model.ghosts.len() || charges.is_null() {
                set_error(
                    error_handle,
                    "fixed D4 charges must contain one value per atom",
                );
                return;
            }
            let values = std::slice::from_raw_parts(charges, count as usize);
            if let Err(message) = model.native.with_fixed_charges(Some(values)) {
                set_error(error_handle, message);
                return;
            }
            model.fixed_charges = Some(values.to_vec());
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_set_charge_model(
    error_handle: Handle,
    model_handle: Handle,
    charge_model: c_int,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if model_handle.is_null() {
                set_error(error_handle, "D4 model is missing");
                return;
            }
            if let Err(message) = model(model_handle).native.set_charge_model(charge_model) {
                set_error(error_handle, message);
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_set_realspace_cutoff(
    error_handle: Handle,
    model_handle: Handle,
    cn: f64,
    disp2: f64,
    disp3: f64,
    width2: f64,
    width3: f64,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if model_handle.is_null() {
                set_error(error_handle, "D4 model is missing");
                return;
            }
            if let Err(message) = model(model_handle).native.set_cutoff(d4::Cutoff {
                cn,
                disp2,
                disp3,
                width2,
                width3,
            }) {
                set_error(error_handle, message);
            }
        });
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_delete_model(handle: *mut Handle) {
    unsafe {
        delete::<Model>(handle);
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_load_param(
    error_handle: Handle,
    method: *mut c_char,
    atm: bool,
) -> Handle {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            let native = (!method.is_null())
                .then(|| CStr::from_ptr(method).to_str().ok())
                .flatten()
                .and_then(|method| d4::load_param(method, atm));
            if let Some(native) = native {
                Box::into_raw(Box::new(Param { native })).cast()
            } else {
                set_error(error_handle, "invalid D4 damping parameters");
                null_mut()
            }
        })
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_new_rational_damping(
    error_handle: Handle,
    s6: f64,
    s8: f64,
    s9: f64,
    a1: f64,
    a2: f64,
    alp: f64,
) -> Handle {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if ![s6, s8, s9, a1, a2, alp]
                .iter()
                .all(|value| value.is_finite())
            {
                set_error(error_handle, "D4 damping parameters must be finite");
                return null_mut();
            }
            Box::into_raw(Box::new(Param {
                native: d4::Param {
                    s6,
                    s8,
                    s9,
                    a1,
                    a2,
                    alpha: alp,
                },
            }))
            .cast()
        })
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_named_parameters(
    error_handle: Handle,
    method: *const c_char,
    values: *mut f64,
) {
    unsafe {
        disprs_d4_get_named_parameters_s9(error_handle, method, std::ptr::null(), values);
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_named_parameters_s9(
    error_handle: Handle,
    method: *const c_char,
    s9: *const f64,
    values: *mut f64,
) {
    unsafe {
        ffi::guard(error_handle, || {
            ffi::clear_error(error_handle);
            if method.is_null() || values.is_null() {
                set_error(error_handle, "D4 parameter name or output is missing");
                return;
            }
            let s9 = s9.as_ref().copied();
            if s9.is_some_and(|value| !value.is_finite()) {
                set_error(error_handle, "D4 ATM scaling must be finite");
                return;
            }
            let native = CStr::from_ptr(method).to_str().ok().and_then(|name| {
                d4::load_param(name, s9.is_none_or(|value| value.abs() > f64::EPSILON))
            });
            if let Some(native) = native {
                let parameters = [
                    native.s6,
                    native.s8,
                    s9.unwrap_or(native.s9),
                    native.a1,
                    native.a2,
                    native.alpha,
                ];
                std::ptr::copy_nonoverlapping(parameters.as_ptr(), values, parameters.len());
            } else {
                set_error(error_handle, "unknown D4 damping parameters");
            }
        });
    }
}
#[cfg(test)]
#[test]
fn setup_panics_are_contained_at_c_boundary() {
    unsafe {
        let mut failure = disprs_d4_new_error();
        let numbers = [6, 8];
        let positions = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        let mut molecule = disprs_d4_new_structure(
            failure,
            2,
            numbers.as_ptr(),
            positions.as_ptr(),
            null_mut(),
            null_mut(),
            null_mut(),
        );
        let mut handle = disprs_d4_new_model(failure, molecule, 0);
        let constructors: [&dyn Fn() -> Handle; 5] = [
            &|| {
                disprs_d4_new_structure(
                    failure,
                    2,
                    numbers.as_ptr(),
                    positions.as_ptr(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                )
            },
            &|| disprs_d4_new_model(failure, molecule, 0),
            &|| disprs_d4_new_custom_model(failure, molecule, 0, 2.0, 1.0, 5.0),
            &|| disprs_d4_load_param(failure, c"pbe".as_ptr().cast_mut(), true),
            &|| disprs_d4_new_rational_damping(failure, 1.0, 1.0, 1.0, 0.4, 4.0, 16.0),
        ];
        for construct in constructors {
            ffi::assert_setup_panic(failure, construct);
        }
        let original = model(handle).native;
        let controls: [&dyn Fn(); 6] = [
            &|| disprs_d4_update_structure(failure, molecule, positions.as_ptr(), null_mut()),
            &|| disprs_d4_set_model_ghost_index(failure, handle, [0].as_ptr(), 1),
            &|| disprs_d4_set_model_work_partition(failure, handle, 1, 2),
            &|| disprs_d4_set_charge_model(failure, handle, 1),
            &|| disprs_d4_set_charge_cutoff(failure, handle, 9.0),
            &|| disprs_d4_set_realspace_cutoff(failure, handle, 9.0, 8.0, 7.0, 1.0, 1.0),
        ];
        for control in controls {
            ffi::assert_setup_panic(failure, control);
            assert!(model(handle).native == original);
            assert_eq!(model(handle).ghosts, [false, false]);
            assert_eq!(structure(molecule).positions, positions);
        }
        let mut values = [42.0; 6];
        ffi::assert_setup_panic(failure, || {
            disprs_d4_get_named_parameters(failure, c"pbe".as_ptr(), values.as_mut_ptr())
        });
        assert_eq!(values, [42.0; 6]);
        ffi::assert_setup_panic(failure, || {
            disprs_d4_get_named_parameters_s9(failure, c"pbe".as_ptr(), &0.0, values.as_mut_ptr())
        });
        assert_eq!(values, [42.0; 6]);
        disprs_d4_get_named_parameters(failure, c"pbe".as_ptr(), values.as_mut_ptr());
        assert_eq!(disprs_d4_check_error(failure), 0);
        assert!(disprs_d4_new_structure(
            failure,
            c_int::MAX,
            numbers.as_ptr(),
            positions.as_ptr(),
            null_mut(),
            null_mut(),
            null_mut()
        )
        .is_null());
        let cell = [f64::MAX, 0.0, 0.0, 0.0, f64::MAX, 0.0, 0.0, 0.0, f64::MAX];
        for periodic in [[true, false, false], [true, true, false], [true; 3]] {
            assert!(disprs_d4_new_structure(
                failure,
                2,
                numbers.as_ptr(),
                positions.as_ptr(),
                null_mut(),
                cell.as_ptr(),
                periodic.as_ptr()
            )
            .is_null());
            assert_eq!(disprs_d4_check_error(failure), 1);
        }
        disprs_d4_update_structure(failure, molecule, positions.as_ptr(), null_mut());
        assert_eq!(disprs_d4_check_error(failure), 0);
        disprs_d4_delete_model(&mut handle);
        disprs_d4_delete_structure(&mut molecule);
        disprs_d4_delete_error(&mut failure);
    }
}

#[cfg(test)]
#[test]
fn element_support_matches_native_models() {
    unsafe {
        let mut failure = disprs_d4_new_error();
        let positions = [0.0, 0.0, 0.0, 3.0, 1.0, 0.0];
        for element in [0, 104, 111, 119] {
            assert!(disprs_d4_new_structure(
                failure,
                2,
                [element, 8].as_ptr(),
                positions.as_ptr(),
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null()
            )
            .is_null());
            assert_eq!(disprs_d4_check_error(failure), 1);
        }
        for element in 112..=118 {
            let mut structure = disprs_d4_new_structure(
                failure,
                2,
                [element, 8].as_ptr(),
                positions.as_ptr(),
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
            );
            assert!(!structure.is_null());
            disprs_d4_update_structure(failure, structure, positions.as_ptr(), std::ptr::null());
            assert_eq!(disprs_d4_check_error(failure), 0);
            let mut charges = [42.0; 2];
            disprs_get_charges(
                failure,
                structure,
                0,
                charges.as_mut_ptr(),
                null_mut(),
                null_mut(),
            );
            assert_eq!(disprs_d4_check_error(failure), 0);
            assert!(charges.iter().all(|value| value.is_finite()));
            charges.fill(42.0);
            disprs_get_charges(
                failure,
                structure,
                1,
                charges.as_mut_ptr(),
                null_mut(),
                null_mut(),
            );
            assert_eq!(disprs_d4_check_error(failure), 1);
            assert_eq!(charges, [42.0; 2]);
            disprs_d4_delete_structure(&mut structure);
        }
        disprs_d4_delete_error(&mut failure);
    }
}

#[cfg(test)]
#[test]
fn ghost_and_partition_controls_validate_before_mutation() {
    unsafe {
        let mut failure = disprs_d4_new_error();
        let numbers = [6, 8];
        let positions = [0.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        let mut structure = disprs_d4_new_structure(
            failure,
            2,
            numbers.as_ptr(),
            positions.as_ptr(),
            null_mut(),
            null_mut(),
            null_mut(),
        );
        let mut handle = disprs_d4_new_model(failure, structure, 0);
        disprs_d4_set_model_ghost_index(failure, handle, [1].as_ptr(), 1);
        assert_eq!(model(handle).ghosts, [false, true]);
        for indices in [[0, 2], [0, -1]] {
            disprs_d4_set_model_ghost_index(failure, handle, indices.as_ptr(), 2);
            assert_eq!(disprs_d4_check_error(failure), 1);
            assert_eq!(model(handle).ghosts, [false, true]);
            error(failure).message = None;
        }
        for (indices, count) in [(std::ptr::null(), 1), (std::ptr::null(), -1)] {
            disprs_d4_set_model_ghost_index(failure, handle, indices, count);
            assert_eq!(disprs_d4_check_error(failure), 1);
            error(failure).message = None;
        }
        disprs_d4_set_model_ghost_index(failure, handle, std::ptr::null(), 0);
        assert_eq!(disprs_d4_check_error(failure), 0);
        disprs_d4_set_model_work_partition(failure, handle, 1, 3);
        let original = model(handle).native;
        disprs_d4_set_model_work_partition(failure, handle, 3, 3);
        assert_eq!(disprs_d4_check_error(failure), 1);
        assert!(model(handle).native == original);
        error(failure).message = None;
        disprs_d4_set_model_work_partition(failure, null_mut(), 0, 1);
        assert_eq!(disprs_d4_check_error(failure), 1);
        error(failure).message = None;
        disprs_d4_set_model_ghost_index(failure, null_mut(), std::ptr::null(), 0);
        assert_eq!(disprs_d4_check_error(failure), 1);
        disprs_d4_delete_model(&mut handle);
        disprs_d4_delete_structure(&mut structure);
        disprs_d4_delete_error(&mut failure);
    }
}

#[cfg(test)]
#[test]
fn model_controls_validate_before_mutation() {
    unsafe {
        let mut failure = disprs_d4_new_error();
        let numbers = [6, 8];
        let positions = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
        let mut structure = disprs_d4_new_structure(
            failure,
            2,
            numbers.as_ptr(),
            positions.as_ptr(),
            null_mut(),
            null_mut(),
            null_mut(),
        );
        let mut handle = disprs_d4_new_custom_model(failure, structure, 0, 2.0, 1.0, 5.0);
        assert!(!handle.is_null());
        disprs_d4_set_realspace_cutoff(failure, handle, 9.0, 8.0, 7.0, 4.0, 3.0);
        disprs_d4_set_charge_model(failure, handle, 1);
        assert_eq!(disprs_d4_check_error(failure), 0);
        let original = model(handle).native;
        disprs_d4_set_charge_model(failure, handle, 2);
        assert_eq!(disprs_d4_check_error(failure), 1);
        assert!(model(handle).native == original);
        error(failure).message = None;
        disprs_d4_set_charge_model(failure, null_mut(), 0);
        assert_eq!(disprs_d4_check_error(failure), 1);
        error(failure).message = None;
        for bad in [f64::NAN, f64::INFINITY, -1.0] {
            assert!(disprs_d4_new_custom_model(failure, structure, 0, bad, 2.0, 6.0).is_null());
            assert_eq!(disprs_d4_check_error(failure), 1);
            error(failure).message = None;
            disprs_d4_set_realspace_cutoff(failure, handle, 9.0, 8.0, 7.0, bad, 3.0);
            assert_eq!(disprs_d4_check_error(failure), 1);
            assert!(model(handle).native == original);
            error(failure).message = None;
        }
        assert!(disprs_d4_new_custom_model(failure, structure, 1, 3.0, 2.0, 5.0).is_null());
        disprs_d4_delete_model(&mut handle);
        disprs_d4_delete_structure(&mut structure);
        disprs_d4_delete_error(&mut failure);
    }
}

#[cfg(test)]
#[test]
fn explicit_damping_matches_named_and_rejects_nonfinite() {
    unsafe {
        let mut error_handle = disprs_d4_new_error();
        let mut values = [42.0; 6];
        disprs_d4_get_named_parameters_s9(
            error_handle,
            c"dftb(3ob)".as_ptr(),
            &0.0,
            values.as_mut_ptr(),
        );
        assert_eq!(values, [1.0, 0.4727337, 0.0, 0.5467502, 4.4955068, 16.0]);
        let previous = values;
        for bad in [f64::NAN, f64::INFINITY] {
            disprs_d4_get_named_parameters_s9(
                error_handle,
                c"pbe".as_ptr(),
                &bad,
                values.as_mut_ptr(),
            );
            assert_eq!(disprs_d4_check_error(error_handle), 1);
            assert_eq!(values, previous);
        }
        disprs_d4_get_named_parameters(error_handle, c"pbe".as_ptr(), values.as_mut_ptr());
        assert_eq!(disprs_d4_check_error(error_handle), 0);
        assert_eq!(values[2], 1.0);
        for atm in [false, true] {
            let expected = d4::load_param("pbe", atm).unwrap();
            let mut handle = disprs_d4_new_rational_damping(
                error_handle,
                expected.s6,
                expected.s8,
                expected.s9,
                expected.a1,
                expected.a2,
                expected.alpha,
            );
            assert!(!handle.is_null());
            for kind in [d4::Model::D4, d4::Model::D4S] {
                let positions = [0.0, 0.0, 0.0, 0.0, 1.4, 1.0, 0.0, -1.4, 1.0];
                let actual =
                    d4::dispersion(&[8, 1, 1], &positions, 0.0, kind, param(handle).native)
                        .unwrap();
                let reference =
                    d4::dispersion(&[8, 1, 1], &positions, 0.0, kind, expected).unwrap();
                assert_eq!(actual.energy, reference.energy);
                assert_eq!(actual.gradient, reference.gradient);
                assert_eq!(actual.virial, reference.virial);
            }
            disprs_d4_delete_param(&mut handle);
            assert!(handle.is_null());
        }
        assert_eq!(disprs_d4_check_error(error_handle), 0);
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                disprs_d4_new_rational_damping(error_handle, 1.0, bad, 1.0, 0.4, 4.0, 16.0)
                    .is_null()
            );
            assert_eq!(disprs_d4_check_error(error_handle), 1);
        }
        disprs_d4_delete_error(&mut error_handle);
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_delete_param(handle: *mut Handle) {
    unsafe {
        delete::<Param>(handle);
    }
}
unsafe fn inputs(
    error_handle: Handle,
    structure_handle: Handle,
    model_handle: Handle,
) -> Option<(&'static mut Structure, d4::Model<'static>)> {
    unsafe {
        ffi::clear_error(error_handle);
        if structure_handle.is_null() || model_handle.is_null() {
            set_error(error_handle, "D4 structure or model is missing");
            return None;
        }
        let structure = structure(structure_handle);
        let model = model(model_handle);
        match model
            .native
            .with_ghosts(&model.ghosts)
            .with_fixed_charges(model.fixed_charges.as_deref())
        {
            Ok(native) => Some((structure, native)),
            Err(message) => {
                set_error(error_handle, message);
                None
            }
        }
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_properties(
    error_handle: Handle,
    structure_handle: Handle,
    model_handle: Handle,
    coordination_numbers: *mut f64,
    partial_charges: *mut f64,
    c6_coefficients: *mut f64,
    polarizabilities: *mut f64,
) {
    unsafe {
        ffi::evaluate(error_handle, || {
            let Some((structure, kind)) = inputs(error_handle, structure_handle, model_handle)
            else {
                return Ok(());
            };
            let result = if structure.periodic.iter().any(|&periodic| periodic) {
                structure.lattice.as_ref().map_or(
                    Err("periodic D4 requires lattice vectors"),
                    |lattice| {
                        d4::periodic_properties(
                            &structure.numbers,
                            &structure.positions,
                            structure.charge,
                            kind,
                            lattice,
                            structure.periodic,
                        )
                    },
                )
            } else {
                d4::properties(
                    &structure.numbers,
                    &structure.positions,
                    structure.charge,
                    kind,
                )
            };
            match result {
                Ok(result) => {
                    ffi::validate_output(&[
                        &result.coordination,
                        &result.charges,
                        &result.c6,
                        &result.polarizabilities,
                    ])?;
                    copy(&result.coordination, coordination_numbers);
                    copy(&result.charges, partial_charges);
                    copy(&result.c6, c6_coefficients);
                    copy(&result.polarizabilities, polarizabilities);
                }
                Err(message) => set_error(error_handle, message),
            }
            Ok(())
        });
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_property_response(
    error_handle: Handle,
    structure_handle: Handle,
    model_handle: Handle,
    coordination: *mut f64,
    charges: *mut f64,
    c6: *mut f64,
    polarizabilities: *mut f64,
    coordination_cartesian: *mut f64,
    coordination_strain: *mut f64,
    charge_cartesian: *mut f64,
    charge_strain: *mut f64,
    c6_cartesian: *mut f64,
    c6_strain: *mut f64,
    polarizability_cartesian: *mut f64,
    polarizability_strain: *mut f64,
) {
    unsafe {
        ffi::evaluate(error_handle, || {
            let Some((structure, model)) = inputs(error_handle, structure_handle, model_handle)
            else {
                return Ok(());
            };
            let result = d4::property_response(
                &structure.numbers,
                &structure.positions,
                structure.charge,
                model,
                structure
                    .lattice
                    .as_ref()
                    .map(|cell| (cell, structure.periodic)),
            )?;
            let outputs = [
                (&result.properties.coordination, coordination),
                (&result.properties.charges, charges),
                (&result.properties.c6, c6),
                (&result.properties.polarizabilities, polarizabilities),
                (&result.coordination_cartesian, coordination_cartesian),
                (&result.coordination_strain, coordination_strain),
                (&result.charge_cartesian, charge_cartesian),
                (&result.charge_strain, charge_strain),
                (&result.c6_cartesian, c6_cartesian),
                (&result.c6_strain, c6_strain),
                (&result.polarizability_cartesian, polarizability_cartesian),
                (&result.polarizability_strain, polarizability_strain),
            ];
            for (values, _) in &outputs {
                ffi::validate_output(&[values])?;
            }
            for (values, output) in outputs {
                copy(values, output);
            }
            Ok(())
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_get_charges(
    error_handle: Handle,
    structure_handle: Handle,
    charge_model: c_int,
    charges: *mut f64,
    cartesian: *mut f64,
    strain: *mut f64,
) {
    unsafe {
        ffi::evaluate(error_handle, || {
            if !ffi::require_handles(error_handle, &[structure_handle]) {
                return Ok(());
            }
            if charges.is_null() {
                set_error(error_handle, "charge output is missing");
                return Ok(());
            }
            let mut model = d4::Model::D4;
            if let Err(message) = model.set_charge_model(charge_model) {
                set_error(error_handle, message);
                return Ok(());
            }
            let structure = structure(structure_handle);
            match d4::charge_response(
                &structure.numbers,
                &structure.positions,
                structure.charge,
                model,
                structure
                    .lattice
                    .as_ref()
                    .map(|cell| (cell, structure.periodic)),
                !cartesian.is_null() || !strain.is_null(),
            ) {
                Ok(result) => {
                    ffi::validate_output(&[&result.charges, &result.cartesian, &result.strain])?;
                    copy(&result.charges, charges);
                    copy(&result.cartesian, cartesian);
                    copy(&result.strain, strain);
                }
                Err(message) => set_error(error_handle, message),
            }
            Ok(())
        });
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_dispersion(
    error_handle: Handle,
    structure_handle: Handle,
    model_handle: Handle,
    param_handle: Handle,
    energy: *mut f64,
    gradient: *mut f64,
    virial: *mut f64,
) {
    unsafe {
        ffi::evaluate(error_handle, || {
            let Some((structure, kind)) = inputs(error_handle, structure_handle, model_handle)
            else {
                return Ok(());
            };
            if param_handle.is_null() || energy.is_null() {
                set_error(error_handle, "D4 parameter or energy output is missing");
                return Ok(());
            }
            let parameter = param(param_handle).native;
            let periodic = structure.periodic.iter().any(|&value| value);
            if gradient.is_null() && virial.is_null() {
                let result = if periodic {
                    structure.lattice.as_ref().map_or(
                        Err("periodic D4 requires lattice vectors"),
                        |lattice| {
                            d4::periodic_energy(
                                &structure.numbers,
                                &structure.positions,
                                structure.charge,
                                kind,
                                parameter,
                                lattice,
                                structure.periodic,
                            )
                        },
                    )
                } else {
                    d4::energy(
                        &structure.numbers,
                        &structure.positions,
                        structure.charge,
                        kind,
                        parameter,
                    )
                };
                match result {
                    Ok(value) => {
                        ffi::validate_output(&[&[value]])?;
                        *energy = value;
                    }
                    Err(message) => set_error(error_handle, message),
                }
            } else {
                let result = if periodic {
                    structure.lattice.as_ref().map_or(
                        Err("periodic D4 requires lattice vectors"),
                        |lattice| {
                            d4::periodic_dispersion(
                                &structure.numbers,
                                &structure.positions,
                                structure.charge,
                                kind,
                                parameter,
                                lattice,
                                structure.periodic,
                            )
                        },
                    )
                } else {
                    d4::dispersion(
                        &structure.numbers,
                        &structure.positions,
                        structure.charge,
                        kind,
                        parameter,
                    )
                };
                match result {
                    Ok(result) => {
                        ffi::validate_output(&[
                            &[result.energy],
                            &result.gradient,
                            &result.virial,
                        ])?;
                        *energy = result.energy;
                        copy(&result.gradient, gradient);
                        copy(&result.virial, virial);
                    }
                    Err(message) => set_error(error_handle, message),
                }
            }
            Ok(())
        });
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_pairwise_dispersion(
    error_handle: Handle,
    structure_handle: Handle,
    model_handle: Handle,
    param_handle: Handle,
    pair2: *mut f64,
    pair3: *mut f64,
) {
    unsafe {
        ffi::evaluate(error_handle, || {
            let Some((structure, kind)) = inputs(error_handle, structure_handle, model_handle)
            else {
                return Ok(());
            };
            if param_handle.is_null() || pair2.is_null() || pair3.is_null() {
                set_error(error_handle, "D4 pairwise inputs are missing");
                return Ok(());
            }
            let parameter = param(param_handle).native;
            let result = if structure.periodic.iter().any(|&value| value) {
                structure
                    .lattice
                    .as_ref()
                    .ok_or("periodic D4 requires lattice vectors")
                    .and_then(|lattice| {
                        d4::periodic_pairwise(
                            &structure.numbers,
                            &structure.positions,
                            structure.charge,
                            kind,
                            parameter,
                            lattice,
                            structure.periodic,
                        )
                    })
            } else {
                d4::pairwise(
                    &structure.numbers,
                    &structure.positions,
                    structure.charge,
                    kind,
                    parameter,
                )
            };
            match result {
                Ok((two_body, three_body)) => {
                    ffi::validate_output(&[&two_body, &three_body])?;
                    copy(&two_body, pair2);
                    copy(&three_body, pair3);
                }
                Err(message) => set_error(error_handle, message),
            }
            Ok(())
        });
    }
}
#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_dispersion_hessian(
    error_handle: Handle,
    structure_handle: Handle,
    model_handle: Handle,
    param_handle: Handle,
    hessian: *mut f64,
) {
    unsafe {
        ffi::evaluate(error_handle, || {
            let Some((structure, kind)) = inputs(error_handle, structure_handle, model_handle)
            else {
                return Ok(());
            };
            if param_handle.is_null() || hessian.is_null() {
                return Err("D4 parameter or Hessian output is missing");
            }
            let cell = if structure.periodic.iter().any(|&active| active) {
                Some((
                    structure
                        .lattice
                        .as_ref()
                        .ok_or("periodic D4 requires lattice vectors")?,
                    structure.periodic,
                ))
            } else {
                None
            };
            let result = d4::hessian(
                &structure.numbers,
                &structure.positions,
                structure.charge,
                kind,
                param(param_handle).native,
                cell,
            )?;
            ffi::validate_output(&[&result])?;
            copy(&result, hessian);
            Ok(())
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d4_get_numerical_hessian(
    error_handle: Handle,
    structure_handle: Handle,
    model_handle: Handle,
    param_handle: Handle,
    hessian: *mut f64,
) {
    unsafe {
        ffi::evaluate(error_handle, || {
            let Some((structure, kind)) = inputs(error_handle, structure_handle, model_handle)
            else {
                return Ok(());
            };
            if param_handle.is_null() || hessian.is_null() {
                set_error(error_handle, "D4 parameter or Hessian output is missing");
                return Ok(());
            }
            if kind.ewald.is_some() {
                return Err("D4 Fourier Hessians are not supported");
            }
            let parameter = param(param_handle).native;
            let evaluate = |positions: &[f64]| {
                if structure.periodic.iter().any(|&active| active) {
                    structure
                        .lattice
                        .as_ref()
                        .ok_or("periodic D4 requires lattice vectors")
                        .and_then(|lattice| {
                            d4::periodic_dispersion(
                                &structure.numbers,
                                positions,
                                structure.charge,
                                kind,
                                parameter,
                                lattice,
                                structure.periodic,
                            )
                        })
                } else {
                    d4::dispersion(
                        &structure.numbers,
                        positions,
                        structure.charge,
                        kind,
                        parameter,
                    )
                }
            };
            let size = structure.positions.len();
            let mut result = vec![0.0; size * size];
            let mut positions = structure.positions.clone();
            let step = 1.0e-4;
            for column in 0..size {
                positions[column] = structure.positions[column] + step;
                let plus = evaluate(&positions);
                positions[column] = structure.positions[column] - step;
                let minus = evaluate(&positions);
                positions[column] = structure.positions[column];
                match (plus, minus) {
                    (Ok(plus), Ok(minus)) => {
                        for row in 0..size {
                            result[column * size + row] =
                                (plus.gradient[row] - minus.gradient[row]) / (2.0 * step);
                        }
                    }
                    (Err(message), _) | (_, Err(message)) => {
                        set_error(error_handle, message);
                        return Ok(());
                    }
                }
            }
            ffi::validate_output(&[&result])?;
            copy(&result, hessian);
            Ok(())
        });
    }
}

unsafe fn copy(values: &[f64], output: *mut f64) {
    unsafe {
        if !output.is_null() {
            std::ptr::copy_nonoverlapping(values.as_ptr(), output, values.len());
        }
    }
}
unsafe fn delete<T>(handle: *mut Handle) {
    unsafe {
        if let Some(handle) = handle.as_mut() {
            if !handle.is_null() {
                drop(Box::from_raw(handle.cast::<T>()));
                *handle = null_mut();
            }
        }
    }
}
