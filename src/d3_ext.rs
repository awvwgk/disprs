//! C ABI for the D3 and gCP corrections.
//!
//! # Safety
//!
//! Every exported function in this module is `unsafe` and shares one contract:
//!
//! - Error and structure handles are shared with the common and D4 APIs.
//!   Other handles must come from the matching `disprs_d3_new_*`/`disprs_d3_load_*`
//!   constructor in this module, must still be live, and must not have been
//!   passed to its `disprs_d3_delete_*` destructor. A null handle is reported
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

use crate::d3::{self, Damping};
use crate::ffi::{self, as_error, as_structure, Handle, Structure};
use std::convert::TryInto;
use std::ffi::{c_char, c_int, CStr};
use std::ptr::null_mut;

struct Model {
    kind: d3::Model,
    partition: d3::WorkPartition,
    ghosts: Vec<bool>,
    cutoff: d3::RealspaceCutoff,
    ewald: d3::EwaldConfig,
    use_ewald: bool,
}

struct Param {
    damping: Damping,
    atm: Option<d3::Atm>,
}

struct Gcp {
    native: d3::Gcp,
    cutoff: d3::GcpCutoff,
    partition: d3::WorkPartition,
}

/// # Safety
/// `handle` must be a live pointer from `disprs_d3_new_model`, and the returned
/// reference must not alias another reference to the same object.
unsafe fn as_model<'a>(handle: Handle) -> &'a mut Model {
    unsafe { &mut *handle.cast() }
}

/// # Safety
/// `handle` must be a live pointer from a `disprs_d3_*_param` constructor, and
/// the returned reference must not alias another reference to the same object.
unsafe fn as_param<'a>(handle: Handle) -> &'a mut Param {
    unsafe { &mut *handle.cast() }
}

/// # Safety
/// `handle` must be a live pointer from `disprs_d3_load_gcp`, and the returned
/// reference must not alias another reference to the same object.
unsafe fn as_gcp<'a>(handle: Handle) -> &'a mut Gcp {
    unsafe { &mut *handle.cast() }
}

unsafe fn dispersion_inputs(
    error: Handle,
    structure: Handle,
    model: Handle,
    param: Handle,
) -> Option<(
    &'static mut Structure,
    &'static mut Model,
    &'static mut Param,
)> {
    unsafe {
        if !ffi::require_handles(error, &[structure, model, param]) {
            return None;
        }
        let structure = as_structure(structure);
        let model = as_model(model);
        if model.ghosts.len() != structure.numbers.len() {
            ffi::set_error(error, "D3 model and structure atom counts differ");
            return None;
        }
        Some((structure, model, as_param(param)))
    }
}

unsafe fn gcp_inputs(
    error: Handle,
    structure: Handle,
    gcp: Handle,
) -> Option<(&'static mut Structure, &'static mut Gcp)> {
    unsafe {
        if !ffi::require_handles(error, &[structure, gcp]) {
            return None;
        }
        let structure = as_structure(structure);
        let gcp = as_gcp(gcp);
        if gcp.native.emiss.len() != structure.numbers.len() {
            ffi::set_error(error, "gCP parameters and structure atom counts differ");
            return None;
        }
        Some((structure, gcp))
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_version() -> c_int {
    10600
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_has_feature(feature: *mut c_char) -> bool {
    unsafe { !feature.is_null() && CStr::from_ptr(feature).to_bytes() == b"native" }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_error() -> Handle {
    ffi::new_error()
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_check_error(error: Handle) -> c_int {
    unsafe { ffi::check_error(error) }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_error(
    error: Handle,
    buffer: *mut c_char,
    size: *const c_int,
) {
    unsafe {
        ffi::get_error(error, buffer, size);
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_delete_error(error: *mut Handle) {
    unsafe {
        ffi::delete_error(error);
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_structure(
    error: Handle,
    natoms: c_int,
    numbers: *const c_int,
    positions: *const f64,
    lattice: *const f64,
    periodic: *const bool,
) -> Handle {
    unsafe {
        let mut handle = ffi::disprs_new_structure(
            error,
            natoms,
            numbers,
            positions,
            std::ptr::null(),
            lattice,
            periodic,
        );
        if !handle.is_null()
            && as_structure(handle)
                .numbers
                .iter()
                .any(|&number| number > 103)
        {
            ffi::set_error(error, "D3 supports atomic numbers 1 through 103");
            ffi::disprs_delete_structure(&mut handle);
        }
        handle
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_update_structure(
    error: Handle,
    structure: Handle,
    positions: *const f64,
    lattice: *const f64,
) {
    unsafe {
        ffi::disprs_update_structure(error, structure, positions, lattice);
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_delete_structure(handle: *mut Handle) {
    unsafe { ffi::disprs_delete_structure(handle) }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_model(error: Handle, structure_handle: Handle) -> Handle {
    unsafe { disprs_d3_new_model_kind(error, structure_handle, 0) }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_model_kind(
    error: Handle,
    structure_handle: Handle,
    kind: c_int,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            ffi::clear_error(error);
            if structure_handle.is_null() {
                if !error.is_null() {
                    as_error(error).message = Some("molecular structure is missing".into());
                }
                return null_mut();
            }
            let structure = as_structure(structure_handle);
            if structure
                .numbers
                .iter()
                .any(|&number| !(1..=103).contains(&number))
            {
                ffi::set_error(error, "D3 supports atomic numbers 1 through 103");
                return null_mut();
            }
            let kind = match kind {
                0 => d3::Model::D3,
                1 => d3::Model::D3S,
                _ => {
                    ffi::set_error(error, "invalid D3 model");
                    return null_mut();
                }
            };
            if kind == d3::Model::D3S
                && structure
                    .numbers
                    .iter()
                    .any(|number| !(1..=94).contains(number))
            {
                ffi::set_error(error, "D3S supports atomic numbers 1 through 94");
                return null_mut();
            }
            Box::into_raw(Box::new(Model {
                kind,
                partition: d3::WorkPartition::SERIAL,
                ghosts: vec![false; structure.numbers.len()],
                cutoff: d3::RealspaceCutoff::default(),
                use_ewald: false,
                ewald: d3::EwaldConfig {
                    rank: 0,
                    tolerance: 1.0e-4,
                    kcut: 0.0,
                    mesh: 0,
                },
            }))
            .cast()
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_smooth_model(error: Handle, structure: Handle) -> Handle {
    unsafe { disprs_d3_new_model_kind(error, structure, 1) }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_model_realspace_cutoff(
    error: Handle,
    model: Handle,
    disp2: f64,
    disp3: f64,
    cn: f64,
) {
    unsafe {
        disprs_d3_set_model_realspace_cutoff_smooth(error, model, disp2, disp3, cn, 0.0, 0.0);
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_model_ghost_index(
    error: Handle,
    model: Handle,
    ghost: *const c_int,
    nidx: c_int,
) {
    unsafe {
        ffi::guard(error, || {
            if !ffi::require_handles(error, &[model]) || nidx < 0 || (nidx > 0 && ghost.is_null()) {
                ffi::set_error(error, "invalid ghost selection");
                return;
            }
            let model = as_model(model);
            if nidx > 0 && !ghost.is_null() {
                if std::slice::from_raw_parts(ghost, nidx as usize)
                    .iter()
                    .any(|&index| index < 0 || index as usize >= model.ghosts.len())
                {
                    ffi::set_error(error, "ghost index is out of range");
                    return;
                }
                for &index in std::slice::from_raw_parts(ghost, nidx as usize) {
                    if let Some(value) = model.ghosts.get_mut(index as usize) {
                        *value = true;
                    }
                }
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_model_realspace_cutoff_smooth(
    error: Handle,
    model: Handle,
    disp2: f64,
    disp3: f64,
    cn: f64,
    width2: f64,
    width3: f64,
) {
    unsafe {
        ffi::guard(error, || {
            if !ffi::require_handles(error, &[model]) {
                return;
            }
            if [disp2, disp3, cn, width2, width3]
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
            {
                ffi::set_error(error, "cutoffs must be finite and nonnegative");
                return;
            }
            let model = as_model(model);
            model.cutoff = d3::RealspaceCutoff {
                disp2,
                disp3,
                cn,
                width2,
                width3,
            };
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_model_ewald(
    error: Handle,
    model: Handle,
    rank: c_int,
    tolerance: f64,
    kcut: f64,
    mesh: c_int,
) {
    unsafe {
        ffi::guard(error, || {
            if !ffi::require_handles(error, &[model]) {
                return;
            }
            if !tolerance.is_finite() || !kcut.is_finite() {
                ffi::set_error(error, "Ewald controls must be finite");
                return;
            }
            if mesh > 0
                && !(mesh as usize)
                    .checked_next_power_of_two()
                    .and_then(|size| size.checked_pow(3))
                    .and_then(|count| count.checked_mul(std::mem::size_of::<[f64; 2]>()))
                    .is_some_and(|bytes| bytes <= isize::MAX as usize)
            {
                ffi::set_error(error, "Ewald mesh size is not representable");
                return;
            }
            let model = as_model(model);
            if model.kind == d3::Model::D3S {
                ffi::set_error(
                    error,
                    "D3S Fourier summation is not supported; use real-space summation",
                );
                return;
            }
            model.use_ewald = true;
            model.ewald = d3::EwaldConfig {
                rank: rank.max(0) as usize,
                tolerance: if tolerance > 0.0 { tolerance } else { 1.0e-4 },
                kcut: kcut.max(0.0),
                mesh,
            };
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_model_work_partition(
    error: Handle,
    model: Handle,
    part: c_int,
    nparts: c_int,
) {
    unsafe {
        ffi::guard(error, || {
            if !ffi::require_handles(error, &[model]) {
                return;
            }
            let model = as_model(model);
            if let Some(partition) = d3::WorkPartition::new(part, nparts) {
                model.partition = partition;
            } else if !error.is_null() {
                as_error(error).message = Some("invalid work partition".into());
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_model_mpi_comm(error: Handle, model: Handle, comm: c_int) {
    unsafe {
        ffi::guard(error, || {
            let _ = (model, comm);
            if !error.is_null() {
                as_error(error).message = Some("MPI support is not available".into());
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_delete_model(handle: *mut Handle) {
    unsafe {
        if let Some(handle) = handle.as_mut() {
            if !handle.is_null() {
                drop(Box::from_raw(handle.cast::<Model>()));
                *handle = null_mut();
            }
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_load_param(
    error: Handle,
    damping: c_int,
    method: *mut c_char,
    atm: bool,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            ffi::clear_error(error);
            let native = (!method.is_null())
                .then(|| CStr::from_ptr(method).to_str().ok())
                .flatten()
                .and_then(|method| d3::load_named(method, damping, atm));
            if let Some(damping_parameter) = native {
                Box::into_raw(Box::new(Param {
                    damping: damping_parameter,
                    atm: atm.then_some(d3::Atm {
                        s9: 1.0,
                        alpha: 16.0,
                    }),
                }))
                .cast()
            } else {
                if !error.is_null() {
                    as_error(error).message = Some("invalid D3 damping parameters".into());
                }
                null_mut()
            }
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_named_parameters(
    error: Handle,
    damping: c_int,
    method: *mut c_char,
    values: *mut f64,
) {
    unsafe {
        ffi::guard(error, || {
            if values.is_null() {
                if !error.is_null() {
                    as_error(error).message = Some("D3 parameter output is missing".into());
                }
                return;
            }
            let mut handle = disprs_d3_load_param(error, damping, method, false);
            if handle.is_null() {
                return;
            }
            let mut result = [1.0, 1.0, 0.0, 1.0, 1.0, 0.4, 5.0, 14.0, 0.0];
            match as_param(handle).damping {
                Damping::Zero {
                    s6,
                    s8,
                    rs6,
                    rs8,
                    alpha,
                } => {
                    result[0] = s6;
                    result[1] = s8;
                    result[3] = rs6;
                    result[4] = rs8;
                    result[7] = alpha;
                }
                Damping::Rational { s6, s8, a1, a2 } => {
                    result[0] = s6;
                    result[1] = s8;
                    result[5] = a1;
                    result[6] = a2;
                }
                Damping::ModifiedZero {
                    s6,
                    s8,
                    rs6,
                    rs8,
                    alpha,
                    beta,
                } => {
                    result[0] = s6;
                    result[1] = s8;
                    result[3] = rs6;
                    result[4] = rs8;
                    result[7] = alpha;
                    result[8] = beta;
                }
                Damping::OptimizedPower {
                    s6,
                    s8,
                    a1,
                    a2,
                    beta,
                } => {
                    result[0] = s6;
                    result[1] = s8;
                    result[5] = a1;
                    result[6] = a2;
                    result[8] = beta;
                }
                Damping::Cso { s6, a1, a2, a3, a4 } => {
                    result[0] = s6;
                    result[3] = a3;
                    result[4] = a4;
                    result[5] = a1;
                    result[6] = a2;
                }
                Damping::Z { s6, s8, a1 } => {
                    result[0] = s6;
                    result[1] = s8;
                    result[5] = a1;
                }
            }
            std::ptr::copy_nonoverlapping(result.as_ptr(), values, result.len());
            disprs_d3_delete_param(&mut handle);
        });
    }
}

#[cfg(test)]
#[test]
fn setup_panics_are_contained_at_c_boundary() {
    unsafe {
        let mut error = disprs_d3_new_error();
        let numbers = [6, 8];
        let positions = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        let mut structure = disprs_d3_new_structure(
            error,
            2,
            numbers.as_ptr(),
            positions.as_ptr(),
            null_mut(),
            null_mut(),
        );
        let mut model = disprs_d3_new_model(error, structure);
        let mut gcp =
            disprs_d3_load_gcp(error, structure, c"pbeh3c".as_ptr().cast_mut(), null_mut());
        let constructors: [&dyn Fn() -> Handle; 11] = [
            &|| {
                disprs_d3_new_structure(
                    error,
                    2,
                    numbers.as_ptr(),
                    positions.as_ptr(),
                    null_mut(),
                    null_mut(),
                )
            },
            &|| disprs_d3_new_model(error, structure),
            &|| disprs_d3_load_param(error, 1, c"pbe".as_ptr().cast_mut(), true),
            &|| disprs_d3_load_gcp(error, structure, c"pbeh3c".as_ptr().cast_mut(), null_mut()),
            &|| disprs_d3_new_zero_damping(error, 1.0, 1.0, 1.0, 1.0, 1.0, 14.0),
            &|| disprs_d3_new_rational_damping(error, 1.0, 1.0, 1.0, 0.4, 4.0, 14.0),
            &|| disprs_d3_new_mzero_damping(error, 1.0, 1.0, 1.0, 1.0, 1.0, 14.0, 0.1),
            &|| disprs_d3_new_mrational_damping(error, 1.0, 1.0, 1.0, 0.4, 4.0, 14.0),
            &|| disprs_d3_new_optimizedpower_damping(error, 1.0, 1.0, 1.0, 0.4, 4.0, 14.0, 2.0),
            &|| disprs_d3_new_cso_damping(error, 1.0, 1.0, 0.4, 4.0, 1.0, 1.0, 14.0),
            &|| disprs_d3_new_z_damping(error, 1.0, 1.0, 1.0, 0.4, 14.0),
        ];
        for construct in constructors {
            ffi::assert_setup_panic(error, construct);
        }
        let controls: [&dyn Fn(); 13] = [
            &|| disprs_d3_update_structure(error, structure, positions.as_ptr(), null_mut()),
            &|| disprs_d3_set_model_ghost_index(error, model, [0].as_ptr(), 1),
            &|| disprs_d3_set_model_realspace_cutoff(error, model, 7.0, 6.0, 5.0),
            &|| disprs_d3_set_model_realspace_cutoff_smooth(error, model, 7.0, 6.0, 5.0, 1.0, 1.0),
            &|| disprs_d3_set_model_ewald(error, model, 1, 0.01, 1.0, 16),
            &|| disprs_d3_set_model_work_partition(error, model, 1, 2),
            &|| disprs_d3_set_model_mpi_comm(error, model, 0),
            &|| {
                disprs_d3_set_gcp_parameters(
                    error,
                    gcp,
                    2,
                    [2.0; 7].as_ptr(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                )
            },
            &|| disprs_d3_set_gcp_controls(error, gcp, &2.0, &false, &false),
            &|| disprs_d3_set_gcp_realspace_cutoff(error, gcp, 7.0, 6.0),
            &|| disprs_d3_set_gcp_work_partition(error, gcp, 1, 2),
            &|| disprs_d3_set_gcp_mpi_comm(error, gcp, 0),
            &|| {
                disprs_d3_get_gcp_parameters(
                    error,
                    gcp,
                    2,
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                    null_mut(),
                )
            },
        ];
        let gcp_controls = as_gcp(gcp).native.controls();
        for control in controls {
            ffi::assert_setup_panic(error, control);
            assert_eq!(as_structure(structure).positions, positions);
            assert_eq!(as_model(model).ghosts, [false, false]);
            assert!(!as_model(model).use_ewald);
            assert_eq!(
                as_model(model).cutoff.disp2,
                d3::RealspaceCutoff::default().disp2
            );
            assert_eq!(as_gcp(gcp).native.controls(), gcp_controls);
            assert_eq!(as_gcp(gcp).native.sigma, 1.0);
        }
        let mut values = [42.0; 9];
        ffi::assert_setup_panic(error, || {
            disprs_d3_get_named_parameters(
                error,
                1,
                c"pbe".as_ptr().cast_mut(),
                values.as_mut_ptr(),
            )
        });
        assert_eq!(values, [42.0; 9]);
        let mut eta = 42.0;
        ffi::assert_setup_panic(error, || {
            disprs_d3_get_gcp_controls(error, gcp, &mut eta, null_mut(), null_mut())
        });
        assert_eq!(eta, 42.0);
        disprs_d3_get_gcp_controls(error, gcp, &mut eta, null_mut(), null_mut());
        assert_eq!(disprs_d3_check_error(error), 0);
        assert_eq!(eta, gcp_controls.0);
        disprs_d3_delete_gcp(&mut gcp);
        disprs_d3_delete_model(&mut model);
        disprs_d3_delete_structure(&mut structure);
        disprs_d3_delete_error(&mut error);
    }
}

#[cfg(test)]
#[test]
fn setup_domains_preserve_state_and_recover() {
    unsafe {
        let mut error = disprs_d3_new_error();
        let numbers = [6, 8];
        let positions = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        assert!(disprs_d3_new_structure(
            error,
            c_int::MAX,
            numbers.as_ptr(),
            positions.as_ptr(),
            null_mut(),
            null_mut()
        )
        .is_null());
        let mut structure = disprs_d3_new_structure(
            error,
            2,
            numbers.as_ptr(),
            positions.as_ptr(),
            null_mut(),
            null_mut(),
        );
        assert_eq!(disprs_d3_check_error(error), 0);
        disprs_d3_update_structure(
            error,
            structure,
            [0.0, 0.0, 0.0, f64::MAX, 0.0, 0.0].as_ptr(),
            null_mut(),
        );
        assert_eq!(disprs_d3_check_error(error), 1);
        assert_eq!(as_structure(structure).positions, positions);
        let mut model = disprs_d3_new_model(error, structure);
        disprs_d3_set_model_ewald(error, model, 0, 1e-4, 0.0, c_int::MAX);
        assert_eq!(disprs_d3_check_error(error), 1);
        assert!(!as_model(model).use_ewald);
        disprs_d3_set_model_ewald(error, model, 0, 1e-4, 0.0, 17);
        assert_eq!(disprs_d3_check_error(error), 0);
        assert_eq!(as_model(model).ewald.mesh, 17);
        let mut gcp =
            disprs_d3_load_gcp(error, structure, c"pbeh3c".as_ptr().cast_mut(), null_mut());
        assert!(disprs_d3_load_gcp(
            error,
            structure,
            c"pbeh3c".as_ptr().cast_mut(),
            [255u8, 0].as_ptr().cast_mut().cast()
        )
        .is_null());
        assert_eq!(disprs_d3_check_error(error), 1);
        let original = as_gcp(gcp).native.slater.clone();
        let controls = as_gcp(gcp).native.controls();
        for bad in [f64::MAX, f64::NAN, -1.0] {
            disprs_d3_set_gcp_controls(error, gcp, &bad, &false, &true);
            assert_eq!(disprs_d3_check_error(error), 1);
            assert_eq!(as_gcp(gcp).native.slater, original);
            assert_eq!(as_gcp(gcp).native.controls(), controls);
            disprs_d3_set_gcp_controls(error, gcp, null_mut(), null_mut(), null_mut());
            assert_eq!(disprs_d3_check_error(error), 0);
        }
        ffi::set_error(error, "stale");
        disprs_d3_get_gcp_parameters(
            error,
            gcp,
            2,
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
        );
        assert_eq!(disprs_d3_check_error(error), 0);
        ffi::set_error(error, "stale");
        disprs_d3_set_gcp_parameters(
            error,
            gcp,
            2,
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
        );
        assert_eq!(disprs_d3_check_error(error), 0);
        disprs_d3_delete_gcp(&mut gcp);
        disprs_d3_delete_model(&mut model);
        disprs_d3_delete_structure(&mut structure);
        disprs_d3_delete_error(&mut error);
    }
}

#[cfg(test)]
#[test]
fn named_parameter_values_match_pbe() {
    unsafe {
        let mut error = disprs_d3_new_error();
        let method = std::ffi::CString::new("pbe").unwrap();
        let mut values = [0.0; 9];
        disprs_d3_get_named_parameters(error, 1, method.as_ptr().cast_mut(), values.as_mut_ptr());
        assert_eq!(
            values,
            [1.0, 0.7875, 0.0, 1.0, 1.0, 0.4289, 4.4407, 14.0, 0.0]
        );
        assert_eq!(disprs_d3_check_error(error), 0);
        disprs_d3_get_named_parameters(error, 6, method.as_ptr().cast_mut(), values.as_mut_ptr());
        assert_eq!(disprs_d3_check_error(error), 1);
        disprs_d3_delete_error(&mut error);
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_zero_damping(
    error: Handle,
    s6: f64,
    s8: f64,
    s9: f64,
    rs6: f64,
    rs8: f64,
    alp: f64,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::finite_parameters(error, &[s6, s8, s9, rs6, rs8, alp]) {
                return null_mut();
            }
            Box::into_raw(Box::new(Param {
                damping: Damping::Zero {
                    s6,
                    s8,
                    rs6,
                    rs8,
                    alpha: alp,
                },
                atm: (s9.abs() >= f64::EPSILON).then_some(d3::Atm {
                    s9,
                    alpha: alp + 2.0,
                }),
            }))
            .cast()
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_rational_damping(
    error: Handle,
    s6: f64,
    s8: f64,
    s9: f64,
    a1: f64,
    a2: f64,
    alp: f64,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::finite_parameters(error, &[s6, s8, s9, a1, a2, alp]) {
                return null_mut();
            }
            Box::into_raw(Box::new(Param {
                damping: Damping::Rational { s6, s8, a1, a2 },
                atm: (s9.abs() >= f64::EPSILON).then_some(d3::Atm {
                    s9,
                    alpha: alp + 2.0,
                }),
            }))
            .cast()
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_mzero_damping(
    error: Handle,
    s6: f64,
    s8: f64,
    s9: f64,
    rs6: f64,
    rs8: f64,
    alp: f64,
    bet: f64,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::finite_parameters(error, &[s6, s8, s9, rs6, rs8, alp, bet]) {
                return null_mut();
            }
            wrap_param(
                s9,
                alp + 2.0,
                Damping::ModifiedZero {
                    s6,
                    s8,
                    rs6,
                    rs8,
                    alpha: alp,
                    beta: bet,
                },
            )
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_mrational_damping(
    error: Handle,
    s6: f64,
    s8: f64,
    s9: f64,
    a1: f64,
    a2: f64,
    alp: f64,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::finite_parameters(error, &[s6, s8, s9, a1, a2, alp]) {
                return null_mut();
            }
            wrap_param(s9, alp + 2.0, Damping::Rational { s6, s8, a1, a2 })
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_optimizedpower_damping(
    error: Handle,
    s6: f64,
    s8: f64,
    s9: f64,
    a1: f64,
    a2: f64,
    alp: f64,
    bet: f64,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::finite_parameters(error, &[s6, s8, s9, a1, a2, alp, bet]) {
                return null_mut();
            }
            wrap_param(
                s9,
                alp + 2.0,
                Damping::OptimizedPower {
                    s6,
                    s8,
                    a1,
                    a2,
                    beta: bet,
                },
            )
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_cso_damping(
    error: Handle,
    s6: f64,
    s9: f64,
    a1: f64,
    a2: f64,
    a3: f64,
    a4: f64,
    alp: f64,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::finite_parameters(error, &[s6, s9, a1, a2, a3, a4, alp]) {
                return null_mut();
            }
            wrap_param(s9, alp + 2.0, Damping::Cso { s6, a1, a2, a3, a4 })
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_new_z_damping(
    error: Handle,
    s6: f64,
    s8: f64,
    s9: f64,
    a1: f64,
    alp: f64,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::finite_parameters(error, &[s6, s8, s9, a1, alp]) {
                return null_mut();
            }
            wrap_param(s9, alp + 2.0, Damping::Z { s6, s8, a1 })
        })
    }
}

fn wrap_param(s9: f64, alpha: f64, native: Damping) -> Handle {
    Box::into_raw(Box::new(Param {
        damping: native,
        atm: (s9.abs() >= f64::EPSILON).then_some(d3::Atm { s9, alpha }),
    }))
    .cast()
}

#[cfg(test)]
#[test]
fn atm_exponent_matches_upstream_constructors_and_named_parameters() {
    unsafe {
        let constructors = [
            disprs_d3_new_zero_damping(null_mut(), 1.0, 1.0, 1.0, 1.0, 1.0, 14.0),
            disprs_d3_new_rational_damping(null_mut(), 1.0, 1.0, 1.0, 0.4, 4.0, 14.0),
            disprs_d3_new_mzero_damping(null_mut(), 1.0, 1.0, 1.0, 1.0, 1.0, 14.0, 0.1),
            disprs_d3_new_mrational_damping(null_mut(), 1.0, 1.0, 1.0, 0.4, 4.0, 14.0),
            disprs_d3_new_optimizedpower_damping(null_mut(), 1.0, 1.0, 1.0, 0.4, 4.0, 14.0, 2.0),
            disprs_d3_new_cso_damping(null_mut(), 1.0, 1.0, 0.4, 4.0, 1.0, 1.0, 14.0),
            disprs_d3_new_z_damping(null_mut(), 1.0, 1.0, 1.0, 0.4, 14.0),
        ];
        for mut handle in constructors {
            assert_eq!(as_param(handle).atm.unwrap().alpha, 16.0);
            disprs_d3_delete_param(&mut handle);
        }
        let method = std::ffi::CString::new("pbe").unwrap();
        for damping in 0..=5 {
            let mut handle =
                disprs_d3_load_param(null_mut(), damping, method.as_ptr().cast_mut(), true);
            assert!(!handle.is_null());
            assert_eq!(as_param(handle).atm.unwrap().alpha, 16.0);
            disprs_d3_delete_param(&mut handle);
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_delete_param(handle: *mut Handle) {
    unsafe {
        if let Some(handle) = handle.as_mut() {
            if !handle.is_null() {
                drop(Box::from_raw(handle.cast::<Param>()));
                *handle = null_mut();
            }
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_properties(
    error: Handle,
    structure: Handle,
    model: Handle,
    coordination: *mut f64,
    c6: *mut f64,
) {
    unsafe {
        ffi::evaluate(error, || {
            if !ffi::require_handles(error, &[structure, model]) {
                return Ok(());
            }
            let structure = as_structure(structure);
            let model = as_model(model);
            if model.ghosts.len() != structure.numbers.len() {
                return Err("D3 model and structure atom counts differ");
            }
            let (native_cn, native_c6) = model.kind.properties(
                &structure.numbers,
                &structure.positions,
                structure.cell()?,
                model.cutoff.cn,
            )?;
            ffi::validate_output(&[&native_cn, &native_c6])?;
            if !coordination.is_null() {
                std::ptr::copy_nonoverlapping(native_cn.as_ptr(), coordination, native_cn.len());
            }
            if !c6.is_null() {
                std::ptr::copy_nonoverlapping(native_c6.as_ptr(), c6, native_c6.len());
            }
            Ok(())
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_property_response(
    error: Handle,
    structure: Handle,
    model: Handle,
    coordination: *mut f64,
    c6: *mut f64,
    coordination_cartesian: *mut f64,
    coordination_strain: *mut f64,
    c6_cartesian: *mut f64,
    c6_strain: *mut f64,
) {
    unsafe {
        ffi::evaluate(error, || {
            if !ffi::require_handles(error, &[structure, model]) {
                return Ok(());
            }
            let structure = as_structure(structure);
            let model = as_model(model);
            if model.ghosts.len() != structure.numbers.len() {
                return Err("D3 model and structure atom counts differ");
            }
            let result = model.kind.property_response(
                &structure.numbers,
                &structure.positions,
                structure.cell()?,
                model.cutoff.cn,
            )?;
            let outputs = [
                (&result.coordination, coordination),
                (&result.c6, c6),
                (&result.coordination_cartesian, coordination_cartesian),
                (&result.coordination_strain, coordination_strain),
                (&result.c6_cartesian, c6_cartesian),
                (&result.c6_strain, c6_strain),
            ];
            for (values, _) in &outputs {
                ffi::validate_output(&[values])?;
            }
            for (values, output) in outputs {
                if !output.is_null() {
                    std::ptr::copy_nonoverlapping(values.as_ptr(), output, values.len());
                }
            }
            Ok(())
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_dispersion(
    error: Handle,
    structure: Handle,
    model: Handle,
    param: Handle,
    energy: *mut f64,
    gradient: *mut f64,
    virial: *mut f64,
) {
    unsafe {
        ffi::evaluate(error, || {
            let Some((structure, model, param)) = dispersion_inputs(error, structure, model, param)
            else {
                return Ok(());
            };
            if energy.is_null() {
                if !error.is_null() {
                    as_error(error).message = Some("energy output is missing".into());
                }
                return Ok(());
            }
            let damping = param.damping;
            if model.kind == d3::Model::D3S || structure.periodic.iter().any(|&periodic| periodic) {
                if !model.use_ewald {
                    let result = model.kind.dispersion(
                        &structure.numbers,
                        &structure.positions,
                        damping,
                        param.atm,
                        model.cutoff,
                        structure.cell()?,
                        &model.ghosts,
                        model.partition,
                    );
                    match result {
                        Ok((result, _, _)) => {
                            ffi::validate_output(&[
                                &[result.energy],
                                &result.gradient,
                                &result.virial,
                            ])?;
                            *energy = result.energy;
                            if !gradient.is_null() {
                                std::ptr::copy_nonoverlapping(
                                    result.gradient.as_ptr(),
                                    gradient,
                                    result.gradient.len(),
                                );
                            }
                            if !virial.is_null() {
                                std::ptr::copy_nonoverlapping(result.virial.as_ptr(), virial, 9);
                            }
                        }
                        Err(message) => {
                            if !error.is_null() {
                                as_error(error).message = Some(message.into());
                            }
                        }
                    }
                    return Ok(());
                }
                let Some(lattice) = structure
                    .lattice
                    .as_ref()
                    .filter(|_| structure.periodic.iter().all(|&periodic| periodic))
                    .filter(|_| param.atm.is_none() && supports_fourier(damping))
                else {
                    if !error.is_null() {
                        as_error(error).message =
                            Some("unsupported periodic D3 configuration".into());
                    }
                    return Ok(());
                };
                if gradient.is_null() && virial.is_null() {
                    let value = d3::periodic_energy_partitioned_with_cutoff(
                        &structure.numbers,
                        &structure.positions,
                        lattice,
                        damping,
                        model.ewald,
                        model.cutoff,
                        &model.ghosts,
                        model.partition,
                    )?;
                    ffi::validate_output(&[&[value]])?;
                    *energy = value;
                } else {
                    let result = d3::periodic_derivatives_partitioned_with_cutoff(
                        &structure.numbers,
                        &structure.positions,
                        lattice,
                        damping,
                        model.ewald,
                        model.cutoff,
                        &model.ghosts,
                        model.partition,
                    )?;
                    ffi::validate_output(&[&[result.energy], &result.gradient, &result.virial])?;
                    *energy = result.energy;
                    if !gradient.is_null() {
                        std::ptr::copy_nonoverlapping(
                            result.gradient.as_ptr(),
                            gradient,
                            result.gradient.len(),
                        );
                    }
                    if !virial.is_null() {
                        std::ptr::copy_nonoverlapping(
                            result.virial.as_ptr(),
                            virial,
                            result.virial.len(),
                        );
                    }
                }
                return Ok(());
            }
            if gradient.is_null() && virial.is_null() {
                let mut native_energy = d3::energy_partitioned_with_cutoff(
                    &structure.numbers,
                    &structure.positions,
                    damping,
                    model.cutoff,
                    &model.ghosts,
                    model.partition,
                )?;
                if let Some(atm) = param.atm {
                    native_energy += d3::atm_energy_with_cutoff(
                        &structure.numbers,
                        &structure.positions,
                        atm,
                        model.cutoff,
                        &model.ghosts,
                        model.partition,
                    )?;
                }
                ffi::validate_output(&[&[native_energy]])?;
                *energy = native_energy;
            } else {
                let (mut native_energy, mut native_gradient, mut native_virial) =
                    d3::gradient_partitioned_with_cutoff(
                        &structure.numbers,
                        &structure.positions,
                        damping,
                        model.cutoff,
                        &model.ghosts,
                        model.partition,
                    )?;
                if let Some(atm) = param.atm {
                    let (atm_energy, atm_gradient, atm_virial) = d3::atm_gradient_with_cutoff(
                        &structure.numbers,
                        &structure.positions,
                        atm,
                        model.cutoff,
                        &model.ghosts,
                        model.partition,
                    )?;
                    native_energy += atm_energy;
                    for (value, contribution) in native_gradient.iter_mut().zip(atm_gradient) {
                        *value += contribution;
                    }
                    for (value, contribution) in native_virial.iter_mut().zip(atm_virial) {
                        *value += contribution;
                    }
                }
                ffi::validate_output(&[&[native_energy], &native_gradient, &native_virial])?;
                *energy = native_energy;
                if !gradient.is_null() {
                    std::ptr::copy_nonoverlapping(
                        native_gradient.as_ptr(),
                        gradient,
                        native_gradient.len(),
                    );
                }
                if !virial.is_null() {
                    std::ptr::copy_nonoverlapping(
                        native_virial.as_ptr(),
                        virial,
                        native_virial.len(),
                    );
                }
            }
            Ok(())
        });
    }
}

fn supports_fourier(damping: Damping) -> bool {
    matches!(damping, Damping::Rational { .. })
        || matches!(
            damping,
            Damping::Zero { alpha, .. }
                if (alpha - alpha.round()).abs() <= f64::EPSILON.sqrt()
        )
}

#[cfg(test)]
#[test]
fn evaluation_failures_preserve_outputs_and_recover() {
    unsafe {
        let mut error = disprs_d3_new_error();
        let mut structure = disprs_d3_new_structure(
            error,
            2,
            [6, 8].as_ptr(),
            [0.0, 0.0, 0.0, 3.0, 1.0, 0.0].as_ptr(),
            std::ptr::null(),
            std::ptr::null(),
        );
        let mut model = disprs_d3_new_model(error, structure);
        let mut param = disprs_d3_new_rational_damping(error, 1.0, 0.8, 1.0, 0.4, 4.5, 14.0);
        let mut gcp =
            disprs_d3_load_gcp(error, structure, c"pbeh3c".as_ptr().cast_mut(), null_mut());
        assert_eq!(disprs_d3_check_error(error), 0);
        as_structure(structure).positions.pop();
        let mut energy = 42.0;
        let mut gradient = [42.0; 6];
        let mut virial = [42.0; 9];
        let mut hessian = [42.0; 36];
        let mut pair2 = [42.0; 4];
        let mut pair3 = [42.0; 4];
        let mut coordination = [42.0; 2];
        disprs_d3_get_properties(
            error,
            structure,
            model,
            coordination.as_mut_ptr(),
            pair2.as_mut_ptr(),
        );
        assert_eq!(disprs_d3_check_error(error), 1);
        assert_eq!(coordination, [42.0; 2]);
        assert_eq!(pair2, [42.0; 4]);
        for derivatives in [false, true] {
            let grad = if derivatives {
                gradient.as_mut_ptr()
            } else {
                null_mut()
            };
            let sigma = if derivatives {
                virial.as_mut_ptr()
            } else {
                null_mut()
            };
            disprs_d3_get_dispersion(error, structure, model, param, &mut energy, grad, sigma);
            assert_eq!(disprs_d3_check_error(error), 1);
            disprs_d3_get_counterpoise(error, structure, gcp, &mut energy, grad, sigma);
            assert_eq!(disprs_d3_check_error(error), 1);
        }
        disprs_d3_get_pairwise_dispersion(
            error,
            structure,
            model,
            param,
            pair2.as_mut_ptr(),
            pair3.as_mut_ptr(),
        );
        assert_eq!(disprs_d3_check_error(error), 1);
        disprs_d3_get_dispersion_hessian(
            error,
            structure,
            model,
            param,
            &mut energy,
            hessian.as_mut_ptr(),
        );
        assert_eq!(disprs_d3_check_error(error), 1);
        disprs_d3_get_counterpoise_hessian(
            error,
            structure,
            gcp,
            &mut energy,
            hessian.as_mut_ptr(),
        );
        assert_eq!(disprs_d3_check_error(error), 1);
        assert_eq!(energy, 42.0);
        assert_eq!(gradient, [42.0; 6]);
        assert_eq!(virial, [42.0; 9]);
        assert_eq!(hessian, [42.0; 36]);
        assert_eq!(pair2, [42.0; 4]);
        assert_eq!(pair3, [42.0; 4]);
        as_structure(structure).positions.push(0.0);
        disprs_d3_get_dispersion(
            error,
            structure,
            model,
            param,
            &mut energy,
            null_mut(),
            null_mut(),
        );
        assert_eq!(disprs_d3_check_error(error), 0);
        assert!(energy.is_finite() && energy != 42.0);
        disprs_d3_get_properties(
            error,
            structure,
            model,
            coordination.as_mut_ptr(),
            pair2.as_mut_ptr(),
        );
        assert_eq!(disprs_d3_check_error(error), 0);
        assert!(coordination
            .iter()
            .chain(&pair2)
            .all(|value| value.is_finite() && *value != 42.0));
        disprs_d3_delete_gcp(&mut gcp);
        disprs_d3_delete_param(&mut param);
        disprs_d3_delete_model(&mut model);
        disprs_d3_delete_structure(&mut structure);
        disprs_d3_delete_error(&mut error);
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_pairwise_dispersion(
    error: Handle,
    structure: Handle,
    model: Handle,
    param: Handle,
    pair2: *mut f64,
    pair3: *mut f64,
) {
    unsafe {
        ffi::evaluate(error, || {
            let Some((structure, model, param)) = dispersion_inputs(error, structure, model, param)
            else {
                return Ok(());
            };
            if (model.kind == d3::Model::D3S || structure.periodic.iter().any(|&periodic| periodic))
                && !pair2.is_null()
                && !pair3.is_null()
                && !model.use_ewald
            {
                let result = model.kind.dispersion(
                    &structure.numbers,
                    &structure.positions,
                    param.damping,
                    param.atm,
                    model.cutoff,
                    structure.cell()?,
                    &model.ghosts,
                    model.partition,
                );
                match result {
                    Ok((_, native2, native3)) => {
                        ffi::validate_output(&[&native2, &native3])?;
                        std::ptr::copy_nonoverlapping(native2.as_ptr(), pair2, native2.len());
                        std::ptr::copy_nonoverlapping(native3.as_ptr(), pair3, native3.len());
                    }
                    Err(message) => {
                        if !error.is_null() {
                            as_error(error).message = Some(message.into());
                        }
                    }
                }
                return Ok(());
            }
            if !structure.periodic.iter().any(|&periodic| periodic)
                && !pair2.is_null()
                && !pair3.is_null()
            {
                let native2 = d3::pairwise_partitioned_with_cutoff(
                    &structure.numbers,
                    &structure.positions,
                    param.damping,
                    model.cutoff,
                    &model.ghosts,
                    model.partition,
                )?;
                let native3 = if let Some(atm) = param.atm {
                    d3::atm_pairwise_with_cutoff(
                        &structure.numbers,
                        &structure.positions,
                        atm,
                        model.cutoff,
                        &model.ghosts,
                        model.partition,
                    )?
                } else {
                    vec![0.0; native2.len()]
                };
                ffi::validate_output(&[&native2, &native3])?;
                std::ptr::copy_nonoverlapping(native2.as_ptr(), pair2, native2.len());
                std::ptr::copy_nonoverlapping(native3.as_ptr(), pair3, native3.len());
            } else {
                if !error.is_null() {
                    as_error(error).message =
                        Some("unsupported periodic D3 pairwise configuration".into());
                }
            }
            Ok(())
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_load_gcp(
    error: Handle,
    structure: Handle,
    method: *mut c_char,
    basis: *mut c_char,
) -> Handle {
    unsafe {
        ffi::guard(error, || {
            if !ffi::require_handles(error, &[structure]) {
                return null_mut();
            }
            let structure = as_structure(structure);
            let mut names = [None, None];
            for (name, pointer) in names.iter_mut().zip([method, basis]) {
                if !pointer.is_null() {
                    match CStr::from_ptr(pointer).to_str() {
                        Ok(value) => *name = Some(value),
                        Err(_) => {
                            ffi::set_error(error, "gCP method and basis must be valid UTF-8");
                            return null_mut();
                        }
                    }
                }
            }
            let [method, basis] = names;
            let Some(native) = d3::load_gcp(&structure.numbers, method, basis) else {
                if !error.is_null() {
                    as_error(error).message = Some("invalid gCP parameters".into());
                }
                return null_mut();
            };
            Box::into_raw(Box::new(Gcp {
                native,
                cutoff: d3::GcpCutoff::default(),
                partition: d3::WorkPartition::SERIAL,
            }))
            .cast()
        })
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_delete_gcp(handle: *mut Handle) {
    unsafe {
        if let Some(handle) = handle.as_mut() {
            if !handle.is_null() {
                drop(Box::from_raw(handle.cast::<Gcp>()));
                *handle = null_mut();
            }
        }
    }
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "effective numbers are element indices below ELEMENTS"
)]
#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_gcp_parameters(
    error: Handle,
    gcp: Handle,
    natoms: c_int,
    scalars: *mut f64,
    flags: *mut bool,
    zeff: *mut c_int,
    emiss: *mut f64,
    xv: *mut f64,
    slater: *mut f64,
    rvdw: *mut f64,
    rvdw_srb: *mut f64,
) {
    unsafe {
        ffi::guard(error, || {
            if gcp.is_null() || natoms < 0 || as_gcp(gcp).native.emiss.len() != natoms as usize {
                if !error.is_null() {
                    as_error(error).message = Some("invalid gCP parameter buffers".into());
                }
                return;
            }
            let param = &as_gcp(gcp).native;
            let values = [
                param.sigma,
                param.alpha,
                param.beta,
                param.dmp_scal,
                param.dmp_exp,
                param.rscal,
                param.qscal,
            ];
            for (source, target) in [
                (values.as_slice(), scalars),
                (&param.emiss, emiss),
                (&param.virtuals, xv),
                (&param.slater, slater),
            ] {
                if !target.is_null() {
                    std::ptr::copy_nonoverlapping(source.as_ptr(), target, source.len());
                }
            }
            if !flags.is_null() {
                std::ptr::copy_nonoverlapping(
                    [param.damp, param.base, param.srb].as_ptr(),
                    flags,
                    3,
                );
            }
            for atom in 0..natoms as usize {
                if !zeff.is_null() {
                    *zeff.add(atom) = param.effective[atom] as c_int + 1;
                }
                for other in 0..natoms as usize {
                    if !rvdw.is_null() {
                        *rvdw.add(atom * natoms as usize + other) =
                            param.radius(atom, other, false);
                    }
                    if !rvdw_srb.is_null() {
                        *rvdw_srb.add(atom * natoms as usize + other) =
                            param.radius(atom, other, true);
                    }
                }
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_gcp_parameters(
    error: Handle,
    gcp: Handle,
    natoms: c_int,
    scalars: *const f64,
    flags: *const bool,
    zeff: *const c_int,
    emiss: *const f64,
    xv: *const f64,
    slater: *const f64,
    rvdw: *const f64,
    rvdw_srb: *const f64,
) {
    unsafe {
        ffi::guard(error, || {
            if gcp.is_null() || natoms < 0 || as_gcp(gcp).native.emiss.len() != natoms as usize {
                if !error.is_null() {
                    as_error(error).message = Some("invalid gCP parameter buffers".into());
                }
                return;
            }
            let count = natoms as usize;
            let mut candidate = as_gcp(gcp).native.clone();
            if !scalars.is_null() {
                let values = std::slice::from_raw_parts(scalars, 7);
                [
                    candidate.sigma,
                    candidate.alpha,
                    candidate.beta,
                    candidate.dmp_scal,
                    candidate.dmp_exp,
                    candidate.rscal,
                    candidate.qscal,
                ] = values.try_into().unwrap();
            }
            if !flags.is_null() {
                let values = std::slice::from_raw_parts(flags, 3);
                [candidate.damp, candidate.base, candidate.srb] = values.try_into().unwrap();
            }
            if !zeff.is_null() {
                candidate.effective = std::slice::from_raw_parts(zeff, count)
                    .iter()
                    .map(|&value| (value as usize).wrapping_sub(1))
                    .collect();
            }
            for (source, target) in [
                (emiss, &mut candidate.emiss),
                (xv, &mut candidate.virtuals),
                (slater, &mut candidate.slater),
            ] {
                if !source.is_null() {
                    *target = std::slice::from_raw_parts(source, count).to_vec();
                }
            }
            for (source, target) in [
                (rvdw, &mut candidate.rvdw),
                (rvdw_srb, &mut candidate.rvdw_srb),
            ] {
                if !source.is_null() {
                    *target = Some(std::slice::from_raw_parts(source, count * count).to_vec());
                }
            }
            match candidate.validate_parameters() {
                Ok(()) => as_gcp(gcp).native = candidate,
                Err(message) => {
                    if !error.is_null() {
                        as_error(error).message = Some(message.into());
                    }
                }
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_gcp_controls(
    error: Handle,
    gcp: Handle,
    eta: *mut f64,
    base: *mut bool,
    srb: *mut bool,
) {
    unsafe {
        ffi::guard(error, || {
            if gcp.is_null() {
                if !error.is_null() {
                    as_error(error).message = Some("gCP parameters are missing".into());
                }
                return;
            }
            let values = as_gcp(gcp).native.controls();
            if let Some(output) = eta.as_mut() {
                *output = values.0;
            }
            if let Some(output) = base.as_mut() {
                *output = values.1;
            }
            if let Some(output) = srb.as_mut() {
                *output = values.2;
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_gcp_controls(
    error: Handle,
    gcp: Handle,
    eta: *const f64,
    base: *const bool,
    srb: *const bool,
) {
    unsafe {
        ffi::guard(error, || {
            let result = if gcp.is_null() {
                Err("gCP parameters are missing")
            } else {
                let param = &mut as_gcp(gcp).native;
                let defaults = param.controls();
                param.set_controls(
                    eta.as_ref().copied().unwrap_or(defaults.0),
                    base.as_ref().copied().unwrap_or(defaults.1),
                    srb.as_ref().copied().unwrap_or(defaults.2),
                )
            };
            if let Err(message) = result {
                if !error.is_null() {
                    as_error(error).message = Some(message.into());
                }
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_gcp_realspace_cutoff(
    error: Handle,
    gcp: Handle,
    bas: f64,
    srb: f64,
) {
    unsafe {
        ffi::guard(error, || {
            ffi::clear_error(error);
            if gcp.is_null()
                || [bas, srb]
                    .iter()
                    .any(|value| !value.is_finite() || *value < 0.0)
            {
                if !error.is_null() {
                    as_error(error).message =
                        Some("gCP cutoffs must be finite and nonnegative".into());
                }
                return;
            }
            let gcp = as_gcp(gcp);
            gcp.cutoff = d3::GcpCutoff { gcp: bas, srb };
        });
    }
}

#[cfg(test)]
#[test]
fn gcp_parameter_update_is_atomic() {
    unsafe {
        let mut failure = disprs_d3_new_error();
        let mut gcp = Box::into_raw(Box::new(Gcp {
            native: d3::load_gcp(&[6, 8], Some("pbeh3c"), None).unwrap(),
            cutoff: d3::GcpCutoff::default(),
            partition: d3::WorkPartition::SERIAL,
        }))
        .cast();
        let mut scalars = [0.0; 7];
        disprs_d3_get_gcp_parameters(
            failure,
            gcp,
            2,
            scalars.as_mut_ptr(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
        );
        assert_eq!(scalars[0], 1.0);
        scalars[0] = 2.0;
        disprs_d3_set_gcp_parameters(
            failure,
            gcp,
            2,
            scalars.as_ptr(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
        );
        assert_eq!(as_gcp(gcp).native.sigma, 2.0);
        let original = as_gcp(gcp).native.slater.clone();
        for bad in [f64::NAN, f64::INFINITY, -1.0] {
            scalars[0] = 3.0;
            scalars[1] = bad;
            disprs_d3_set_gcp_parameters(
                failure,
                gcp,
                2,
                scalars.as_ptr(),
                null_mut(),
                null_mut(),
                null_mut(),
                null_mut(),
                null_mut(),
                null_mut(),
                null_mut(),
            );
            assert_eq!(disprs_d3_check_error(failure), 1);
            assert_eq!(as_gcp(gcp).native.sigma, 2.0);
            assert_eq!(as_gcp(gcp).native.slater, original);
            as_error(failure).message = None;
        }
        let radii = [5.0, 6.0, 5.0, 5.0];
        disprs_d3_set_gcp_parameters(
            failure,
            gcp,
            2,
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            null_mut(),
            radii.as_ptr(),
            null_mut(),
        );
        assert_eq!(disprs_d3_check_error(failure), 1);
        assert!(as_gcp(gcp).native.rvdw.is_none());
        disprs_d3_delete_gcp(&mut gcp);
        disprs_d3_delete_error(&mut failure);
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_gcp_work_partition(
    error: Handle,
    gcp: Handle,
    part: c_int,
    nparts: c_int,
) {
    unsafe {
        ffi::guard(error, || {
            if !ffi::require_handles(error, &[gcp]) {
                return;
            }
            let gcp = as_gcp(gcp);
            if let Some(partition) = d3::WorkPartition::new(part, nparts) {
                gcp.partition = partition;
            } else if !error.is_null() {
                as_error(error).message = Some("invalid work partition".into());
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_set_gcp_mpi_comm(error: Handle, gcp: Handle, comm: c_int) {
    unsafe {
        ffi::guard(error, || {
            let _ = (gcp, comm);
            if !error.is_null() {
                as_error(error).message = Some("MPI support is not available".into());
            }
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_counterpoise(
    error: Handle,
    structure: Handle,
    gcp: Handle,
    energy: *mut f64,
    gradient: *mut f64,
    virial: *mut f64,
) {
    unsafe {
        ffi::evaluate(error, || {
            let Some((structure, gcp)) = gcp_inputs(error, structure, gcp) else {
                return Ok(());
            };
            if !energy.is_null() {
                let param = &gcp.native;
                if !gradient.is_null() || !virial.is_null() {
                    let result = d3::gcp_derivatives(
                        &structure.numbers,
                        &structure.positions,
                        structure.lattice.as_ref(),
                        structure.periodic,
                        param,
                        gcp.cutoff,
                        gcp.partition,
                    )?;
                    ffi::validate_output(&[&[result.energy], &result.gradient, &result.virial])?;
                    *energy = result.energy;
                    if !gradient.is_null() {
                        std::ptr::copy_nonoverlapping(
                            result.gradient.as_ptr(),
                            gradient,
                            result.gradient.len(),
                        );
                    }
                    if !virial.is_null() {
                        std::ptr::copy_nonoverlapping(
                            result.virial.as_ptr(),
                            virial,
                            result.virial.len(),
                        );
                    }
                } else {
                    let value = d3::gcp_energy(
                        &structure.numbers,
                        &structure.positions,
                        structure.lattice.as_ref(),
                        structure.periodic,
                        param,
                        gcp.cutoff,
                        gcp.partition,
                    )?;
                    ffi::validate_output(&[&[value]])?;
                    *energy = value;
                }
            } else {
                if !error.is_null() {
                    as_error(error).message = Some("energy output is missing".into());
                }
            }
            Ok(())
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_dispersion_hessian(
    error: Handle,
    structure: Handle,
    model: Handle,
    param: Handle,
    energy: *mut f64,
    hessian: *mut f64,
) {
    unsafe {
        ffi::evaluate(error, || {
            let Some((structure, model, param)) = dispersion_inputs(error, structure, model, param)
            else {
                return Ok(());
            };
            if (model.kind == d3::Model::D3S || structure.periodic.iter().any(|&value| value))
                && !model.use_ewald
                && !energy.is_null()
                && !hessian.is_null()
            {
                let result = model.kind.hessian(
                    &structure.numbers,
                    &structure.positions,
                    param.damping,
                    param.atm,
                    model.cutoff,
                    structure.cell()?,
                    &model.ghosts,
                    model.partition,
                );
                match result {
                    Ok((value, matrix)) => {
                        ffi::validate_output(&[&[value], &matrix])?;
                        *energy = value;
                        std::ptr::copy_nonoverlapping(matrix.as_ptr(), hessian, matrix.len());
                    }
                    Err(message) => {
                        if !error.is_null() {
                            as_error(error).message = Some(message.into());
                        }
                    }
                }
                return Ok(());
            }
            if !structure.periodic.iter().any(|&periodic| periodic)
                && !energy.is_null()
                && !hessian.is_null()
            {
                let (mut native_energy, mut native_hessian) = d3::hessian_partitioned_with_cutoff(
                    &structure.numbers,
                    &structure.positions,
                    param.damping,
                    model.cutoff,
                    &model.ghosts,
                    model.partition,
                )?;
                if let Some(atm) = param.atm {
                    let (atm_energy, atm_hessian) = d3::atm_hessian_with_cutoff(
                        &structure.numbers,
                        &structure.positions,
                        atm,
                        model.cutoff,
                        &model.ghosts,
                        model.partition,
                    )?;
                    native_energy += atm_energy;
                    for (value, contribution) in native_hessian.iter_mut().zip(atm_hessian) {
                        *value += contribution;
                    }
                }
                ffi::validate_output(&[&[native_energy], &native_hessian])?;
                *energy = native_energy;
                std::ptr::copy_nonoverlapping(
                    native_hessian.as_ptr(),
                    hessian,
                    native_hessian.len(),
                );
            } else {
                if !error.is_null() {
                    as_error(error).message =
                        Some("unsupported periodic D3 Hessian configuration".into());
                }
            }
            Ok(())
        });
    }
}

#[no_mangle]
pub unsafe extern "C" fn disprs_d3_get_counterpoise_hessian(
    error: Handle,
    structure: Handle,
    gcp: Handle,
    energy: *mut f64,
    hessian: *mut f64,
) {
    unsafe {
        ffi::evaluate(error, || {
            let Some((structure, gcp)) = gcp_inputs(error, structure, gcp) else {
                return Ok(());
            };
            if !energy.is_null() && !hessian.is_null() {
                let param = &gcp.native;
                let result = d3::gcp_hessian(
                    &structure.numbers,
                    &structure.positions,
                    structure.lattice.as_ref(),
                    structure.periodic,
                    param,
                    gcp.cutoff,
                    gcp.partition,
                )?;
                ffi::validate_output(&[&[result.energy], &result.hessian])?;
                *energy = result.energy;
                std::ptr::copy_nonoverlapping(
                    result.hessian.as_ptr(),
                    hessian,
                    result.hessian.len(),
                );
            } else {
                if !error.is_null() {
                    as_error(error).message = Some("gCP Hessian output is missing".into());
                }
            }
            Ok(())
        });
    }
}
