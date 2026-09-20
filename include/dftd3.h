#pragma once

/** @file
 * Source-compatibility names for the supported simple-dftd3 C API subset.
 * Units, buffer sizes, ownership, and errors follow disprs.h; handles cannot
 * be exchanged with a separately loaded upstream library. Aliases do not add
 * validation or functionality. Prefer disprs.h for disprs-specific extensions.
 */

#include <disprs.h>

#ifdef __cplusplus
#include <cstddef>
extern "C" {
#else
#include <stddef.h>
#endif

/** Shared error handle; interchangeable with disprs_error and dftd4_error. */
typedef struct disprs_error_handle *dftd3_error;
/** Shared geometry; interchangeable with disprs_structure and dftd4_structure. */
typedef struct disprs_structure_handle *dftd3_structure;
/** Owned model handle; same contract as disprs_d3_model. */
typedef struct _dftd3_model *dftd3_model;
/** Owned damping handle; same contract as disprs_d3_param. */
typedef struct _dftd3_param *dftd3_param;
/** Owned gCP handle; same contract as disprs_d3_gcp. */
typedef struct _dftd3_gcp *dftd3_gcp;

/** Alias of disprs_d3_get_version. */
#define dftd3_get_version disprs_d3_get_version
/** Alias of disprs_d3_has_feature; only "native" is supported. */
#define dftd3_has_feature disprs_d3_has_feature
/** Alias of disprs_d3_new_error. */
#define dftd3_new_error disprs_d3_new_error
/** Alias of disprs_d3_check_error. */
#define dftd3_check_error disprs_d3_check_error
/** Alias of disprs_d3_get_error. */
#define dftd3_get_error disprs_d3_get_error
/** Alias of disprs_d3_new_structure; coordinates are in bohr. */
#define dftd3_new_structure disprs_d3_new_structure
/** Alias of disprs_d3_update_structure. */
#define dftd3_update_structure disprs_d3_update_structure
/** Alias of disprs_d3_new_model; selects standard D3. */
#define dftd3_new_d3_model disprs_d3_new_model
/** Alias of disprs_d3_set_model_realspace_cutoff. */
#define dftd3_set_model_realspace_cutoff disprs_d3_set_model_realspace_cutoff
/** Alias of disprs_d3_set_model_ghost_index; indices are zero-based and additive. */
#define dftd3_set_model_ghost_index disprs_d3_set_model_ghost_index
/** Alias of disprs_d3_set_model_realspace_cutoff_smooth. */
#define dftd3_set_model_realspace_cutoff_smooth disprs_d3_set_model_realspace_cutoff_smooth
/** Alias of disprs_d3_set_model_ewald; subject to the same Fourier restrictions. */
#define dftd3_set_model_ewald disprs_d3_set_model_ewald
/** Alias of disprs_d3_set_model_work_partition; performs no communication. */
#define dftd3_set_model_work_partition disprs_d3_set_model_work_partition
/** Alias of disprs_d3_set_model_mpi_comm; always reports unsupported MPI. */
#define dftd3_set_model_mpi_comm disprs_d3_set_model_mpi_comm
/** Alias of disprs_d3_new_zero_damping. */
#define dftd3_new_zero_damping disprs_d3_new_zero_damping
/** Alias of disprs_d3_new_rational_damping. */
#define dftd3_new_rational_damping disprs_d3_new_rational_damping
/** Alias of disprs_d3_new_mzero_damping. */
#define dftd3_new_mzero_damping disprs_d3_new_mzero_damping
/** Alias of disprs_d3_new_mrational_damping. */
#define dftd3_new_mrational_damping disprs_d3_new_mrational_damping
/** Alias of disprs_d3_new_optimizedpower_damping. */
#define dftd3_new_optimizedpower_damping disprs_d3_new_optimizedpower_damping
/** Alias of disprs_d3_new_cso_damping. */
#define dftd3_new_cso_damping disprs_d3_new_cso_damping
/** Alias of disprs_d3_new_z_damping. */
#define dftd3_new_z_damping disprs_d3_new_z_damping
/** Alias of disprs_d3_load_gcp; method/basis may be NULL if a supported set resolves. */
#define dftd3_load_gcp_param disprs_d3_load_gcp
/** Alias of disprs_d3_set_gcp_realspace_cutoff. */
#define dftd3_set_gcp_realspace_cutoff disprs_d3_set_gcp_realspace_cutoff
/** Alias of disprs_d3_set_gcp_work_partition. */
#define dftd3_set_gcp_work_partition disprs_d3_set_gcp_work_partition
/** Alias of disprs_d3_set_gcp_mpi_comm; always reports unsupported MPI. */
#define dftd3_set_gcp_mpi_comm disprs_d3_set_gcp_mpi_comm
/** Alias of disprs_d3_get_dispersion; energy is required, gradient/virial optional. */
#define dftd3_get_dispersion disprs_d3_get_dispersion
/** Alias of disprs_d3_get_pairwise_dispersion; both matrices are required. */
#define dftd3_get_pairwise_dispersion disprs_d3_get_pairwise_dispersion
/** Alias of disprs_d3_get_dispersion_hessian; energy and Hessian are required. */
#define dftd3_get_dispersion_hessian disprs_d3_get_dispersion_hessian
/** Alias of disprs_d3_get_counterpoise. */
#define dftd3_get_counterpoise disprs_d3_get_counterpoise
/** Alias of disprs_d3_get_counterpoise_hessian. */
#define dftd3_get_counterpoise_hessian disprs_d3_get_counterpoise_hessian

/** @brief Load named zero-damping parameters.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] method Required read-only NUL-terminated UTF-8 functional name.
 * @param[in] atm Enable ATM with s9=1 and exponent 16; false disables ATM.
 * @return Owned parameters, or NULL for unknown names.
 */
static inline dftd3_param dftd3_load_zero_damping(dftd3_error error, char *method, bool atm) {
    return disprs_d3_load_param(error, DISPRS_D3_ZERO, method, atm);
}

/** Load named BJ parameters (RATIONAL selector).
 * @copydetails dftd3_load_zero_damping
 */
static inline dftd3_param dftd3_load_rational_damping(dftd3_error error, char *method, bool atm) {
    return disprs_d3_load_param(error, DISPRS_D3_RATIONAL, method, atm);
}

/** Load named modified-zero parameters (MODIFIED_ZERO selector).
 * @copydetails dftd3_load_zero_damping
 */
static inline dftd3_param dftd3_load_mzero_damping(dftd3_error error, char *method, bool atm) {
    return disprs_d3_load_param(error, DISPRS_D3_MODIFIED_ZERO, method, atm);
}

/** Load named modified-BJ parameters (MODIFIED_RATIONAL selector).
 * @copydetails dftd3_load_zero_damping
 */
static inline dftd3_param dftd3_load_mrational_damping(dftd3_error error, char *method, bool atm) {
    return disprs_d3_load_param(error, DISPRS_D3_MODIFIED_RATIONAL, method, atm);
}

/** Load named optimized-power parameters (OPTIMIZED_POWER selector).
 * @copydetails dftd3_load_zero_damping
 */
static inline dftd3_param dftd3_load_optimizedpower_damping(dftd3_error error, char *method,
                                                            bool atm) {
    return disprs_d3_load_param(error, DISPRS_D3_OPTIMIZED_POWER, method, atm);
}

/** Load named CSO parameters (CSO selector).
 * @copydetails dftd3_load_zero_damping
 */
static inline dftd3_param dftd3_load_cso_damping(dftd3_error error, char *method, bool atm) {
    return disprs_d3_load_param(error, DISPRS_D3_CSO, method, atm);
}

/** Compatibility loader; no named Z parameters exist. Use dftd3_new_z_damping instead.
 * @param[in,out] error Receives failure status; NULL discards diagnostics.
 * @param[in] method Required read-only NUL-terminated UTF-8 functional name.
 * @param[in] atm Requested ATM flag; no parameter set is available for either value.
 * @return Always NULL.
 */
static inline dftd3_param dftd3_load_z_damping(dftd3_error error, char *method, bool atm) {
    return disprs_d3_load_param(error, DISPRS_D3_Z, method, atm);
}

#define DFTD3_DELETE_ADAPTER(name, compat_type, disprs_type, target)                               \
    static inline void name(compat_type *handle) {                                                 \
        if (handle == NULL)                                                                        \
            return;                                                                                \
        disprs_type value = *handle;                                                               \
        target(&value);                                                                            \
        *handle = NULL;                                                                            \
    }

/** Free shared error state; see disprs_d3_delete_error.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD3_DELETE_ADAPTER(dftd3_delete_error, dftd3_error, disprs_d3_error, disprs_d3_delete_error)
/** Free shared geometry; see disprs_d3_delete_structure.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD3_DELETE_ADAPTER(dftd3_delete_structure, dftd3_structure, disprs_d3_structure,
                     disprs_d3_delete_structure)
/** Free D3/D3S settings; see disprs_d3_delete_model.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD3_DELETE_ADAPTER(dftd3_delete_model, dftd3_model, disprs_d3_model, disprs_d3_delete_model)
/** Free D3 damping parameters; see disprs_d3_delete_param.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD3_DELETE_ADAPTER(dftd3_delete_param, dftd3_param, disprs_d3_param, disprs_d3_delete_param)
/** Free gCP parameters; see disprs_d3_delete_gcp.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD3_DELETE_ADAPTER(dftd3_delete_gcp, dftd3_gcp, disprs_d3_gcp, disprs_d3_delete_gcp)

#undef DFTD3_DELETE_ADAPTER

/** C11 generic destruction using the matching typed destructor.
 * Not usable in C++; use the typed dftd3_delete_* functions with &ptr instead.
 * @param[in,out] ptr Modifiable handle variable (not &ptr), set to NULL;
 * already-null handles allowed.
 */
#define dftd3_delete(ptr)                                                                          \
    _Generic((ptr),                                                                                \
        dftd3_error: dftd3_delete_error,                                                           \
        dftd3_structure: dftd3_delete_structure,                                                   \
        dftd3_model: dftd3_delete_model,                                                           \
        dftd3_param: dftd3_delete_param,                                                           \
        dftd3_gcp: dftd3_delete_gcp)(&(ptr))

#ifdef __cplusplus
}
#endif