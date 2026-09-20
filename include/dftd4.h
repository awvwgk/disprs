#pragma once

/** @file
 * Source-compatibility names for the supported dftd4 C API subset.
 * Units, buffer sizes, ownership, and errors follow disprs.h; do not mix handles
 * with a separately loaded upstream library. Prefer disprs.h for native extensions.
 */

#include <disprs.h>

#ifdef __cplusplus
#include <cstddef>
extern "C" {
#else
#include <stddef.h>
#endif

/** Shared error handle; interchangeable with disprs_error and dftd3_error. */
typedef struct disprs_error_handle *dftd4_error;
/** Shared geometry; interchangeable with disprs_structure and dftd3_structure. */
typedef struct disprs_structure_handle *dftd4_structure;
/** Owned model handle; same contract as disprs_d4_model. */
typedef struct _dftd4_model *dftd4_model;
/** Owned damping handle; same contract as disprs_d4_param. */
typedef struct _dftd4_param *dftd4_param;

/** Alias of disprs_d4_get_version. */
#define dftd4_get_version disprs_d4_get_version
/** Alias of disprs_d4_new_error. */
#define dftd4_new_error disprs_d4_new_error
/** Alias of disprs_d4_check_error. */
#define dftd4_check_error disprs_d4_check_error
/** Alias of disprs_d4_get_error. */
#define dftd4_get_error disprs_d4_get_error
/** Alias of disprs_d4_new_structure; charge is optional, coordinates use bohr. */
#define dftd4_new_structure disprs_d4_new_structure
/** Alias of disprs_d4_update_structure. */
#define dftd4_update_structure disprs_d4_update_structure
/** Alias of disprs_d4_new_rational_damping. */
#define dftd4_new_rational_damping disprs_d4_new_rational_damping
/** Alias of disprs_d4_load_param. */
#define dftd4_load_rational_damping disprs_d4_load_param
/** Alias of disprs_d4_get_properties; every output is optional. */
#define dftd4_get_properties disprs_d4_get_properties
/** Alias of disprs_d4_get_dispersion; energy required, gradient/virial optional. */
#define dftd4_get_dispersion disprs_d4_get_dispersion
/** Alias of disprs_d4_get_pairwise_dispersion; both matrices are required. */
#define dftd4_get_pairwise_dispersion disprs_d4_get_pairwise_dispersion
/** Alias of disprs_d4_get_numerical_hessian; fixed 1e-4 bohr displacement. */
#define dftd4_get_numerical_hessian disprs_d4_get_numerical_hessian

/** Create D4 settings; defaults follow disprs_d4_new_model with DISPRS_D4.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] mol Live shared geometry with D4-supported elements; not retained.
 * @return Owned model, or NULL on error.
 */
static inline dftd4_model dftd4_new_d4_model(dftd4_error error, dftd4_structure mol) {
    return disprs_d4_new_model(error, mol, DISPRS_D4);
}

/** Create D4S settings; defaults follow disprs_d4_new_model with DISPRS_D4S.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] mol Live shared geometry with D4-supported elements; not retained.
 * @return Owned model, or NULL on error.
 */
static inline dftd4_model dftd4_new_d4s_model(dftd4_error error, dftd4_structure mol) {
    return disprs_d4_new_model(error, mol, DISPRS_D4S);
}

/** Create custom D4 settings; other defaults follow disprs_d4_new_custom_model.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] mol Live shared geometry with D4-supported elements; not retained.
 * @param[in] ga Finite positive charge-scaling amplitude in atomic-unit conventions.
 * @param[in] gc Finite positive hardness multiplier in the charge-scaling exponent.
 * @param[in] wf Finite positive Gaussian CN width.
 * @return Owned model, or NULL on error.
 */
static inline dftd4_model dftd4_custom_d4_model(dftd4_error error, dftd4_structure mol, double ga,
                                                double gc, double wf) {
    return disprs_d4_new_custom_model(error, mol, DISPRS_D4, ga, gc, wf);
}

/** Create custom D4S settings with wf=6; other defaults follow disprs_d4_new_custom_model.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] mol Live shared geometry with D4-supported elements; not retained.
 * @param[in] ga Finite positive charge-scaling amplitude in atomic-unit conventions.
 * @param[in] gc Finite positive hardness multiplier in the charge-scaling exponent.
 * @return Owned model, or NULL on error.
 */
static inline dftd4_model dftd4_custom_d4s_model(dftd4_error error, dftd4_structure mol, double ga,
                                                 double gc) {
    return disprs_d4_new_custom_model(error, mol, DISPRS_D4S, ga, gc, 6.0);
}

/** Set hard cutoffs; both switching widths reset to zero.
 * Argument order differs from disprs_d4_set_realspace_cutoff.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings to update.
 * @param[in] disp2 Finite nonnegative two-body radius in bohr.
 * @param[in] disp3 Finite nonnegative ATM radius in bohr.
 * @param[in] cn Finite nonnegative CN radius in bohr.
 */
static inline void dftd4_set_model_realspace_cutoff(dftd4_error error, dftd4_model model,
                                                    double disp2, double disp3, double cn) {
    disprs_d4_set_realspace_cutoff(error, model, cn, disp2, disp3, 0.0, 0.0);
}

/** Set cutoffs with switching as in disprs_d4_set_realspace_cutoff; argument order differs.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings to update.
 * @param[in] disp2 Finite nonnegative two-body radius in bohr.
 * @param[in] disp3 Finite nonnegative ATM radius in bohr.
 * @param[in] cn Finite nonnegative CN radius in bohr.
 * @param[in] width2 Finite nonnegative two-body switching width in bohr, capped at disp2.
 * @param[in] width3 Finite nonnegative ATM switching width in bohr, capped at disp3.
 */
static inline void dftd4_set_model_realspace_cutoff_smooth(dftd4_error error, dftd4_model model,
                                                           double disp2, double disp3, double cn,
                                                           double width2, double width3) {
    disprs_d4_set_realspace_cutoff(error, model, cn, disp2, disp3, width2, width3);
}

#define DFTD4_DELETE_ADAPTER(name, compat_type, disprs_type, target)                               \
    static inline void name(compat_type *handle) {                                                 \
        if (handle == NULL)                                                                        \
            return;                                                                                \
        disprs_type value = *handle;                                                               \
        target(&value);                                                                            \
        *handle = NULL;                                                                            \
    }

/** Free shared error state; see disprs_d4_delete_error.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD4_DELETE_ADAPTER(dftd4_delete_error, dftd4_error, disprs_d4_error, disprs_d4_delete_error)
/** Free shared geometry; see disprs_d4_delete_structure.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD4_DELETE_ADAPTER(dftd4_delete_structure, dftd4_structure, disprs_d4_structure,
                     disprs_d4_delete_structure)
/** Free D4/D4S settings; see disprs_d4_delete_model.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD4_DELETE_ADAPTER(dftd4_delete_model, dftd4_model, disprs_d4_model, disprs_d4_delete_model)
/** Free D4 damping parameters; see disprs_d4_delete_param.
 * @param[in,out] handle Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
DFTD4_DELETE_ADAPTER(dftd4_delete_param, dftd4_param, disprs_d4_param, disprs_d4_delete_param)
#undef DFTD4_DELETE_ADAPTER

/** C11 generic destruction using the matching typed destructor.
 * Not usable in C++; use the typed dftd4_delete_* functions with &ptr instead.
 * @param[in,out] ptr Modifiable handle variable (not &ptr), set to NULL;
 * already-null handles allowed.
 */
#define dftd4_delete(ptr)                                                                          \
    _Generic((ptr),                                                                                \
        dftd4_error: dftd4_delete_error,                                                           \
        dftd4_structure: dftd4_delete_structure,                                                   \
        dftd4_model: dftd4_delete_model,                                                           \
        dftd4_param: dftd4_delete_param)(&(ptr))

#ifdef __cplusplus
}
#endif