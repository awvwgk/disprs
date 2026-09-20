#pragma once

/** @file
 * Native dispersion and gCP C API.
 *
 * Handles must be live, of the declared type, and created by this library.
 * Do not share a handle between concurrent calls or alias arguments/outputs.
 * Constructors return owned handles (NULL on failure); delete each exactly once.
 * Destructors accept NULL slots or NULL handles and set a live slot to NULL.
 * Error and structure handles are interchangeable across the common, D3, D4,
 * and gCP APIs. Their D3/D4 lifecycle names remain compatibility aliases.
 * Models and damping parameters are model-family-specific: never pass D3
 * handles to D4 or vice versa. gCP parameters are a separate handle type.
 * Input arrays/strings are borrowed for the call; retained data is copied.
 * Strings are NUL-terminated UTF-8, read-only even where declared char *.
 *
 * An error argument receives call status; NULL discards diagnostics. Check it
 * after each call: computational/configuration calls reset the previous status.
 * Do not inspect calculation outputs after failure; calculation buffers are
 * preserved on reported errors. Invalid/dangling pointers cannot be validated.
 * Non-NULL buffers must have the documented capacity; outputs must not overlap.
 *
 * N is the structure atom count. Coordinates are 3*N doubles, consecutive xyz
 * triples, in bohr. Lattices are 9 doubles, three consecutive lattice vectors
 * (column-major). Energies are hartree; gradients are dE/dR, not forces, in
 * hartree/bohr. Virials are 9 strain derivatives in hartree, row + 3*column.
 * Hessians use column-major (3*N,3*N), in hartree/bohr^2.
 * Property Jacobians are property-major: coordinate 3*atom+axis or strain
 * row+3*column varies fastest. Strain transforms both coordinates and lattice.
 * Derivatives assume hard-cutoff/image-selection boundaries are not crossed.
 */

#ifdef __cplusplus
#include <cstdbool>
extern "C" {
#else
#include <stdbool.h>
#endif

/** Native library major version. */
#define DISPRS_VERSION_MAJOR 0
/** Native library minor version. */
#define DISPRS_VERSION_MINOR 1
/** Native library patch version. */
#define DISPRS_VERSION_PATCH 0

/** Owned error state shared by D3, D4, and gCP; query after each operation. */
typedef void *disprs_error;
/** Owned geometry shared by D3, D4, and gCP; D3/gCP ignore its total charge.
 * Atomic identities, total charge, and periodic flags are fixed at creation.
 * A selected model must support the structure's elements.
 */
typedef void *disprs_structure;
/** Compatibility alias of disprs_error; interchangeable with D4 error handles. */
typedef disprs_error disprs_d3_error;
/** Compatibility alias of disprs_structure; D3 creation sets total charge to zero. */
typedef disprs_structure disprs_d3_structure;
/** Owned D3/D3S settings for a fixed atom count; not a D4 model. */
typedef void *disprs_d3_model;
/** Owned D3 damping and optional ATM parameters; not interchangeable with D4. */
typedef void *disprs_d3_param;
/** Owned gCP parameters tied to ordered atomic numbers; distinct from dispersion
 * model/damping handles. Accepts the shared geometry type, independent of D3/D4.
 */
typedef void *disprs_d3_gcp;
/** Compatibility alias of disprs_error; interchangeable with D3 error handles. */
typedef disprs_error disprs_d4_error;
/** Compatibility alias of disprs_structure; interchangeable with D3 geometries. */
typedef disprs_structure disprs_d4_structure;
/** Owned D4/D4S settings, including charge controls; not a D3 model. */
typedef void *disprs_d4_model;
/** Owned D4 BJ/ATM parameters; not D3 parameters (different ATM exponent convention). */
typedef void *disprs_d4_param;

/** D3 two-body damping selector; modified BJ uses the BJ functional form. */
typedef enum {
    DISPRS_D3_ZERO,              /**< Original zero damping. */
    DISPRS_D3_RATIONAL,          /**< Becke-Johnson (BJ) damping. */
    DISPRS_D3_MODIFIED_ZERO,     /**< Modified zero damping. */
    DISPRS_D3_MODIFIED_RATIONAL, /**< Reparameterized BJ damping. */
    DISPRS_D3_OPTIMIZED_POWER,   /**< Optimized-power damping. */
    DISPRS_D3_CSO,               /**< C6-only sigmoid-scaled damping. */
    DISPRS_D3_Z                  /**< C6-dependent Z damping; no named parameter sets. */
} disprs_d3_damping;

/** D3 interpolation model. */
typedef enum {
    DISPRS_D3, /**< Standard interpolation, atomic numbers 1..103. */
    DISPRS_D3S /**< Pair-width interpolation, atomic numbers 1..94. */
} disprs_d3_model_kind;
/** D4 interpolation model. */
typedef enum {
    DISPRS_D4, /**< Standard coordination-number weights. */
    DISPRS_D4S /**< Pair-dependent coordination-number weights. */
} disprs_d4_model_kind;
/** Charge-equilibration model. */
typedef enum {
    DISPRS_EEQ,  /**< EEQ charge equilibration. */
    DISPRS_EEQBC /**< Bond-capacity EEQ, atomic numbers 1..103. */
} disprs_d4_charge_model;

/** Return a static version string; never modify or free the pointer. */
const char *disprs_get_version(void);
/** Allocate shared, initially clear error state; release with disprs_delete_error. */
disprs_error disprs_new_error(void);
/** Query status without clearing it.
 * @param[in] error Error handle, or NULL.
 * @return 1 for a pending error or NULL handle, otherwise 0.
 */
int disprs_check_error(disprs_error error);
/** Copy error text without clearing status; truncates with NUL termination.
 * NULL error/buffer or nonpositive capacity writes nothing.
 * @param[in] error Error handle to query, or NULL.
 * @param[out] buffer Writable message buffer, or NULL to skip output.
 * @param[in] buffer_size Capacity in bytes, including NUL; NULL means 512 bytes.
 */
void disprs_get_error(disprs_error error, char *buffer, const int *buffer_size);
/** Free an error from any common, D3, or D4 constructor.
 * @param[in,out] error Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
void disprs_delete_error(disprs_error *error);
/** Copy a geometry into an owned shared structure.
 * Model-specific element limits are checked when constructing/evaluating a model.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] natoms Atom count N>0, with representable array sizes.
 * @param[in] numbers Required N atomic numbers in 1..118.
 * @param[in] positions Required 3*N finite xyz values in bohr; separations >=1e-6 bohr.
 * @param[in] charge Optional finite total charge in electron-charge units; NULL means
 * zero. D3/gCP ignore charge.
 * @param[in] lattice Optional 9 finite column-major cell values in bohr;
 * active vectors must be independent.
 * @param[in] periodic Optional 3 active-vector flags; NULL means all with a cell,
 * otherwise none. Active periodicity requires lattice.
 * @return Owned structure, or NULL on invalid input.
 */
disprs_structure disprs_new_structure(disprs_error error, int natoms, const int *numbers,
                                      const double *positions, const double *charge,
                                      const double *lattice, const bool *periodic);
/** Replace geometry; identities, count, charge, and periodic flags remain fixed.
 * Invalid input preserves geometry; constraints follow disprs_new_structure.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] structure Live common, D3, or D4 structure to update.
 * @param[in] positions Required 3*N finite xyz values in bohr; separations >=1e-6 bohr.
 * @param[in] lattice Optional 9 finite column-major cell values in bohr;
 * active vectors must be independent. NULL retains the cell.
 */
void disprs_update_structure(disprs_error error, disprs_structure structure,
                             const double *positions, const double *lattice);
/** Free a common, D3, or D4 structure; all three destructor names are equivalent.
 * @param[in,out] structure Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
void disprs_delete_structure(disprs_structure *structure);
/** Return D3 compatibility version 10600 (1.6.0), not the native library version. */
int disprs_d3_get_version(void);
/** Test a feature name.
 * @param[in] feature Read-only NUL-terminated UTF-8 name, or NULL.
 * @return True only for "native"; NULL/unknown returns false.
 */
bool disprs_d3_has_feature(char *feature);

/** Allocate an initially clear error state; release with disprs_d3_delete_error. */
disprs_d3_error disprs_d3_new_error(void);
/** @copydoc disprs_check_error */
int disprs_d3_check_error(disprs_d3_error error);
/** @copydoc disprs_get_error */
void disprs_d3_get_error(disprs_d3_error error, char *buffer, const int *buffer_size);
/** @copydoc disprs_delete_error */
void disprs_d3_delete_error(disprs_d3_error *error);

/** Copy a geometry into an owned structure with zero total charge.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] natoms Atom count N>0, with representable array sizes.
 * @param[in] numbers Required N atomic numbers in 1..103.
 * @param[in] positions Required 3*N finite xyz values in bohr; separations >=1e-6 bohr.
 * @param[in] lattice Optional 9 finite column-major cell values in bohr;
 * active vectors must be independent.
 * @param[in] periodic Optional 3 active-vector flags; NULL means all with a cell,
 * otherwise none. Active periodicity requires lattice.
 * @return Owned structure, or NULL on invalid input; no input buffer is retained.
 */
disprs_d3_structure disprs_d3_new_structure(disprs_d3_error error, int natoms, const int *numbers,
                                            const double *positions, const double *lattice,
                                            const bool *periodic);
/** @copydoc disprs_update_structure */
void disprs_d3_update_structure(disprs_d3_error error, disprs_d3_structure structure,
                                const double *positions, const double *lattice);
/** @copydoc disprs_delete_structure */
void disprs_d3_delete_structure(disprs_d3_structure *structure);

/** Create standard D3 settings; defaults follow disprs_d3_new_model_kind.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared structure with Z=1..103; not retained.
 * @return Owned D3 model, or NULL on error.
 */
disprs_d3_model disprs_d3_new_model(disprs_d3_error error, disprs_d3_structure structure);
/** Create D3 or D3S settings.
 * Defaults: CN/pair/triplet cutoffs 40/60/40 bohr, no switching or ghosts,
 * serial work, real-space summation. Model and evaluation structure must have N atoms.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared structure; Z=1..103 for D3, Z=1..94 for D3S.
 * The structure handle is not retained.
 * @param[in] kind Interpolation selector: DISPRS_D3 or DISPRS_D3S.
 * @return Owned model, or NULL on error.
 */
disprs_d3_model disprs_d3_new_model_kind(disprs_d3_error error, disprs_d3_structure structure,
                                         disprs_d3_model_kind kind);
/** Create D3S settings; defaults follow disprs_d3_new_model_kind.
 * Uses standard D3 damping; supports molecular/periodic real space, not Fourier.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared structure with Z=1..94; not retained.
 * @return Owned D3S model, or NULL on error.
 */
disprs_d3_model disprs_d3_new_smooth_model(disprs_d3_error error, disprs_d3_structure structure);
/** Free D3/D3S model settings.
 * @param[in,out] model Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
void disprs_d3_delete_model(disprs_d3_model *model);
/** Set hard cutoffs; resets both switching widths to zero.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D3/D3S settings to update.
 * @param[in] dispersion2 Finite nonnegative two-body radius in bohr.
 * @param[in] dispersion3 Finite nonnegative ATM radius in bohr.
 * @param[in] coordination_number Finite nonnegative CN radius in bohr.
 */
void disprs_d3_set_model_realspace_cutoff(disprs_d3_error error, disprs_d3_model model,
                                          double dispersion2, double dispersion3,
                                          double coordination_number);
/** Add atoms to the ghost mask; duplicates are allowed.
 * Ghosts remain in CN calculations but not explicit dispersion pairs/triplets.
 * Does not clear previous selections; invalid input preserves the mask.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D3/D3S settings for N atoms.
 * @param[in] ghost Array of count zero-based indices in [0,N); NULL only if count=0.
 * @param[in] count Number of indices, >=0; zero leaves the mask unchanged.
 */
void disprs_d3_set_model_ghost_index(disprs_d3_error error, disprs_d3_model model, const int *ghost,
                                     int count);
/** Set cutoffs and switching widths. A width of zero gives a hard cutoff;
 * 0<width<radius gives quintic switching; width>=radius disables switching in D3.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D3/D3S settings to update.
 * @param[in] dispersion2 Finite nonnegative two-body radius in bohr.
 * @param[in] dispersion3 Finite nonnegative ATM radius in bohr.
 * @param[in] coordination_number Finite nonnegative CN radius in bohr.
 * @param[in] width2 Finite nonnegative two-body switching width in bohr.
 * @param[in] width3 Finite nonnegative ATM switching width in bohr.
 */
void disprs_d3_set_model_realspace_cutoff_smooth(disprs_d3_error error, disprs_d3_model model,
                                                 double dispersion2, double dispersion3,
                                                 double coordination_number, double width2,
                                                 double width3);
/** Enable full-3D Fourier two-body summation.
 * No disable flag: create a new model to restore real space. Periodic evaluation requires
 * full periodicity and BJ/modified-BJ or zero damping with an integer exponent,
 * with no ATM. Only the CN cutoff applies; pair switching is not applied.
 * Periodic pair matrices and Hessians are unavailable. Molecular D3 ignores this setting.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D3 settings; D3S is unsupported.
 * @param[in] rank Reference-C6 expansion cap; <=0 selects by tolerance.
 * @param[in] tolerance Finite factorization tolerance; <=0 means 1e-4.
 * @param[in] reciprocal_cutoff Finite reciprocal radius in bohr^-1; <=0 is automatic.
 * @param[in] mesh Negative for direct summation, zero for automatic SPME, positive
 * for SPME rounded up to a power of two; storage must be representable.
 */
void disprs_d3_set_model_ewald(disprs_d3_error error, disprs_d3_model model, int rank,
                               double tolerance, double reciprocal_cutoff, int mesh);
/** Assign model work to one worker.
 * No communication occurs. Sum energy/derivative/pair outputs over all workers;
 * unpartitioned properties remain complete. (part,parts)=(0,1) restores serial work.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D3/D3S settings to update.
 * @param[in] part Zero-based worker index, 0<=part<parts.
 * @param[in] parts Positive worker count.
 */
void disprs_d3_set_model_work_partition(disprs_d3_error error, disprs_d3_model model, int part,
                                        int parts);
/** Compatibility stub: always reports "MPI support is not available" via error.
 * @param[in,out] error Receives unsupported-MPI status; NULL discards diagnostics.
 * @param[in] model Ignored; no state changes.
 * @param[in] communicator Ignored Fortran MPI communicator token.
 */
void disprs_d3_set_model_mpi_comm(disprs_d3_error error, disprs_d3_model model, int communicator);

/** Load named D3 damping parameters.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] damping Selector ZERO through CSO; Z has no named sets.
 * @param[in] method Required read-only NUL-terminated UTF-8 functional name.
 * @param[in] atm Enable ATM with s9=1 and exponent 16; false disables ATM.
 * @return Owned parameters, or NULL for unknown names/selectors.
 */
disprs_d3_param disprs_d3_load_param(disprs_d3_error error, disprs_d3_damping damping, char *method,
                                     bool atm);
/** Query named parameters without allocating; output is unchanged on error.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] damping Selector ZERO through CSO; Z has no named sets.
 * @param[in] method Required read-only NUL-terminated UTF-8 functional name.
 * @param[out] values Required 9 values: [s6,s8,s9,rs6,rs8,a1,a2,alpha,beta]; s9=0.
 * Unused fields have defaults; CSO stores a3/a4 in rs6/rs8 and does not use s8.
 */
void disprs_d3_get_named_parameters(disprs_d3_error error, disprs_d3_damping damping, char *method,
                                    double values[9]);
/** Create zero damping. Positive radii/exponents are caller constraints;
 * only scalar finiteness is checked.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] s6 Finite dimensionless C6 scale.
 * @param[in] s8 Finite dimensionless C8 scale.
 * @param[in] s9 Finite dimensionless ATM scale; |s9|<DBL_EPSILON disables ATM.
 * @param[in] rs6 Positive finite dimensionless C6 vdW-radius multiplier.
 * @param[in] rs8 Positive finite dimensionless C8 vdW-radius multiplier.
 * @param[in] alpha Positive finite C6 exponent; C8 and ATM use alpha+2.
 * @return Owned parameters, or NULL on error.
 */
disprs_d3_param disprs_d3_new_zero_damping(disprs_d3_error error, double s6, double s8, double s9,
                                           double rs6, double rs8, double alpha);
/** Create BJ damping. The radius a1*sqrt(C8/C6)+a2 must be positive
 * (caller constraint); only scalar finiteness is checked.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] s6 Finite dimensionless C6 scale.
 * @param[in] s8 Finite dimensionless C8 scale.
 * @param[in] s9 Finite dimensionless ATM scale; |s9|<DBL_EPSILON disables ATM.
 * @param[in] a1 Finite dimensionless damping-radius slope.
 * @param[in] a2 Finite damping-radius offset in bohr.
 * @param[in] alpha Finite ATM exponent control; exponent is alpha+2, normally alpha=14.
 * @return Owned parameters, or NULL on error.
 */
disprs_d3_param disprs_d3_new_rational_damping(disprs_d3_error error, double s6, double s8,
                                               double s9, double a1, double a2, double alpha);
/** Create modified-zero damping. Shifted power bases must stay positive;
 * only scalar finiteness is checked.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] s6 Finite dimensionless C6 scale.
 * @param[in] s8 Finite dimensionless C8 scale.
 * @param[in] s9 Finite dimensionless ATM scale; |s9|<DBL_EPSILON disables ATM.
 * @param[in] rs6 Positive finite dimensionless C6 vdW-radius multiplier.
 * @param[in] rs8 Positive finite dimensionless C8 vdW-radius multiplier.
 * @param[in] alpha Positive finite C6 exponent; C8 and ATM use alpha+2.
 * @param[in] beta Finite inverse-bohr shift in R/(rs*RvdW)+beta*RvdW.
 * @return Owned parameters, or NULL on error.
 */
disprs_d3_param disprs_d3_new_mzero_damping(disprs_d3_error error, double s6, double s8, double s9,
                                            double rs6, double rs8, double alpha, double beta);
/** Create modified BJ damping; fitted parameters distinguish it from BJ.
 * @copydetails disprs_d3_new_rational_damping
 */
disprs_d3_param disprs_d3_new_mrational_damping(disprs_d3_error error, double s6, double s8,
                                                double s9, double a1, double a2, double alpha);
/** Create optimized-power damping. The BJ radius a1*sqrt(C8/C6)+a2 must be
 * positive; only scalar finiteness is checked.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] s6 Finite dimensionless C6 scale.
 * @param[in] s8 Finite dimensionless C8 scale.
 * @param[in] s9 Finite dimensionless ATM scale; |s9|<DBL_EPSILON disables ATM.
 * @param[in] a1 Finite dimensionless damping-radius slope.
 * @param[in] a2 Finite damping-radius offset in bohr.
 * @param[in] alpha Finite ATM exponent control; exponent is alpha+2.
 * @param[in] beta Finite additional dimensionless power, normally nonnegative.
 * @return Owned parameters, or NULL on error.
 */
disprs_d3_param disprs_d3_new_optimizedpower_damping(disprs_d3_error error, double s6, double s8,
                                                     double s9, double a1, double a2, double alpha,
                                                     double beta);
/** Create CSO damping with no explicit C8 contribution. All scalars must be
 * finite (checked); the damping radius must be positive (caller constraint).
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] s6 Finite dimensionless long-range C6 scale.
 * @param[in] s9 Finite dimensionless ATM scale; |s9|<DBL_EPSILON disables ATM.
 * @param[in] a1 Finite sigmoid amplitude.
 * @param[in] a2 Finite multiplier of sqrt(C8/C6) at the sigmoid midpoint.
 * @param[in] a3 Finite dimensionless damping-radius slope.
 * @param[in] a4 Finite damping-radius offset in bohr.
 * @param[in] alpha Finite ATM exponent control; exponent is alpha+2.
 * @return Owned parameters, or NULL on error.
 */
disprs_d3_param disprs_d3_new_cso_damping(disprs_d3_error error, double s6, double s9, double a1,
                                          double a2, double a3, double a4, double alpha);
/** Create Z damping; all scalar inputs must be finite (checked).
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] s6 Finite dimensionless C6 scale.
 * @param[in] s8 Finite dimensionless C8 scale.
 * @param[in] s9 Finite dimensionless ATM scale; |s9|<DBL_EPSILON disables ATM.
 * @param[in] a1 Positive finite C6-dependent damping coefficient in atomic units.
 * Positivity is a caller constraint.
 * @param[in] alpha Finite ATM exponent control; exponent is alpha+2.
 * @return Owned parameters, or NULL on error.
 */
disprs_d3_param disprs_d3_new_z_damping(disprs_d3_error error, double s6, double s8, double s9,
                                        double a1, double alpha);
/** Free D3 damping parameters.
 * @param[in,out] param Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
void disprs_d3_delete_param(disprs_d3_param *param);

/** Compute CN and C6, including ghosts; results are not work-partitioned.
 * Does not compute dynamic polarizabilities. NULL outputs skip copying.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D3/D3S settings for N atoms; supplies the CN cutoff.
 * @param[out] coordination Optional N dimensionless coordination numbers.
 * @param[out] c6 Optional symmetric N*N C6 matrix in hartree bohr^6.
 */
void disprs_d3_get_properties(disprs_d3_error error, disprs_d3_structure structure,
                              disprs_d3_model model, double *coordination, double *c6);
/** Compute properties and full geometry responses, including ghosts.
 * Atomic/C6 Jacobians have Fortran shapes (3,N,N)/(3,N,N,N) and
 * (3,3,N)/(3,3,N,N); each property is a row with derivative index fastest.
 * Cartesian units are property units/bohr; strain retains property units.
 * Includes CN response; outputs are not work-partitioned. NULL outputs skip copying.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D3/D3S settings for N atoms; supplies the CN cutoff.
 * @param[out] coordination Optional N dimensionless coordination numbers.
 * @param[out] c6 Optional symmetric N*N C6 matrix in hartree bohr^6.
 * @param[out] coordination_cartesian Optional 3*N*N Cartesian CN derivatives.
 * @param[out] coordination_strain Optional 9*N strain CN derivatives.
 * @param[out] c6_cartesian Optional 3*N*N*N Cartesian C6 derivatives.
 * @param[out] c6_strain Optional 9*N*N strain C6 derivatives.
 */
void disprs_d3_get_property_response(disprs_d3_error error, disprs_d3_structure structure,
                                     disprs_d3_model model, double *coordination, double *c6,
                                     double *coordination_cartesian, double *coordination_strain,
                                     double *c6_cartesian, double *c6_strain);
/** Evaluate two-body plus enabled ATM terms with model cutoffs, ghosts, and partition.
 * Supports molecules and 1D/2D/3D real-space periodicity; Fourier restrictions
 * are documented at disprs_d3_set_model_ewald. Layouts follow the file contract.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D3/D3S settings for N atoms.
 * @param[in] param Live D3 damping/ATM parameters.
 * @param[out] energy Required scalar energy in hartree.
 * @param[out] gradient Optional 3*N xyz energy derivatives in hartree/bohr; NULL skips output.
 * @param[out] virial Optional 9 strain derivatives in hartree; NULL skips output.
 */
void disprs_d3_get_dispersion(disprs_d3_error error, disprs_d3_structure structure,
                              disprs_d3_model model, disprs_d3_param param, double *energy,
                              double *gradient, double *virial);
/** Decompose dispersion energy into symmetric pair matrices.
 * Summing all entries gives the partition's two-body/ATM energies; pair3 is
 * zero if ATM is off. Periodic diagonal entries can contain image contributions.
 * Supports real-space molecular/periodic evaluation; no periodic Fourier pairs.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D3/D3S settings for N atoms.
 * @param[in] param Live D3 damping/ATM parameters.
 * @param[out] pair2 Required N*N two-body energy matrix in hartree.
 * @param[out] pair3 Required N*N ATM energy matrix in hartree.
 */
void disprs_d3_get_pairwise_dispersion(disprs_d3_error error, disprs_d3_structure structure,
                                       disprs_d3_model model, disprs_d3_param param, double *pair2,
                                       double *pair3);
/** Analytical fixed-lattice Hessian including enabled ATM and CN response.
 * Uses model cutoffs/ghosts/partition; no periodic Fourier Hessians.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D3/D3S settings for N atoms.
 * @param[in] param Live D3 damping/ATM parameters.
 * @param[out] energy Required scalar energy in hartree.
 * @param[out] hessian Required column-major (3*N,3*N) matrix in hartree/bohr^2.
 */
void disprs_d3_get_dispersion_hessian(disprs_d3_error error, disprs_d3_structure structure,
                                      disprs_d3_model model, disprs_d3_param param, double *energy,
                                      double *hessian);

/** Load gCP parameters for ordered atomic numbers.
 * Names must resolve a supported set, e.g. method="pbeh3c", basis=NULL.
 * Defaults: basis/SRB cutoffs 60 bohr, serial work. No structure pointer is retained.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry whose elements must support the chosen set.
 * @param[in] method Optional read-only NUL-terminated UTF-8 method name.
 * @param[in] basis Optional read-only NUL-terminated UTF-8 basis name.
 * @return Owned gCP parameters, or NULL for unsupported elements/names.
 */
disprs_d3_gcp disprs_d3_load_gcp(disprs_d3_error error, disprs_d3_structure structure, char *method,
                                 char *basis);
/** Free gCP parameters.
 * @param[in,out] gcp Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
void disprs_d3_delete_gcp(disprs_d3_gcp *gcp);
/** Read gCP parameters; NULL outputs skip copying.
 * sigma scales gCP; alpha/beta define exp(-alpha*R^beta); dmp_scal/dmp_exp
 * define damping via dmp_scal*(R/RvdW)^dmp_exp. rscal scales the short-range
 * decay rate and qscal its energy amplitude; formulas use atomic units.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] gcp Live gCP parameters for N atoms.
 * @param[in] natoms Must equal the parameter atom count N.
 * @param[out] scalars Optional 7 values: sigma,alpha,beta,dmp_scal,dmp_exp,rscal,qscal.
 * @param[out] flags Optional 3 booleans: damp,base,srb (damping/base/SRB enabled).
 * @param[out] zeff Optional N effective atomic numbers in 1..36.
 * @param[out] emiss Optional N missing-basis energies in hartree.
 * @param[out] xv Optional N virtual-orbital counts.
 * @param[out] slater Optional N Slater exponents in bohr^-1.
 * @param[out] rvdw Optional symmetric N*N gCP pair radii in bohr, including diagonal.
 * @param[out] rvdw_srb Optional symmetric N*N SRB pair radii in bohr, including diagonal.
 */
void disprs_d3_get_gcp_parameters(disprs_d3_error error, disprs_d3_gcp gcp, int natoms,
                                  double *scalars, bool *flags, int *zeff, double *emiss,
                                  double *xv, double *slater, double *rvdw, double *rvdw_srb);
/** Replace selected gCP parameters; layouts/meaning follow disprs_d3_get_gcp_parameters.
 * Inputs are copied; the entire update is rejected if any field is invalid.
 * Direct slater/zeff overrides do not update eta; changing eta later regenerates slater.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] gcp Live gCP parameters for N atoms.
 * @param[in] natoms Must equal the parameter atom count N.
 * @param[in] scalars Optional 7 finite nonnegative scalar controls; NULL retains values.
 * @param[in] flags Optional 3 booleans: damp,base,srb; NULL retains values.
 * @param[in] zeff Optional N effective atomic numbers in 1..36; NULL retains values.
 * @param[in] emiss Optional N finite missing-basis energies in hartree; NULL retains values.
 * @param[in] xv Optional N finite virtual-orbital counts; NULL retains values.
 * @param[in] slater Optional N finite nonnegative Slater exponents in bohr^-1; NULL retains values.
 * @param[in] rvdw Optional N*N gCP radii in bohr; finite, positive, exactly symmetric.
 * NULL retains values.
 * @param[in] rvdw_srb Optional N*N SRB radii in bohr; finite, positive, exactly symmetric.
 * NULL retains values.
 */
void disprs_d3_set_gcp_parameters(disprs_d3_error error, disprs_d3_gcp gcp, int natoms,
                                  const double *scalars, const bool *flags, const int *zeff,
                                  const double *emiss, const double *xv, const double *slater,
                                  const double *rvdw, const double *rvdw_srb);
/** Read gCP controls; NULL outputs skip copying.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] gcp Live gCP parameters.
 * @param[out] eta Optional scalar Slater-exponent scale.
 * @param[out] base Optional base-correction enable flag.
 * @param[out] srb Optional short-range-basis correction enable flag.
 */
void disprs_d3_get_gcp_controls(disprs_d3_error error, disprs_d3_gcp gcp, double *eta, bool *base,
                                bool *srb);
/** Update gCP controls; invalid input preserves all controls.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] gcp Live gCP parameters to update.
 * @param[in] eta Optional finite nonnegative Slater scale; changing it regenerates
 * exponents, which must remain finite. NULL retains the value.
 * @param[in] base Optional base-correction enable flag; NULL retains the value.
 * @param[in] srb Optional short-range-basis correction enable flag; NULL retains the value.
 */
void disprs_d3_set_gcp_controls(disprs_d3_error error, disprs_d3_gcp gcp, const double *eta,
                                const bool *base, const bool *srb);
/** Set gCP hard-cutoff radii.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] gcp Live gCP parameters to update.
 * @param[in] basis Finite nonnegative basis-correction radius in bohr.
 * @param[in] short_range_basis Finite nonnegative SRB radius in bohr.
 */
void disprs_d3_set_gcp_realspace_cutoff(disprs_d3_error error, disprs_d3_gcp gcp, double basis,
                                        double short_range_basis);
/** Assign gCP terms to one worker.
 * No communication occurs; sum all workers' results. (0,1) restores serial work.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] gcp Live gCP parameters to update.
 * @param[in] part Zero-based worker index, 0<=part<parts.
 * @param[in] parts Positive worker count.
 */
void disprs_d3_set_gcp_work_partition(disprs_d3_error error, disprs_d3_gcp gcp, int part,
                                      int parts);
/** Compatibility stub: always reports unsupported MPI via error.
 * @param[in,out] error Receives unsupported-MPI status; NULL discards diagnostics.
 * @param[in] gcp Ignored; no state changes.
 * @param[in] communicator Ignored Fortran MPI communicator token.
 */
void disprs_d3_set_gcp_mpi_comm(disprs_d3_error error, disprs_d3_gcp gcp, int communicator);
/** Evaluate gCP/SRB for molecular or periodic geometry.
 * Uses gcp cutoffs/partition, independently of dispersion models and ghost masks.
 * Layouts follow the file contract.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms.
 * @param[in] gcp Live gCP parameters with identical ordered elements to structure.
 * @param[out] energy Required scalar energy in hartree.
 * @param[out] gradient Optional 3*N xyz energy derivatives in hartree/bohr; NULL skips output.
 * @param[out] virial Optional 9 strain derivatives in hartree; NULL skips output.
 */
void disprs_d3_get_counterpoise(disprs_d3_error error, disprs_d3_structure structure,
                                disprs_d3_gcp gcp, double *energy, double *gradient,
                                double *virial);
/** Analytical fixed-lattice gCP/SRB Hessian with gcp cutoffs/partition.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms.
 * @param[in] gcp Live gCP parameters with identical ordered elements to structure.
 * @param[out] energy Required scalar energy in hartree.
 * @param[out] hessian Required column-major (3*N,3*N) matrix in hartree/bohr^2.
 */
void disprs_d3_get_counterpoise_hessian(disprs_d3_error error, disprs_d3_structure structure,
                                        disprs_d3_gcp gcp, double *energy, double *hessian);

/** Return D4 compatibility version 40200 (4.2.0), not the native library version. */
int disprs_d4_get_version(void);
/** Allocate an initially clear error state; release with disprs_d4_delete_error. */
disprs_d4_error disprs_d4_new_error(void);
/** @copydoc disprs_check_error */
int disprs_d4_check_error(disprs_d4_error error);
/** @copydoc disprs_get_error */
void disprs_d4_get_error(disprs_d4_error error, char *buffer, const int *buffer_size);
/** @copydoc disprs_delete_error */
void disprs_d4_delete_error(disprs_d4_error *error);

/** Copy a geometry into an owned shared structure, checking D4 element support.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] natoms Atom count N>0, with representable array sizes.
 * @param[in] numbers Required N atomic numbers in 1..103 or 112..118.
 * @param[in] positions Required 3*N finite xyz values in bohr; separations >=1e-6 bohr.
 * @param[in] charge Optional finite total charge in electron-charge units; NULL means zero.
 * @param[in] lattice Optional 9 finite column-major cell values in bohr;
 * active vectors must be independent.
 * @param[in] periodic Optional 3 active-vector flags; NULL means all with a cell,
 * otherwise none. Active periodicity requires lattice.
 * @return Owned structure, or NULL on invalid input; all supplied data is copied.
 */
disprs_d4_structure disprs_d4_new_structure(disprs_d4_error error, int natoms, const int *numbers,
                                            const double *positions, const double *charge,
                                            const double *lattice, const bool *periodic);
/** @copydoc disprs_update_structure */
void disprs_d4_update_structure(disprs_d4_error error, disprs_d4_structure structure,
                                const double *positions, const double *lattice);
/** @copydoc disprs_delete_structure */
void disprs_d4_delete_structure(disprs_d4_structure *structure);

/** Create D4 or D4S settings.
 * Defaults: ga=3, gc=2, wf=6, EEQ, CN/pair/triplet cutoffs 30/60/40 bohr,
 * charge cutoff 60 bohr, no switching/ghosts/fixed charges, serial real-space work.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with D4-supported elements;
 * retains the atom count, not the structure pointer.
 * @param[in] model Interpolation selector DISPRS_D4 or DISPRS_D4S.
 * @return Owned model, or NULL on error.
 */
disprs_d4_model disprs_d4_new_model(disprs_d4_error error, disprs_d4_structure structure,
                                    disprs_d4_model_kind model);
/** Free D4/D4S model settings.
 * @param[in,out] model Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
void disprs_d4_delete_model(disprs_d4_model *model);
/** Add atoms to the ghost mask; duplicates are allowed.
 * Ghosts remain in CN/charge equilibration but not explicit dispersion terms.
 * Does not clear previous selections; invalid input preserves the mask.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings for N atoms.
 * @param[in] ghost Array of count zero-based indices in [0,N); NULL only if count=0.
 * @param[in] count Number of indices, >=0; zero leaves the mask unchanged.
 */
void disprs_d4_set_model_ghost_index(disprs_d4_error error, disprs_d4_model model, const int *ghost,
                                     int count);
/** Assign model terms to one worker.
 * Sum all workers' energy/derivative/pair outputs; properties remain complete.
 * No communication occurs. (0,1) restores serial work.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings to update.
 * @param[in] part Zero-based worker index, 0<=part<parts.
 * @param[in] parts Positive worker count.
 */
void disprs_d4_set_model_work_partition(disprs_d4_error error, disprs_d4_model model, int part,
                                        int parts);
/** Enable full-3D two-body Fourier dispersion.
 * ATM, Hessians, pair matrices and pair smoothing are unsupported.
 * No disable flag: create a new model to restore real-space summation.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings to update.
 * @param[in] rank Per-species-pair reference-C6 expansion cap; <=0 selects by tolerance.
 * @param[in] tolerance Finite factorization tolerance; <=0 means 1e-4.
 * @param[in] kcut Finite reciprocal radius in bohr^-1; <=0 is automatic.
 * @param[in] mesh Negative for direct summation, zero for automatic SPME, positive
 * for SPME rounded up to a power of two; storage must be representable.
 */
void disprs_d4_set_model_ewald(disprs_d4_error error, disprs_d4_model model, int rank,
                               double tolerance, double kcut, int mesh);
/** Create custom charge-scaling and CN settings; other defaults match disprs_d4_new_model.
 * All controls use atomic-unit conventions.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with D4-supported elements; not retained.
 * @param[in] model Interpolation selector DISPRS_D4 or DISPRS_D4S.
 * @param[in] ga Finite positive charge-scaling amplitude.
 * @param[in] gc Finite positive hardness multiplier in the charge-scaling exponent.
 * @param[in] wf Finite positive Gaussian CN width; D4S requires 6 and uses tabulated
 * pair widths instead.
 * @return Owned model, or NULL on invalid settings.
 */
disprs_d4_model disprs_d4_new_custom_model(disprs_d4_error error, disprs_d4_structure structure,
                                           disprs_d4_model_kind model, double ga, double gc,
                                           double wf);
/** Select charge equilibration and matching reference charges.
 * Fixed charges override the solver, not the reference-table selection.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings to update.
 * @param[in] charge_model DISPRS_EEQ or DISPRS_EEQBC; EEQBC supports Z=1..103.
 */
void disprs_d4_set_charge_model(disprs_d4_error error, disprs_d4_model model,
                                disprs_d4_charge_model charge_model);
/** Copy fixed charges without renormalization; derivatives are zero.
 * Invalid input preserves the previous setting.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings for N atoms.
 * @param[in] charges N finite geometry-independent charges in electron-charge units;
 * NULL with count=0 restores EEQ/EEQBC.
 * @param[in] count Must equal N for supplied charges, or zero with charges=NULL.
 */
void disprs_d4_set_fixed_charges(disprs_d4_error error, disprs_d4_model model,
                                 const double *charges, int count);
/** Set the periodic charge-summation cutoff.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings to update.
 * @param[in] cutoff Finite positive charge-summation radius in bohr.
 */
void disprs_d4_set_charge_cutoff(disprs_d4_error error, disprs_d4_model model, double cutoff);
/** Set real-space cutoffs. Width zero gives a hard cutoff; positive
 * widths give quintic switching and are limited to the corresponding radius.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in,out] model Live D4/D4S settings to update.
 * @param[in] cn Finite nonnegative CN radius in bohr.
 * @param[in] disp2 Finite nonnegative two-body radius in bohr.
 * @param[in] disp3 Finite nonnegative ATM radius in bohr.
 * @param[in] width2 Finite nonnegative two-body switching width in bohr.
 * @param[in] width3 Finite nonnegative ATM switching width in bohr.
 */
void disprs_d4_set_realspace_cutoff(disprs_d4_error error, disprs_d4_model model, double cn,
                                    double disp2, double disp3, double width2, double width3);

/** Load named BJ parameters, also usable with D4S.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] method Required read-only NUL-terminated UTF-8 functional name (before any '/').
 * @param[in] atm True selects the ATM fit; false prefers the two-body fit and sets s9=0.
 * @return Owned parameters, or NULL for unknown names.
 */
disprs_d4_param disprs_d4_load_param(disprs_d4_error error, char *method, bool atm);
/** Create BJ/ATM parameters. The radius a1*sqrt(C8/C6)+a2 and ATM exponent
 * must be positive (caller constraints); scalar finiteness is checked.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] s6 Finite dimensionless C6 scale.
 * @param[in] s8 Finite dimensionless C8 scale.
 * @param[in] s9 Finite dimensionless ATM scale; |s9|<DBL_EPSILON disables ATM.
 * @param[in] a1 Finite dimensionless damping-radius slope.
 * @param[in] a2 Finite damping-radius offset in bohr.
 * @param[in] alp Positive finite ATM exponent, normally 16; used directly, not D3's alpha+2.
 * @return Owned parameters, or NULL on error.
 */
disprs_d4_param disprs_d4_new_rational_damping(disprs_d4_error error, double s6, double s8,
                                               double s9, double a1, double a2, double alp);
/** Free D4 damping parameters.
 * @param[in,out] param Handle slot, set to NULL; NULL/already-cleared slots allowed.
 */
void disprs_d4_delete_param(disprs_d4_param *param);
/** Query the ATM-inclusive fit; unknown names leave output unchanged.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] method Required read-only NUL-terminated UTF-8 functional name.
 * @param[out] values Required 6 values: [s6,s8,s9,a1,a2,alp];
 * units/meaning follow disprs_d4_new_rational_damping.
 */
void disprs_d4_get_named_parameters(disprs_d4_error error, const char *method, double *values);
/** Query a fit with an optional ATM-scale override; output is preserved on error.
 * Returned s9 is the override; other parameters are fitted.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] method Required read-only NUL-terminated UTF-8 functional name.
 * @param[in] s9 Optional finite scalar override; NULL uses the fitted ATM scale.
 * Near-zero selects the two-body fit; otherwise selects the ATM fit.
 * @param[out] values Required 6 values: [s6,s8,s9,a1,a2,alp];
 * units/meaning follow disprs_d4_new_rational_damping.
 */
void disprs_d4_get_named_parameters_s9(disprs_d4_error error, const char *method, const double *s9,
                                       double *values);
/** Solve charges independently of any dispersion model settings.
 * Uses default charge cutoffs; omitting both Jacobians avoids derivative work.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms and total charge.
 * @param[in] charge_model DISPRS_EEQ or DISPRS_EEQBC; EEQBC requires Z<=103.
 * @param[out] charges Required N charges in electron-charge units.
 * @param[out] cartesian Optional 3*N*N property-major charge derivatives in charge/bohr;
 * NULL skips output.
 * @param[out] strain Optional 9*N property-major charge derivatives in charge units;
 * NULL skips output.
 */
void disprs_get_charges(disprs_d4_error error, disprs_d4_structure structure, int charge_model,
                        double *charges, double *cartesian, double *strain);

/** Compute properties; NULL outputs skip copying.
 * Includes ghost atoms; properties are not work-partitioned. D4S static
 * polarizabilities use same-element pair widths.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D4/D4S settings for N atoms.
 * @param[out] coordination_numbers Optional N dimensionless coordination numbers.
 * @param[out] partial_charges Optional N charges in electron-charge units.
 * @param[out] c6_coefficients Optional symmetric N*N C6 matrix in hartree bohr^6.
 * @param[out] polarizabilities Optional N static polarizabilities in bohr^3.
 */
void disprs_d4_get_properties(disprs_d4_error error, disprs_d4_structure structure,
                              disprs_d4_model model, double *coordination_numbers,
                              double *partial_charges, double *c6_coefficients,
                              double *polarizabilities);
/** Compute properties and their full geometry responses; NULL outputs skip copying.
 * Layout: property-major as in disprs_d3_get_property_response; units are those
 * of disprs_d4_get_properties divided by bohr for Cartesian derivatives.
 * Includes EEQ/EEQBC relaxation unless charges are fixed; no work partitioning.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D4/D4S settings for N atoms.
 * @param[out] coordination Optional N dimensionless coordination numbers.
 * @param[out] charges Optional N charges in electron-charge units.
 * @param[out] c6 Optional symmetric N*N C6 matrix in hartree bohr^6.
 * @param[out] polarizabilities Optional N static polarizabilities in bohr^3.
 * @param[out] coordination_cartesian Optional 3*N*N Cartesian CN derivatives.
 * @param[out] coordination_strain Optional 9*N strain CN derivatives.
 * @param[out] charge_cartesian Optional 3*N*N Cartesian charge derivatives.
 * @param[out] charge_strain Optional 9*N strain charge derivatives.
 * @param[out] c6_cartesian Optional 3*N*N*N Cartesian C6 derivatives.
 * @param[out] c6_strain Optional 9*N*N strain C6 derivatives.
 * @param[out] polarizability_cartesian Optional 3*N*N Cartesian polarizability derivatives.
 * @param[out] polarizability_strain Optional 9*N strain polarizability derivatives.
 */
void disprs_d4_get_property_response(disprs_d4_error error, disprs_d4_structure structure,
                                     disprs_d4_model model, double *coordination, double *charges,
                                     double *c6, double *polarizabilities,
                                     double *coordination_cartesian, double *coordination_strain,
                                     double *charge_cartesian, double *charge_strain,
                                     double *c6_cartesian, double *c6_strain,
                                     double *polarizability_cartesian,
                                     double *polarizability_strain);
/** Evaluate two-body plus enabled ATM terms.
 * Uses model cutoffs, ghosts, charges, and partition; charge response is included
 * unless charges are fixed. Supports molecular/1D/2D/3D real space; Fourier
 * restrictions follow disprs_d4_set_model_ewald. Layouts follow the file contract.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D4/D4S settings for N atoms.
 * @param[in] param Live D4 BJ/ATM parameters.
 * @param[out] energy Required scalar energy in hartree.
 * @param[out] gradient Optional 3*N xyz energy derivatives in hartree/bohr; NULL skips output.
 * @param[out] virial Optional 9 strain derivatives in hartree; NULL skips output.
 */
void disprs_d4_get_dispersion(disprs_d4_error error, disprs_d4_structure structure,
                              disprs_d4_model model, disprs_d4_param param, double *energy,
                              double *gradient, double *virial);
/** Decompose dispersion energy into symmetric pair matrices.
 * Their sums give the partition's two-body/ATM energies; disabled ATM yields zeros.
 * Periodic diagonals may contain image contributions; Fourier mode is unsupported.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D4/D4S settings for N atoms.
 * @param[in] param Live D4 BJ/ATM parameters.
 * @param[out] pair2 Required N*N two-body energy matrix in hartree.
 * @param[out] pair3 Required N*N ATM energy matrix in hartree.
 */
void disprs_d4_get_pairwise_dispersion(disprs_d4_error error, disprs_d4_structure structure,
                                       disprs_d4_model model, disprs_d4_param param, double *pair2,
                                       double *pair3);
/** Central-difference Hessian at fixed lattice, step=1e-4 bohr.
 * Geometry is preserved. Includes charge relaxation unless fixed; no Fourier mode.
 * Displacements must not cross cutoff/image boundaries. Prefer the analytical routine.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D4/D4S settings for N atoms.
 * @param[in] param Live D4 BJ/ATM parameters.
 * @param[out] hessian Required column-major (3*N,3*N) matrix in hartree/bohr^2.
 */
void disprs_d4_get_numerical_hessian(disprs_d4_error error, disprs_d4_structure structure,
                                     disprs_d4_model model, disprs_d4_param param, double *hessian);
/** Analytical fixed-lattice Hessian.
 * Includes enabled ATM, CN response, and charge relaxation unless charges are fixed.
 * Uses model cutoffs/ghosts/partition; supports molecular/periodic real space,
 * not Fourier. Output is preserved on error.
 * @param[in,out] error Receives status; NULL discards diagnostics.
 * @param[in] structure Live shared geometry with N atoms supported by model.
 * @param[in] model Live D4/D4S settings for N atoms.
 * @param[in] param Live D4 BJ/ATM parameters.
 * @param[out] hessian Required column-major (3*N,3*N) matrix in hartree/bohr^2.
 */
void disprs_d4_get_dispersion_hessian(disprs_d4_error error, disprs_d4_structure structure,
                                      disprs_d4_model model, disprs_d4_param param,
                                      double *hessian);

#ifdef __cplusplus
}
#endif
