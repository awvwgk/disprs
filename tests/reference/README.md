# Upstream reference coverage

This is partial upstream coverage, not certification that disprs passes both
complete upstream suites. Unsupported tests are not counted as passing.

## D3S Reference

`cargo test --lib d3::smooth` includes all 41 D3S(BJ)/BLYP energies from
[the authors' implementation](https://github.com/ntkachenko95/D3S/tree/a09648fb5c2042c6b06df0346d3c6b71807b7d40).
`d3s_scan.xyz` and `d3s_scan_energies.txt` are unmodified copies of
`Test_geometries/H2_NiKur_BLYP_1D_scan.xyz` and `D3S_BJ.txt` at that revision.
The reference geometries use angstrom; the test uses the authors' conversion
of 0.52917726 angstrom/bohr. The library retains its existing CODATA-derived
covalent radii, and the published SI pair widths are rounded to ten decimal
places. Energy tolerance is 1e-7 hartree; existing upstream tolerances are
unchanged. Expected energies are never regenerated from the native code.

Additional native tests check both CN derivative orders, Cartesian forces and
Hessians, C6 responses, strain/virial derivatives, ghost exclusion, work
partition sums, invalid elements, and molecular/periodic cases. Zero, BJ,
modified-zero/BJ, optimized-power, CSO and Z damping reuse the existing kernels;
ATM tests exercise an optional extension beyond the authors' two-body code.

## NN-Dx PyTorch Oracle

The optional [Rust comparison](../nndx_reference.rs) executes
[nndx.py](nndx.py), which imports the unmodified author fork at
`599e840148b8f018fceaed585be3e87a300e4724`. It checks the checkout revision and
tracked-file cleanliness before generating any values. This is a physics-backend
comparison with deterministic test residuals and the fork's example isolated-atom
spectra, not a reproduction of a trained NN-Dx functional or its benchmarks.

Use Python 3.12: the fork declares Python `<3.13` and NumPy `<2`.

```sh
git clone https://github.com/Yufan121/tad-dftd4.git /tmp/disprs-nndx-reference
git -C /tmp/disprs-nndx-reference checkout 599e840148b8f018fceaed585be3e87a300e4724
python3.12 -m venv /tmp/disprs-nndx-py312
/tmp/disprs-nndx-py312/bin/python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
/tmp/disprs-nndx-py312/bin/python -m pip install /tmp/disprs-nndx-reference \
  'tad-mctc==0.5.3' 'tad-multicharge==0.3.3' 'qcelemental==0.51.2' 'pydantic>=2,<3'
DISPRS_NNDX_PYTHON=/tmp/disprs-nndx-py312/bin/python \
DISPRS_NNDX_REFERENCE=/tmp/disprs-nndx-reference \
cargo test --test nndx_reference -- --ignored --nocapture
```

The reference test is explicitly ignored in ordinary Cargo runs, not counted
as passing when dependencies are absent. All reference values are produced live
by the fork and PyTorch autograd; no runtime reference dependency or stored
native-generated expected values are used. It checks separate two-/three-body
energies, C6/C8, scaled spectra, coordinate/charge derivatives, all 30 atomic
input derivatives, and native EEQ charges and response. Cases include neutral
and charged molecules, zero/nonzero shifts, disabled ATM, zero/tiny spectra,
the charge floor, heavy-element tables, and hard-cutoff boundaries. The pinned
fork's two-edge ordered ATM mask is tested literally, including its fractional
contribution for triangles with exactly two in-range edges.

The fork does not declare its `pydantic` dependency, and its `rcov` coordination
argument is incompatible with current `tad-mctc`; the compatible charge/toolkit
versions above are enforced by the oracle. Float64 is selected before importing
these packages so their module-level EEQ tables are not rounded to float32.
The NN-Dx EEQ path matches their CN floor and radius conversion; it does not
change standard D4's Fortran conventions. Validated with Python 3.12.14,
PyTorch 2.14.0+cpu, NumPy 1.26.4 and Pydantic 2.13.5: all 35 oracle cases pass.

Energy tolerance is `100*epsilon`; charge tolerance is `1e-12`; coefficient and
spectrum tolerance is `1e-11 + 1e-12*abs(reference)`; derivative tolerance is
`1e-12 + 1e-10*abs(reference)`. Native unit tests independently check finite
differences for every atomic input and EEQ-composed coordinate derivatives, plus
input validation. No D3/D4/gCP reference tolerances are modified.

## Reproduce

```sh
meson setup build-reference -Dfortran=enabled \
  -Dd3_reference=/path/to/simple-dftd3 \
  -Dd4_reference=/path/to/dftd4
meson test -C build-reference --suite reference --print-errorlogs
```

Either source path may be omitted. test-drive 0.5.0 and mstore 0.3.0 are fetched
through commit-pinned wraps when not installed. Production never links upstream
dispersion or charge libraries. The optional Fortran API oracle described below
links them in a separate reference executable only.

Tested source revisions:

- simple-dftd3 1.6.0: `6004cfd62ad2a248ad32167392bd08c674da6161`.
- dftd4 4.2.0: `6e1f59c3f39d919a2dbef0601d2576727c8b30e8`.

CI checks out these exact public release snapshots. A newer suite needs a fresh
inventory, not an assumption of coverage.

## Executed cases

- D3 pairwise: all 12 upstream cases, including crystals, ATM and cutoff
  smoothing. The adapter translates constant damping records into equivalent
  owning-handle constructors; inputs and the `100*epsilon` energy/pair-sum
  assertions are unchanged. These are consistency checks, not independent
  pair-matrix reference values.
- D3 regression: the original 90-atom Ta `dftbplus#871` energy/gradient-norm
  regression passes unchanged at `100*epsilon`.
- D3/gCP Hessians: all 10 dispersion and all 9 gCP Hessian cases, including
  periodic gCP, retain original finite-difference steps, symmetry checks and
  tolerances. Only test/error and gCP module imports are adapted.
- D3 partitions: 5 of 9 original cases (invalid, serial, dispersion, Hessian,
  counterpoise) retain their numerical assertions, including all six damping
  families and atomic energies. Library errors are translated to test-drive;
  expected local errors retain their library type. CN-object, on-demand C6,
  reducer and unsupported Fourier+ATM cases are explicitly excluded.
- D3 Fourier: 9 of 32 original cases cover fixed/automatic SPME energy and
  derivatives, BJ/zero real-space convergence, reciprocal-cutoff convergence,
  supercell/translation invariance, Cartesian and strain derivatives. These
  tests only use atomic energies through their sum, so the adapter calls the
  scalar API instead; all comparisons, tolerances and steps are retained.
- D3 model: 7 of 9 original cases cover four literal reference-weight arrays
  (including AmF3), two CN-weight derivative checks and all-ghost energy/force/
  virial norms. The two symbol-versus-number table-overload checks are excluded
  object-API tests.
- D4/D4S model: 18 of 23 original cases cover literal EEQ/EEQBC reference
  weights, CN/charge weight derivatives, polarizabilities and their derivatives,
  and unsupported-element rejection. Three GFN2-reference-charge cases remain
  unsupported; two Fortran polymorphic-factory cases are excluded object-API
  tests. No EEQ substitution is made for GFN2.

The model runner includes the production Rust sources in a test-only executable,
without exporting new public model APIs. Its extractor requires the exact
registered case counts and reports every exclusion. Original tolerances remain
`100*epsilon` for D3 weights/ghost norms and `sqrt(epsilon)` for derivatives and
D4 weights/polarizabilities; finite differences retain the `1e-6` step. Weight
derivative checks are consistency tests, while stored upstream weight and
polarizability arrays provide independent reference values. The combined model
runner requires both source options; the D3 Fortran suites need only D3.

- D3: all 44 cases in `test/unit/test_dftd3.f90`, compiled as Fortran against
  disprs. Only `use mctc_env_testing` is replaced with `use testdrive` in the
  build directory. Assertions, inputs, expected energies, and tolerances remain
  unchanged. Includes molecular energy, Cartesian/strain finite differences,
  ATM, all seven damping families, actinides, and cutoff smoothing.
- D3 parameters: all 21 groups in `test/unit/test_param.f90`, including named
  parameter energy arrays, explicit Z parameters and seven expected-error
  groups. The build adapter translates mctc library errors into test-drive
  errors without altering inputs, energy assertions or tolerances.
- D3 periodic: the 1D, 2D, 3D and periodic-ATM suites are retained, with
  original inputs, energies and tolerances. Ten of fourteen original cases
  pass. Four energy cases (`gh185`, acetic, adaman and cyanamide) are explicit
  expected failures due to upstream image-box truncation, not counted as parity.
  Three additional complete-image energy oracles from upstream's own kernels
  pass at the original `100*epsilon` tolerance. All original Cartesian and
  strain derivative assertions pass without changed tolerances or steps.
- gCP: all 66 active cases in `test/unit/test_gcp.f90`, retained as Fortran
  with test-drive. Module imports are adapted and the unused formatting import
  removed; inputs, allocation checks, mutations and numerical assertions are
  unchanged. The two cases commented out upstream are not counted.
- D4: all 41 D4/D4S cases (34 MB16-43 including EEQBC, two AmF3, two actinide, two custom-model,
  one smooth-cutoff) from
  `test/unit/test_dftd4.f90`, run in Rust against the native implementation.
  Parameters and reference energies are extracted from the source; geometries
  are exported by Fortran/mstore or extracted from explicit upstream arrays.
  Energy tolerance remains `100*epsilon` and
  derivative tolerance `sqrt(epsilon)`, with the original `1e-6` step.
  All three checks in each AmF3 case and smoothing pairwise/gradient/strain
  checks are retained, giving 47 executable records. The extractor requires
  exactly 41 cases.
- D4 parameters: 118 named-energy checks and 67 libxc alias pairs from
  `test/unit/test_param.f90`. The named-energy loop retains upstream's
  15-bohr ATM cutoff and `100*epsilon` tolerance; its unused trailing reference
  value is not treated as an additional test.
- D4 pairwise: all six cases in `test/unit/test_pairwise.f90`, including two
  crystals, retain the original total-energy versus pair-sum assertion at
  `100*epsilon`. These are consistency checks, not independent pair-matrix values.
- D4 periodic: all ten cases in `test/unit/test_periodic.f90` are retained.
  Six Cartesian/strain derivative cases pass with the original `1e-6`/`1e-7`
  steps and `sqrt(epsilon)`/`100*sqrt(epsilon)` tolerances. Four original energy
  assertions are explicit expected failures due to the image defects below,
  not counted as parity. Four independent corrected-image energy checks pass
  at `100*epsilon`. Overall the D4 runner reports 177 passing original records,
  four expected failures, four corrected-image checks and 67 alias pairs.
- C APIs: both complete pinned programs, simple-dftd3 `test/api/api-test.c` and
  dftd4 `test/api/example.c`, compile and run unchanged against the installed-style
  compatibility headers and native disprs library. No upstream dispersion library
  is linked. Assertions are explicitly enabled even in release builds.
  The local `tests/c_compat.c` additionally checks all seven D3 constructors,
  all named loaders (including the expected missing Z entry), gCP optional outputs,
  and all four D4 model constructors across molecular/1D/2D/3D geometry. D4
  Hessians, gradients and virials are checked by finite differences; properties,
  charge conservation, pair sums, exact state preservation, buffer bounds, invalid
  handles/inputs, deletion and error recovery are also checked.

Source files and fixtures remain in their upstream projects with their original
license notices. Generated test inputs are build artifacts, not replacement
reference values obtained from disprs.

### Fortran API Oracle

With `d4_reference` enabled and installed pkg-config dependencies `dftd4=4.2.0`,
`s-dftd3=1.6.0` and `multicharge=0.5.0`, `fortran_api.f90` is compiled twice: once
against upstream and once against disprs. The upstream executable generates
build-local reference values; the native consumer compares eight combinations of
D4/D4S, EEQ/EEQBC and molecular/3D boundaries, including custom models, smooth
cutoffs, properties, pair matrices, native analytical versus upstream numerical Hessians and analytical charge
derivatives. Energy and pair comparisons use `100*epsilon`; all other tolerances
are explicit in the consumer. Missing reference libraries omit this optional
test; the local Fortran tests always run when Fortran is enabled.

The same consumer compares D3 atomic energies, pair matrices, gradients, virials
and Hessians for 16 combinations of BJ/zero damping, ATM on/off and molecular,
1D, 2D and 3D boundaries. It checks all 118 D4 integer parameter IDs plus invalid
IDs, with default, zero, epsilon and fractional ATM scaling. All six damping
values and allocation status are compared, including the separate DFTB two-body fits.
It also compares properties, energies and gradients for D4/D4S with each element
112-118, plus standalone EEQ charges. Local Rust/Python tests retain invalid-element
and EEQBC-limit checks. Python tests additionally cover analytical-Hessian partition
sums and independent Cartesian/strain charge responses across all periodicities.

D3 property checks compare CN and undamped C6 for molecular and skew 1D/2D/3D
cells at three CN cutoffs, updating coordinates and the cell between queries.
The `DISPRS_NATIVE` preprocessing branch calls the disprs property extension;
the upstream branch directly calls its CN, reference-weight and atomic-C6
methods. Upstream has no corresponding top-level property query. Comparisons
retain `1e-12` CN and `1e-10` C6 tolerances; no upstream runtime fallback is used.
Local binding tests additionally cover optional outputs, symmetry, ghost/partition
invariance, shape errors, failed-output preservation and recovery.

Sixteen fixed-charge D4/D4S cases compare C6, polarizabilities, energies,
gradients and virials for EEQ/EEQBC reference data and molecular/1D/2D/3D cells.
The upstream branch supplies charges directly to its model and damping kernels,
retaining neutral ATM and CN derivatives without charge relaxation. Its explicit
image-repetition overload uses three repetitions on active axes and zero on
inactive axes, sufficient for this fixed geometry and all physical cutoffs.
The cutoff-based upstream helper ignores directional flags once any is true;
using explicit repetitions avoids that defect without changing tolerances or
existing upstream expectations. No upstream code is linked in production.

The self-consistent D4/charge 1D/2D cases use local finite differences and conservation checks instead of
claiming parity with upstream's different finite-image EEQ treatment. Neither
this high-level consumer nor the existing D3 tests establish complete public
model/type-bound API coverage, which is explicitly out of scope. See the Fortran compatibility section in the
top-level README for the remaining interface gaps.

## Remaining coverage and blockers

| Upstream suite | Current status |
| --- | --- |
| D3 molecular dispersion | 44/44 executed |
| D3 model | 7/9 numerical cases ported and passing; 2 symbol-overload API checks excluded |
| D3 parameters | 21/21 groups executed; optional citation outputs still absent. Upstream 1.6.0 also has no named-Z entries |
| D3 pairwise | All 12 upstream cases retained and passing, plus 16 independent Fortran atomic/pairwise/derivative comparisons |
| D3 periodic | 10 original cases plus 3 independent complete-image oracles pass; 4 legacy image-box energy assertions are explicit expected failures. Native 1D/2D/3D derivatives, Hessians, all damping families, ATM and translations tested |
| D3 gCP | 66/66 active cases retained; parameter arrays, base/SRB switches and eta exposed |
| D3 regression | Original Ta crystal energy/gradient regression retained and passing |
| D3/gCP Hessians | All 10 D3 and 9 gCP upstream cases retained and passing |
| D3 partitions | 5/9 retained and passing; CN-object, on-demand C6, reducer and Fourier+ATM cases excluded |
| D3 Fourier | 9/32 ported and passing; remaining low-level eigensolver/kernel, low-rank object/atomic-output, citation, cutoff-estimator and rejection cases are not claimed as upstream coverage |
| D4 Fourier | Local direct/SPME extension, not an upstream Fourier-suite port: converged 600-bohr native real-space comparison, mesh refinement, analytical response differences, charge variants, custom models, ghosts/partitions and binding checks; existing upstream real-space oracles remain unchanged |
| D4 molecular dispersion | 41/41 executed, including both EEQBC cases |
| D4 model | 18/23 numerical/error cases ported and passing; 3 GFN2-reference-charge cases unsupported, 2 object-factory cases excluded |
| D4 parameters | 118 named energies, 67 alias pairs and all 118 Fortran integer IDs with four ATM settings; low-level types remain incomplete |
| D4 pairwise | All 6 upstream cases retained and passing, including 2 periodic cases |
| D4 periodic | All 10 upstream cases retained: 6 derivative cases pass, 4 original energies are explicit expected failures; 4 independent corrected-image oracles pass. Local directional EEQ/EEQBC and partial EEQ Ewald checks also retained |
| C API programs | Both pinned upstream C programs retained unchanged and passing; extensive local C compatibility checks supplement them |
| Fortran export examples | Existing disprs tests; not a port of every upstream example |
| Python framework/interface tests | Existing disprs tests only, not complete upstream Python suites |
| D3 CLI validation | Not executable against disprs: no upstream-compatible CLI |
| D3 MPI/output | MPI backend and upstream citation/formatting APIs are not implemented; suites excluded |

Passing the complete suites requires implementing the missing functionality and
porting the remaining checks, not merely changing the test harness. Do not relax
tolerances, replace expected values, or mark unsupported cases as successes.

The local Rust test `periodic_eeqbc_matches_dftd4` retains independently
generated dftd4 4.2.0 / multicharge 0.5.0 values for a charged C/O skew cell,
for both D4 and D4S: energy, charges, C6, polarizabilities, gradients and virial.
It is not counted as a ported upstream periodic suite. Local finite differences
also check charged EEQBC responses and custom-model Cartesian/lattice derivatives.

For periodic EEQBC, compatibility includes a multicharge 0.5.0 self-image
indexing defect: `get_pairs` compacts distances after excluding the origin,
then indexes the uncompressed translation list. The EEQBC kernel reproduces
that predecessor-image selection; changing it changes periodic reference
charges. Revisit this behavior when updating the pinned multicharge version.
The existing EEQ image selection is unchanged.

Directional EEQBC filters inactive lattice directions and does not reproduce
the 3D-only predecessor-image defect. Partial-periodic EEQ uses Gaussian
slab/wire Ewald sums, not a vacuum-padded 3D Ewald calculation. Kernel tests
check neutral direct-image convergence, splitting invariance, self terms,
special-function values, far-field stability and coordinate/lattice/width
derivatives. The Python test changes the real-space splitting hint and
requires invariant charges. See [the derivation](../../EWALD.md).

## Defects found

The retained C API programs exposed unchecked D3 null handles, missing rejection
of coincident structures, stale errors after successful calls, a null error handle
incorrectly reporting success, and missing default 512-byte error-buffer semantics.
These are fixed at the shared FFI/D3 entry points. Structure updates validate before
mutation; both models reject nonfinite geometry and D3 rejects nonfinite explicit
parameters. gCP now honors gradient-only and virial-only requests independently.
Both structure constructors now default to full periodicity when a lattice is
supplied without periodic flags, matching the upstream constructor contract.
The D4 numerical-Hessian API uses native gradients and upstream's exact
step/layout; no upstream fallback is involved.
Fallible C constructors, configuration, parameter queries and numerical evaluators
contain unwinding panics and propagate native errors. Exported setup calls have
injected-panic checks for null returns, preserved state/outputs and recovery.
The C compatibility suite checks gCP control overflow/recovery and nonfinite D3/D4
energy, derivative, pairwise and Hessian results without changing output buffers.
Rust checks also cover unrepresentable atom counts, derived geometry and FFT mesh
overflow, and malformed gCP names. These checks do not cover invalid/wrong-type
pointers, undersized buffers, allocation aborts or `panic=abort` builds.

The periodic suites exposed missing scalar-PBC broadcasting and missing
lattice updates in the compatible Fortran D3 adapter; both are fixed.
The `gh185` reference differs by 1.83812e-9 hartree because upstream's
origin-centered CN image box omits images for its unwrapped skew-cell geometry.
A diagnostic reproducing that box restored the original value, but would
make results depend on which cell contains an atom. Production instead uses
complete pair-centered image bounds and tests translation by thirteen cells.
The original `gh185` assertion remains intact and is marked as an expected
failure; the adapter adds an early error return to avoid overwriting it.
No energy reference or tolerance is replaced, and this is not reported as
full upstream periodic parity.

Three X23 energy fixtures have the same image-completeness issue: acetic
(8.32124e-10 hartree), adaman (4.90565e-9) and cyanamide (1.58433e-8).
The standalone `d3_images.f90` oracle evaluates these structures with upstream
1.6.0's own CN, C6 and damping kernels, generating translations out to 120 bohr
while keeping the actual CN and dispersion cutoffs at 30 and 60 bohr. Its
complete sums agree with native results at the original strict tolerance.
`d3_tests.py` keeps the original assertions as expected failures and adds
separate cases for these independently generated values. The oracle is not
linked into production or the retained-suite executables. To regenerate with
installed upstream s-dftd3 and mstore development packages:

```sh
gfortran $(pkg-config --cflags s-dftd3 mstore) tests/reference/d3_images.f90 \
  $(pkg-config --libs s-dftd3 mstore) -o /tmp/disprs-d3-images
OMP_NUM_THREADS=1 /tmp/disprs-d3-images
```

Tight periodic ATM differences exposed accumulation noise from the much larger
real-space traversal. Compensated energy summation fixes all three original
periodic-ATM derivative checks without weakening their tolerances.

The D4 acetic/adaman energy cases also omit dispersion images. Acetic additionally
exposes two multicharge 0.5.0 issues: the self-image indexing defect described
above, and a nearest-image search restricted to 27 cells applied to unwrapped
coordinates (one atom lies more than two cells along the short lattice axis).
Native EEQ already wraps coordinates and selects actual nearest self-images;
no production behavior was changed to reproduce these defects.

The independent `d4_images.f90` oracle uses upstream D4/D4S and multicharge
kernels. It prints four energies per case: original, complete dispersion images,
complete CN plus dispersion images, and additionally corrected EEQ images.
It generates translations to 120 bohr but retains the original physical cutoffs.
For the last column only, it wraps a copy of the charge geometry and repairs
the upstream cache's self-image indices before upstream matrix assembly and
LAPACK solution. Corrected acetic matrices agree with native within `5e-16`,
RHS within `1.4e-15` and charges within `2.2e-14`. All four final energies agree
at the original `100*epsilon` tolerance. Those last-column values are separate
regression checks in `d4.rs`; the original extracted values remain unchanged,
and unexpected passes require revisiting the defect classification.

To regenerate with upstream development packages (never linked into production
or the retained-suite runner):

```sh
gfortran -ffree-line-length-none -J/tmp \
  $(pkg-config --cflags dftd4 mctc-lib multicharge mstore) \
  tests/reference/d4_images.f90 \
  $(pkg-config --libs dftd4 mctc-lib multicharge mstore) -o /tmp/disprs-d4-images
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /tmp/disprs-d4-images
```

With Meson's fallback mstore, replace its pkg-config entry with
`-Ibuild-reference/subprojects/mstore/libmstore.so.0.3.0.p` and
`build-reference/subprojects/mstore/libmstore.a` before the upstream libraries.

The retained D3 suite exposed missing `+2` in the ATM exponent for zero,
modified-zero, optimized-power, CSO and Z constructors, also affecting named loading.
It also exposed rounded Bohr conversion in the D3 data generator. Regenerated
radii now use the same CODATA-derived conversion as mctc-lib; the unchanged
upstream suite passes its original strict tolerance.

The parameter suite additionally exposed missing aliases, missing BOP,
revPBE0, PBE38, M06-HF and HCTH120 zero-damping records, and OPBE zero-damping
`s8=2.033` instead of upstream `2.055`. These are fixed in the shared lookup
and embedded table, so every binding receives the same corrections.

The retained gCP suite exposed stale mSVP Li/Be basis counts (10 instead of
the pinned 1.6.0 values of 9). Regenerating the gCP asset from the pinned source
fixes both mSVP and PBEh-3c reference cases without changing tolerances.