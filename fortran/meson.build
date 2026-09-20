project(
  'disprs-fortran',
  'fortran',
  version: '0.1.0',
  license: 'Apache-2.0',
  meson_version: '>=1.3',
  default_options: ['default_library=shared'],
)

disprs_dep = dependency(
  'disprs',
  method: 'pkg-config',
  version: ['>=0.1.0', '<0.2.0'],
)
mctc_dep = dependency(
  'mctc-lib',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['json=disabled', 'openmp=false'],
)
fc = meson.get_compiler('fortran')

fortran_api = library(
  'disprs_fortran',
  ['disprs_ext.f90', 'dftd3_ext.f90', 'dftd4_ext.f90', 'multicharge_ext.f90'],
  dependencies: [disprs_dep, mctc_dep],
  install: true,
)

module_id = meson.project_name() / fc.get_id() + '-' + fc.version()
meson.add_install_script(
  find_program('python3'),
  files('install-mod.py'),
  get_option('includedir') / module_id,
)

pkg = import('pkgconfig')
pkg.generate(
  fortran_api,
  filebase: meson.project_name(),
  name: meson.project_name(),
  description: 'Fortran API for disprs',
  requires: ['disprs', 'mctc-lib'],
  subdirs: ['', module_id],
)

fortran_smoke = executable(
  'fortran-smoke',
  'test.f90',
  link_with: fortran_api,
  dependencies: [disprs_dep, mctc_dep],
)
test('Fortran API smoke test', fortran_smoke)

fortran_compat = executable(
  'fortran-compat',
  'test_dftd3.f90',
  link_with: fortran_api,
  dependencies: [disprs_dep, mctc_dep],
)
test('s-dftd3 Fortran compatibility', fortran_compat)
foreach failure : ['model', 'damping', 'gcp-constructor', 'gcp-cutoff', 'gcp-hessian']
  test('s-dftd3 Fortran invalid ' + failure, fortran_compat,
    args: ['--invalid-' + failure], should_fail: true)
endforeach
fortran_d4 = executable('fortran-d4', 'test_dftd4.f90',
  link_with: fortran_api, dependencies: [disprs_dep, mctc_dep])
test('dftd4 Fortran compatibility', fortran_d4, timeout: 120)
fortran_charge = executable('fortran-charge', 'test_multicharge.f90',
  link_with: fortran_api, dependencies: [disprs_dep, mctc_dep])
test('multicharge Fortran compatibility', fortran_charge, timeout: 120)