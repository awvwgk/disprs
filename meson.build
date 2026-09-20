project(
  'disprs',
  ['rust', 'c'],
  version: '0.1.0',
  license: ['Apache-2.0'],
  meson_version: '>=1.11',
  default_options: [
    'buildtype=release',
    'c_std=c11',
    'b_ndebug=true',
    'warning_level=3',
    'default_library=shared',
  ],
)
python = import('python').find_installation(required: get_option('python'), pure: false)
python_package = get_option('python_package')
if python_package and (not python.found() or not get_option('fortran').disabled())
  error('python_package requires Python enabled and Fortran disabled')
endif
library_install_dir = get_option('libdir')
if (get_option('d3_reference') != '' or get_option('d4_reference') != '') and not get_option('fortran').enabled()
  error('Reference tests require -Dfortran=enabled')
endif
if python_package
  library_install_dir = python.get_install_dir() / 'disprs'
endif
rust_dependencies = import('rust').workspace().package().dependencies()
libdisprs = library(
  'disprs',
  'src/lib.rs',
  dependencies: rust_dependencies,
  rust_args: ['--edition=2021'],
  rust_crate_type: 'cdylib',
  install: true,
  install_dir: library_install_dir,
)
disprs_dep = declare_dependency(
  link_with: libdisprs,
  include_directories: include_directories('include'),
)

pkg = import('pkgconfig')
if not python_package
  install_headers('include/disprs.h')
  if get_option('d3_compat')
    install_headers('include/dftd3.h')
  endif
  if get_option('d4_compat')
    install_headers('include/dftd4.h')
  endif
  pkg.generate(
    libdisprs,
    name: 'disprs',
    description: 'D3, gCP, and D4 dispersion corrections',
  )
  install_data('LICENSE', install_dir: get_option('datadir') / 'licenses' / 'disprs')
endif

c_smoke = executable(
  'c-smoke',
  'tests/c_smoke.c',
  dependencies: disprs_dep,
  c_args: ['-UNDEBUG'],
)
test('C API smoke test', c_smoke)

c_compat = executable(
  'c-compat',
  'tests/c_compat.c',
  dependencies: disprs_dep,
  c_args: ['-UNDEBUG'],
)
test('D3/D4 C compatibility', c_compat, timeout: 120)

if python.found()
  python.install_sources(
    'python/disprs/__init__.py',
    'python/disprs/_framework.py',
    'python/disprs/ase_ext.py',
    'python/disprs/pyscf_ext.py',
    'python/disprs/qcschema_ext.py',
    subdir: 'disprs',
  )
  python_test_env = environment()
  python_test_env.set('DISPRS_LIBRARY', libdisprs.full_path())
  python_test_env.set('PYTHONPATH', meson.current_source_dir() / 'python')
  python_test_env.set('OMP_NUM_THREADS', '1')
  python_test_env.set('OPENBLAS_NUM_THREADS', '1')
  test(
    'Python API smoke test',
    python,
    args: files('python/test.py'),
    env: python_test_env,
  )
  test(
    'Python framework interfaces',
    python,
    args: [files('python/test_frameworks.py'), '-v'],
    env: python_test_env,
    timeout: 120,
  )
endif

if add_languages('fortran', required: get_option('fortran'))
  fc = meson.get_compiler('fortran')
  install_python = find_program('python3')
  mctc_dep = dependency(
    'mctc-lib',
    fallback: ['mctc-lib', 'mctc_dep'],
    default_options: ['json=disabled', 'openmp=false'],
  )
  fortran_api = library(
    'disprs_fortran',
    ['fortran/disprs_ext.f90', 'fortran/dftd3_ext.f90', 'fortran/dftd4_ext.f90', 'fortran/multicharge_ext.f90'],
    dependencies: [disprs_dep, mctc_dep],
    install: true,
  )
  fortran_module_id = 'disprs-fortran' / fc.get_id() + '-' + fc.version()
  meson.add_install_script(
    install_python,
    files('fortran/install-mod.py'),
    get_option('includedir') / fortran_module_id,
  )
  pkg.generate(
    fortran_api,
    filebase: 'disprs-fortran',
    name: 'disprs-fortran',
    description: 'Fortran API for disprs',
    requires: ['disprs', 'mctc-lib'],
    subdirs: ['', fortran_module_id],
  )
  fortran_smoke = executable(
    'fortran-smoke',
    'fortran/test.f90',
    link_with: fortran_api,
    dependencies: [disprs_dep, mctc_dep],
  )
  test('Fortran API smoke test', fortran_smoke)
  fortran_compat = executable(
    'fortran-compat',
    'fortran/test_dftd3.f90',
    link_with: fortran_api,
    dependencies: [disprs_dep, mctc_dep],
  )
  test('s-dftd3 Fortran compatibility', fortran_compat)
  foreach failure : ['model', 'damping', 'gcp-constructor', 'gcp-cutoff', 'gcp-hessian']
    test('s-dftd3 Fortran invalid ' + failure, fortran_compat,
      args: ['--invalid-' + failure], should_fail: true)
  endforeach
  fortran_d4 = executable('fortran-d4', 'fortran/test_dftd4.f90',
    link_with: fortran_api, dependencies: [disprs_dep, mctc_dep])
  test('dftd4 Fortran compatibility', fortran_d4, timeout: 120)
  fortran_charge = executable('fortran-charge', 'fortran/test_multicharge.f90',
    link_with: fortran_api, dependencies: [disprs_dep, mctc_dep])
  test('multicharge Fortran compatibility', fortran_charge, timeout: 120)
  if get_option('d3_reference') != '' or get_option('d4_reference') != ''
    subdir('tests/reference')
  endif
endif
