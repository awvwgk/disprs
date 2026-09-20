testdrive_dep = dependency('test-drive', fallback: ['test-drive', 'testdrive_dep'])
mstore_dep = dependency('mstore', fallback: ['mstore', 'mstore_dep'])
reference_python = import('python').find_installation()
if get_option('d3_reference') != ''
d3_c_api = executable('upstream-d3-c-api',
  get_option('d3_reference') / 'test/api/api-test.c',
  dependencies: disprs_dep, c_args: ['-UNDEBUG'])
test('upstream-d3-c-api', d3_c_api, suite: 'reference')
d3_tests = []
foreach source : ['test_dftd3.f90', 'test_param.f90', 'test_gcp.f90', 'test_pairwise.f90', 'test_regression.f90', 'test_hessian.f90', 'test_gcp_hessian.f90', 'test_partition.f90', 'test_fourier.f90', 'test_periodic_1d.f90', 'test_periodic_2d.f90', 'test_periodic_3d.f90', 'test_periodic_atm.f90']
d3_tests += configure_file(
  input: get_option('d3_reference') / 'test/unit' / source,
  output: source,
  command: [reference_python, files('d3_tests.py'), '@INPUT@', '@OUTPUT@'],
)
endforeach
reference_d3 = executable(
  'reference-d3',
  ['main.f90', d3_tests],
  link_with: fortran_api,
  dependencies: [disprs_dep, mctc_dep, mstore_dep, testdrive_dep],
)
test('upstream-d3-molecular-and-parameters', reference_d3, suite: 'reference', timeout: 600)
test('upstream-d3-regression', reference_d3, args: ['regression'], suite: 'reference', timeout: 120)
test('upstream-d3-hessians', reference_d3, args: ['hessians'], suite: 'reference', timeout: 300)
test('upstream-d3-partitions', reference_d3, args: ['partitions'], suite: 'reference', timeout: 120)
test('upstream-d3-fourier', reference_d3, args: ['fourier'], suite: 'reference', timeout: 300)
endif

if get_option('d4_reference') != ''
  upstream_d4_dep = dependency('dftd4', version: '4.2.0', required: false)
  upstream_charge_dep = dependency('multicharge', version: '0.5.0', required: false)
  upstream_d3_dep = dependency('s-dftd3', version: '1.6.0', required: false)
  if upstream_d4_dep.found() and upstream_charge_dep.found() and upstream_d3_dep.found()
    upstream_fortran = executable('upstream-fortran-api', 'fortran_api.f90',
      dependencies: [upstream_d4_dep, upstream_charge_dep, upstream_d3_dep, mctc_dep],
      fortran_args: ['-cpp'])
    fortran_values = custom_target('upstream-fortran-api-values',
      output: 'fortran-api-values.txt', command: [upstream_fortran], capture: true)
    native_fortran = executable('native-fortran-api', 'fortran_api.f90',
      link_with: fortran_api, dependencies: [disprs_dep, mctc_dep],
      fortran_args: ['-cpp', '-DDISPRS_NATIVE'])
    test('upstream-fortran-api', native_fortran, args: [fortran_values], suite: 'reference', timeout: 120)
  endif
  d4_c_api = executable('upstream-d4-c-api',
    get_option('d4_reference') / 'test/api/example.c',
    dependencies: disprs_dep, c_args: ['-UNDEBUG'])
  test('upstream-d4-c-api', d4_c_api, suite: 'reference')
  exporter = executable('reference-geometries', 'geometries.f90', dependencies: [mctc_dep, mstore_dep])
  geometries = custom_target(
    'reference-geometries-data',
    output: 'geometries.txt',
    command: [exporter],
    capture: true,
  )
  cases = custom_target(
    'reference-d4-cases',
    input: [get_option('d4_reference') / 'test/unit/test_dftd4.f90', geometries,
      get_option('d4_reference') / 'test/unit/test_param.f90',
      get_option('d4_reference') / 'test/unit/test_pairwise.f90',
      get_option('d4_reference') / 'test/unit/test_periodic.f90'],
    output: 'd4-cases.txt',
    command: [reference_python, files('d4_cases.py'), '@INPUT0@', '@INPUT1@', '@INPUT2@', '@INPUT3@', '@INPUT4@'],
    capture: true,
  )
  reference_d4 = executable('reference-d4', 'd4.rs',
    dependencies: rust_dependencies, rust_args: ['--edition=2021'])
  test('upstream-d4-molecular-and-parameters', reference_d4, args: [cases], suite: 'reference', timeout: 120)
  if get_option('d3_reference') != ''
    model_cases = custom_target('upstream-model-cases',
      input: [geometries, get_option('d3_reference') / 'test/unit/test_model.f90',
        get_option('d4_reference') / 'test/unit/test_model.f90'],
      output: 'model-cases.txt',
      command: [reference_python, files('model_cases.py'), '@INPUT@'], capture: true)
    reference_models = executable('reference-models', 'models.rs',
      dependencies: rust_dependencies, rust_args: ['--edition=2021'])
    test('upstream-models', reference_models, args: [model_cases], suite: 'reference')
  endif
endif