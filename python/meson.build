project(
  'disprs-python',
  version: '0.1.0',
  license: 'Apache-2.0',
  meson_version: '>=1.3',
)

disprs_dep = dependency(
  'disprs',
  method: 'pkg-config',
  version: ['>=0.1.0', '<0.2.0'],
)
python = import('python').find_installation('python3', required: true)
python.install_sources(
  'disprs/__init__.py',
  'disprs/_framework.py',
  'disprs/ase_ext.py',
  'disprs/pyscf_ext.py',
  'disprs/qcschema_ext.py',
  subdir: 'disprs',
)

test_env = environment()
test_env.set('PYTHONPATH', meson.current_source_dir())
test_env.set('OMP_NUM_THREADS', '1')
test_env.set('OPENBLAS_NUM_THREADS', '1')
test('Python API smoke test', python, args: files('test.py'), env: test_env)
test('Python framework interfaces', python, args: [files('test_frameworks.py'), '-v'], env: test_env, timeout: 120)