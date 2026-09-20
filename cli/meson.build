project(
  'disprs-cli',
  'rust',
  version: '0.1.0',
  meson_version: '>=1.12',
  default_options: ['buildtype=release'],
)

rust = import('rust')
package = rust.workspace().package()
disprs_dep = dependency('disprs', method: 'pkg-config', version: ['>=0.1.0', '<0.2.0'])
library_name = host_machine.system() == 'darwin' ? 'libdisprs.dylib' : 'libdisprs.so'
library_path = disprs_dep.get_variable(pkgconfig: 'libdir') / library_name
config = configuration_data()
config.set_quoted('VERSION', meson.project_version())
config.set_quoted('LIBRARY', library_path)
build_config = configure_file(
  input: 'src/build_config.rs.in',
  output: 'build_config.rs',
  configuration: config
)
cli = package.executable(
  structured_sources(
    [
      'src/main.rs',
      'src/disprs_ext.rs',
      'src/qcschema.rs',
      build_config,
    ]
  ),
  install: true,
  rust_args: ['--cfg', 'meson'],
)

if not meson.is_cross_build()
  subdir('assets')
  test_env = environment()
  test_env.set('DISPRS_LIBRARY', library_path)
  test_env.set('OPENBLAS_NUM_THREADS', '1')
  test_env.set('OMP_NUM_THREADS', '1')
  rust.test('CLI', cli, env: test_env, timeout: 180)
  test('CLI installed-library lookup', cli, args: ['parameters', 'pbe', '--name', 'd3(op)'])
endif