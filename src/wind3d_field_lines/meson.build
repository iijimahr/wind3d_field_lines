py.install_sources(
  [
    '__init__.py',
    'demo_arcade.py',
    'demo_bipolar.py',
    'integrator.py',
    'potential_field.py',
    'types.py',
  ],
  subdir: 'wind3d_field_lines',
)

f2py_wrappers = custom_target(
  'bbtobln_wrappers',
  input: 'fortran/field_line_integrator.f90',
  output: [
    '_bbtoblnmodule.c',
    '_bbtobln-f2pywrappers.f',
  ],
  command: [
    py,
    '-m',
    'numpy.f2py',
    '@INPUT@',
    '-m',
    '_bbtobln',
    '--lower',
    '--build-dir',
    '@OUTDIR@',
  ],
)

py.extension_module(
  '_bbtobln',
  [
    'fortran/field_line_integrator.f90',
    fortranobject_c,
    f2py_wrappers,
  ],
  include_directories: numpy_includes,
  dependencies: py_dep,
  subdir: 'wind3d_field_lines',
  install: true,
)
