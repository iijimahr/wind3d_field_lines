project(
  'wind3d_field_lines',
  ['c', 'fortran'],
  version: '0.2.0',
  meson_version: '>=1.5.0',
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

numpy_include_dir = run_command(
  py,
  [
    '-c',
    'import numpy; print(numpy.get_include())',
  ],
  check: true,
).stdout().strip()

f2py_include_dir = run_command(
  py,
  [
    '-c',
    'import numpy.f2py; print(numpy.f2py.get_include())',
  ],
  check: true,
).stdout().strip()

fortranobject_c = run_command(
  py,
  [
    '-c',
    'import numpy.f2py, pathlib; print(pathlib.Path(numpy.f2py.get_include()) / "fortranobject.c")',
  ],
  check: true,
).stdout().strip()

numpy_includes = include_directories(numpy_include_dir, f2py_include_dir)

subdir('src/wind3d_field_lines')
