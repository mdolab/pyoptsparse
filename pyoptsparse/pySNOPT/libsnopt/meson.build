project(
  'snopt-module',
  'c',
  meson_version: '>= 1.1',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'fortran_std=f2008',
  ],
)
sourcedir = get_option('sourcedir')
pyos_sourcedir = '..' / 'source'

fs = import('fs')
if sourcedir == ''
  assert(false, 'sourcedir must be specified')
elif fs.is_file(sourcedir)
  is_dyn_lib=true
  lib_name=fs.name(sourcedir)
  if lib_name.startswith('lib')
    lib_name = lib_name.substring(3)
  endif
  lib_dir=fs.parent(sourcedir)
elif fs.is_dir(sourcedir)
  add_languages('fortran')
  is_dyn_lib=false
else
  assert(false, 'sourcedir must be either a path to a file or directory')
endif
# Get Python
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()

# Get NumPy include directory
incdir_numpy = run_command(py3,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true
).stdout().strip()
inc_np = include_directories(incdir_numpy)

# Get f2py include directory
incdir_f2py = run_command(py3,
  ['-c', 'import os; import numpy; print(os.path.join(numpy.get_include(), "..", "..", "f2py", "src"))'],
  check: true
).stdout().strip()
inc_f2py = include_directories(incdir_f2py)

# Copy fortranobject.c locally
run_command(py3,
  ['-c', 'import os, shutil, numpy; shutil.copy(os.path.join(numpy.get_include(), "..", "..", "f2py", "src", "fortranobject.c"), ".")'],
  check: true
)
fortranobject_c = 'fortranobject.c'

if is_dyn_lib
  # Find precompiled SNOPT library
  cc = meson.get_compiler('c')
  snopt_lib = cc.find_library(lib_name, dirs: ['{lib_dir}'], required: true)
else
  # Get list of all Fortran files
  result = run_command(pyos_sourcedir / 'grab-all-fortran-files.py', sourcedir, check: true)
  snopt_source_files = result.stdout().strip().splitlines()
  snopt_source_files += [pyos_sourcedir / 'closeunit.f', pyos_sourcedir / 'openunit.f']
endif

# Generate C wrapper from f2py interface
snopt_source = custom_target('snoptmodule.c',
  input: [pyos_sourcedir /'f2py' / 'snopt.pyf'],
  output: 'snoptmodule.c',
  command: [py3, '-m', 'numpy.f2py', '@INPUT@', '--lower', '--build-dir', '.']
)

# Build extension module
deps = [py3_dep]
sources = [snopt_source, fortranobject_c]
if is_dyn_lib
  deps += [snopt_lib]
else
  sources += [snopt_source_files]
endif

py3.extension_module('snopt',
  sources,
  include_directories: [inc_np, inc_f2py],
  dependencies: deps,
  install: true,
)
