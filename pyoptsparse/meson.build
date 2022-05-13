# Platform detection
is_windows = host_machine.system() == 'windows'
is_mingw = is_windows and cc.get_id() == 'gcc'
is_msvc = is_windows and cc.get_id() == 'msvc'

cython_c_args = []
if is_windows
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  if is_mingw
    add_global_link_arguments(gcc_link_args, language: ['c', 'cpp'])
    # Force gcc to float64 long doubles for compatibility with MSVC
    # builds, for C only.
    add_global_arguments('-mlong-double-64', language: 'c')
    # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
    add_global_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c', 'cpp'])
    # Manual add of MS_WIN64 macro when not using MSVC.
    # https://bugs.python.org/issue28267
    bitness = run_command('_build_utils/gcc_build_bitness.py').stdout().strip()
    if bitness == '64'
      add_global_arguments('-DMS_WIN64', language: ['c', 'cpp', 'fortran'])
    endif
    # Silence warnings emitted by PyOS_snprintf for (%zd), see
    # https://github.com/rgommers/scipy/issues/118.
    # Use as c_args for extensions containing Cython code
    cython_c_args += ['-Wno-format-extra-args', '-Wno-format']
  endif
  if is_msvc
    # Keeps msvc from complaining about database file access
    add_global_arguments('/MP', language: 'c')
  endif
  if meson.get_compiler('fortran').get_id() == 'gcc'
    add_global_link_arguments(gcc_link_args, language: ['fortran'])
    # Flag needed to work around BLAS and LAPACK Gfortran dependence on
    # undocumented C feature when passing single character string
    # arguments.
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
    # https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
    add_global_arguments('-fno-optimize-sibling-calls',
      language: ['fortran'])
  endif
endif

blas = dependency(get_option('blas'))
lapack = dependency(get_option('lapack'))

# NumPy include directory - needed in all submodules
incdir_numpy = run_command(py3,
  [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
  ],
  check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
# numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

python_sources = [
    '__init__.py',
    'pyOpt_MPI.py',
    'pyOpt_constraint.py',
    'pyOpt_error.py',
    'pyOpt_gradient.py',
    'pyOpt_history.py',
    'pyOpt_objective.py',
    'pyOpt_optimization.py',
    'pyOpt_optimizer.py',
    'pyOpt_solution.py',
    'pyOpt_utils.py',
    'pyOpt_variable.py',
    'types.py'
]

py3.install_sources(
  python_sources,
  pure: true,
  subdir: 'pyoptsparse'
)

subdir('pySNOPT')
subdir('pyIPOPT')
subdir('pySLSQP')
subdir('pyCONMIN')
subdir('pyNLPQLP')
subdir('pyNSGA2')
subdir('pyPSQP')
subdir('pyALPSO')
subdir('pyParOpt')
subdir('postprocessing')
