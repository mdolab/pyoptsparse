# Platform detection
is_windows = host_machine.system() == 'windows'
is_msvc = is_windows and cc.get_id() == 'msvc'

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
