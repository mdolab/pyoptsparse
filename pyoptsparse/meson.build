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

# Meson prefixes filenames of intermediate compiled objects with their filepath. This poses a problem for conda builds
# since conda ensures the host environment directory has 255 characters so the meson object filenames then exceed 255
# characters. To remedy this, the fortranobject.c file from numpy is copied into pyoptsparse so that the meson build
# uses a relative path, rather than an absolute path, thus reducing the auto generated object filename
run_command(py3,
  [
    '-c',
    'import os; os.chdir(".."); import shutil; import numpy; inc = numpy.get_include(); shutil.copy(os.path.join(inc, "..", "..", "f2py", "src", "fortranobject.c"), "pyoptsparse")'
  ],
  check: true
)

inc_np = include_directories(incdir_numpy)


# TODO: pyoptsparse supports numpy>=1.16 but numpy.f2py.get_include() wasnt added until later, raise numpy version?
#incdir_f2py = run_command(py3,
#  [
#    '-c',
#    'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'
#  ],
#  check : true
#).stdout().strip()
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)

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

# test imports
envdata = environment()
python_paths = [join_paths(meson.current_build_dir(), '..')]
envdata.prepend('PYTHONPATH', python_paths)

progs = [['SLSQP',  'pySLSQP',  'slsqp'],
         ['CONMIN', 'pyCONMIN', 'conmin'],
         ['PSQP',   'pyPSQP',   'psqp'],
         ['NSGA2',  'pyNSGA2',  'nsga2']]


foreach p : progs
  import_command = 'from pyoptsparse.' + p[1] + ' import '+p[2]+'; print('+p[2]+'.__file__)'
  test(
    'import test for '+p[0],
    py3,
    args: ['-c', import_command],
    env: envdata
  )
endforeach
