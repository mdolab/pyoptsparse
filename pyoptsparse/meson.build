# NumPy include directory - needed in all submodules
incdir_numpy = get_option('incdir_numpy')
if incdir_numpy == ''
    incdir_numpy = run_command(py3_target,
      [
        '-c',
        'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
      ],
      check: true
    ).stdout().strip()
endif
# this creates a raw string which is useful for Windows use of '\' for paths
incdir_numpy = '''@0@'''.format(incdir_numpy)

# HACK: Meson prefixes filenames of intermediate compiled objects with their filepath. This poses a problem for conda builds
# since conda ensures the host environment directory has 255 characters so the meson object filenames then exceed 255
# characters. To remedy this, the fortranobject.c file from numpy is copied into pyoptsparse so that the meson build
# uses a relative path, rather than an absolute path, thus reducing the auto generated object filename
# see for example https://github.com/mesonbuild/meson/issues/4226
run_command(py3_command,
  [
    '-c',
    'import os; os.chdir(".."); import shutil; shutil.copy(os.path.join(r"' + incdir_numpy + '", "..", "..", "f2py", "src", "fortranobject.c"), "pyoptsparse")'
  ],
  check: true
)

inc_np = include_directories(incdir_numpy)


# TODO: pyoptsparse supports numpy>=1.16 but numpy.f2py.get_include() wasnt added until later, raise numpy version?
#incdir_f2py = run_command(py3_target,
#  [
#    '-c',
#    'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'
#  ],
#  check : true
#).stdout().strip()
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)


# TODO: this is kept in here so that when meson becomes pep517-compliant
# we can uncomment these to have meson install the source files

#python_sources = [
#    '__init__.py',
#    'pyOpt_MPI.py',
#    'pyOpt_constraint.py',
#    'pyOpt_error.py',
#    'pyOpt_gradient.py',
#    'pyOpt_history.py',
#    'pyOpt_objective.py',
#    'pyOpt_optimization.py',
#    'pyOpt_optimizer.py',
#    'pyOpt_solution.py',
#    'pyOpt_utils.py',
#    'pyOpt_variable.py',
#    'types.py'
#]

#py3_target.install_sources(
#  python_sources,
#  pure: true,
#  subdir: 'pyoptsparse'
#)

subdir('pySNOPT')
subdir('pySLSQP')
subdir('pyCONMIN')
subdir('pyNLPQLP')
subdir('pyNSGA2')
subdir('pyPSQP')
#subdir('pyALPSO')
#subdir('pyParOpt')
#subdir('postprocessing')

# test imports
# envdata = environment()
# python_paths = [join_paths(meson.current_build_dir(), '..')]
# envdata.prepend('PYTHONPATH', python_paths)

# progs = [['SLSQP',  'pySLSQP',  'slsqp'],
#          ['CONMIN', 'pyCONMIN', 'conmin'],
#          ['PSQP',   'pyPSQP',   'psqp'],
#          ['NSGA2',  'pyNSGA2',  'nsga2']]


# foreach p : progs
#   import_command = 'from pyoptsparse.' + p[1] + ' import '+p[2]+'; print('+p[2]+'.__file__)'
#   test(
#     'import test for '+p[0],
#     py3_command,
#     args: ['-c', import_command],
#     env: envdata
#   )
# endforeach
