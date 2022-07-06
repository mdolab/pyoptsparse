# Much of this is from SciPy

project(
  'pyoptsparse',
  'c', 'cpp',
# unnecessary metadata commented out until Meson supports PEP517 and installation with pip
#  version: 'x.x.x',
#  license: 'GPL-3',
  meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'cpp_std=c++14',
  ],
)

fortranobject_c = '../fortranobject.c'

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)

# https://mesonbuild.com/Python-module.html
# Here we differentiate from the python used by meson, py3_command, and that python target, py3_target. This is useful
# when cross compiling like on conda-forge
py_mod = import('python')
py3_command = py_mod.find_installation()
if get_option('python_target') != ''
    py3_target = py_mod.find_installation(get_option('python_target'))
else
    py3_target = py3_command
endif
py3_dep = py3_target.dependency()

subdir('pyoptsparse')