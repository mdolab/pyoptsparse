project(
  'pyoptsparse',
  'c', 'cpp',
  # Note that the git commit hash cannot be added dynamically here (it is added
  # in the dynamically generated and installed `scipy/version.py` though - see
  # tools/version_utils.py
  version: '2.8.2',
  license: 'GPL-3',
  meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'cpp_std=c++14',
    'blas=openblas',
    'lapack=openblas'
  ],
)

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
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()

subdir('pyoptsparse')