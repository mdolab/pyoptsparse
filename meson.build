project(
  'pyoptsparse',
  'c', 'cpp',
  meson_version: '>= 0.64',
  default_options: [
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
    'c_std=c17',
    'cpp_std=c++17',
  ],
)

# <!-- from scipy meson.build

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.get_id() == 'gcc'
  # -std=legacy is not supported by all Fortran compilers, but very useful with
  # gfortran since it avoids a ton of warnings that we don't care about.
  # Needs fixing in Meson, see https://github.com/mesonbuild/meson/issues/11633.
  add_project_arguments('-std=legacy', language: 'fortran')
endif

if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

if host_machine.system() == 'darwin'
  if cc.has_link_argument('-Wl,-dead_strip')
    # Allow linker to strip unused symbols
    add_project_link_arguments('-Wl,-dead_strip', language : ['c', 'cpp', 'fortran'])
  endif
endif

# --!>

# install python sources
install_subdir('pyoptsparse', install_dir: py3.get_install_dir())

# install non-python sources
subdir('pyoptsparse')
