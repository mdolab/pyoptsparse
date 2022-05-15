# This is done like this so conda-build succeeds
py_path = run_command(py3,
  [
    'find_python_path.py'
  ],
  check: true
).stdout().strip()

# The following first tries
if is_windows
    ipopt_dirs = []
    ipopt_idir = ''
    if get_option('ipopt_dir') != ''
        ipopt_dirs = [get_option('ipopt_dir') / 'Library'/ 'lib']
        ipopt_idir = get_option('ipopt_dir') / 'include' / 'coin-or'
    endif
    ipopt_dep = cc.find_library('ipopt-3', required: false)
else
    ipopt_dirs = []
    ipopt_idir = ''
    if get_option('ipopt_dir') != ''
        ipopt_dirs = [get_option('ipopt_dir') / 'lib']
        ipopt_idir = get_option('ipopt_dir') / 'include' / 'coin-or'
    endif
    ipopt_dep = cc.find_library('ipopt', required: false, dirs: ipopt_dirs)
endif

if ipopt_dep.found()
    fs = import('fs')

    if fs.is_dir(ipopt_idir)

        ipopt_inc = include_directories(ipopt_idir)
        py3.extension_module('pyipoptcore',
                        'src/callback.c',
                        'src/pyipoptcoremodule.c',
                        include_directories: [inc_np, 'src', ipopt_inc],
                        dependencies : [py3_dep, ipopt_dep],
                        subdir: 'pyoptsparse/pyIPOPT',
                        link_language: 'c',
                        install : true)
    endif
endif

python_sources = [
    '__init__.py',
    'pyIPOPT.py',
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'pyoptsparse/pyIPOPT'
)