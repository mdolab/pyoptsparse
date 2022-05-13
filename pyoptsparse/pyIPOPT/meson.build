# This is done like this so conda-build succeeds
py_path = run_command(py3,
  [
    'find_python_path.py'
  ],
  check: true
).stdout().strip()

if is_windows
    ipopt_dep = cc.find_library('ipopt-3', required: false)
else
    ipopt_dep = dependency('ipopt', required: false)
endif

if ipopt_dep.found()

    if is_windows
        coinor_inc = include_directories(py_path / 'Library' / 'include' / 'coin-or')
    else
        coinor_inc = include_directories(py_path / '..' / 'include' / 'coin-or')
    endif

    py3.extension_module('pyipoptcore',
                    'src/callback.c',
                    'src/pyipoptcoremodule.c',
                    include_directories: [inc_np, 'src', coinor_inc],
                    dependencies : [py3_dep, ipopt_dep],
                    subdir: 'pyoptsparse/pyIPOPT',
                    link_language: 'c',
                    install : true)
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