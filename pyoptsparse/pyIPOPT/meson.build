fs = import('fs')

if get_option('ipopt_dir') != '' or fs.is_dir('Ipopt')

    if is_windows
        ipopt_dirs = []
        ipopt_idir = ''
        ipopt_dirs = [get_option('ipopt_dir') / 'Library'/ 'lib']
        ipopt_idir = get_option('ipopt_dir') / 'include' / 'coin-or'

        ipopt_dep = cc.find_library('ipopt-3', required: false)
    else
        ipopt_dirs = []
        ipopt_idir = ''

        # Ipopt installs differently on some systems (i.e. Fedora)
        if fs.is_dir(get_option('ipopt_dir') / 'lib')
            ipopt_dirs = [get_option('ipopt_dir') / 'lib']
        elif fs.is_dir(get_option('ipopt_dir') / 'lib64')
            ipopt_dirs = [get_option('ipopt_dir') / 'lib64']
        endif

        if fs.is_dir(get_option('ipopt_dir') / 'include' / 'coin-or')
            ipopt_idir = get_option('ipopt_dir') / 'include' / 'coin-or'
        elif fs.is_dir(get_option('ipopt_dir') / 'include' / 'coin')
            ipopt_idir = get_option('ipopt_dir') / 'include' / 'coin'
        elif not fs.is_dir('Ipopt')
            error('IPOPT include directory not found in: ', get_option('ipopt_dir') / 'include')
        endif

        # We also include the special case for Azure
        if fs.is_dir('Ipopt')
            ipopt_dirs = ['Ipopt' / 'lib']
            ipopt_idir = 'Ipopt' / 'include' / 'coin-or'
        endif
    endif

    ipopt_dep = cc.find_library('ipopt', required: true, dirs: ipopt_dirs)
    if fs.is_dir(ipopt_idir)
        ipopt_inc = include_directories(ipopt_idir)
    else
        error('IPOPT include directory not found: ', ipopt_inc)
    endif

    py3.extension_module('pyipoptcore',
                    'src/callback.c',
                    'src/pyipoptcoremodule.c',
                    include_directories: [inc_np, 'src', ipopt_inc],
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