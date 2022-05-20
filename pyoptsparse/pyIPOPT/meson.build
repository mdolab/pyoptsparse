fs = import('fs')

if get_option('ipopt_dir') != '' or fs.is_dir('Ipopt')

    ipopt_dir = ''

    if get_option('ipopt_dir') != ''
        ipopt_dir = get_option('ipopt_dir')
    elif fs.is_dir('Ipopt')
        ipopt_dir = fs.is_dir('Ipopt')
    endif

    ipopt_lib = []
    ipopt_idir = ''

    # Ipopt installs differently on some systems (i.e. Fedora)
    if fs.is_dir(ipopt_dir / 'lib')
        ipopt_lib = [ipopt_dir / 'lib']
    elif fs.is_dir(ipopt_dir / 'lib64')
        ipopt_lib = [ipopt_dir / 'lib64']
    endif


    if fs.is_dir(ipopt_dir / 'include' / 'coin-or')
        ipopt_idir = ipopt_dir / 'include' / 'coin-or'
    elif fs.is_dir(ipopt_dir / 'include' / 'coin')
        ipopt_idir = ipopt_dir / 'include' / 'coin'
    endif

    ipopt_dep = cc.find_library('ipopt-3', required: false, dirs: ipopt_lib)   # only relevant on windows
    if not ipopt_dep.found()
        ipopt_dep = cc.find_library('ipopt', required: true, dirs: ipopt_lib)
    endif

    if fs.is_dir(ipopt_idir)
        ipopt_inc = include_directories(ipopt_idir)
    else
        error('IPOPT include directory not found: ', ipopt_dir)
    endif

    py3_target.extension_module('pyipoptcore',
                    'src/callback.c',
                    'src/pyipoptcoremodule.c',
                    include_directories: [inc_np, 'src', ipopt_inc],
                    dependencies : [py3_dep, ipopt_dep],
                    subdir: 'pyoptsparse/pyIPOPT',
                    link_language: 'c',
                    install : false)
endif

#python_sources = [
#    '__init__.py',
#    'pyIPOPT.py',
#]
#
#py3_target.install_sources(
#  python_sources,
#  pure: false,
#  subdir: 'pyoptsparse/pyIPOPT'
#)
