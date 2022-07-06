conmin_source = custom_target('conminmodule.c',
                            input : ['source/f2py/conmin.pyf',
                                    ],
                            output : ['conminmodule.c', 'conmin-f2pywrappers.f'],
                            command: [py3_command, '-m', 'numpy.f2py', '@INPUT@',
                            '--lower', '--build-dir', 'pyoptsparse/pyCONMIN']
                            )

py3_target.extension_module('conmin',
                    'source/openunit.f',
                    'source/cnmn00.f',
                    'source/cnmn01.f',
                    'source/cnmn02.f',
                    'source/cnmn03.f',
                    'source/cnmn04.f',
                    'source/cnmn05.f',
                    'source/cnmn06.f',
                    'source/cnmn07.f',
                    'source/cnmn08.f',
                    'source/cnmn09.f',
                    'source/conmin.f',
                    'source/closeunit.f',
                conmin_source,
                fortranobject_c,
                include_directories: [inc_np, inc_f2py],
                dependencies : py3_dep,
                subdir: 'pyoptsparse/pyCONMIN',
                install : false,
                build_rpath: '')

#python_sources = [
#    '__init__.py',
#    'pyCONMIN.py',
#    'LICENSE'
#]
#
#py3_target.install_sources(
#  python_sources,
#  pure: false,
#  subdir: 'pyoptsparse/pyCONMIN'
#)
