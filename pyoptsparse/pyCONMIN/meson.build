conmin_source = custom_target('conminmodule.c',
                            input : ['source/f2py/conmin.pyf',
                                    ],
                            output : ['conminmodule.c', 'conmin-f2pywrappers.f'],
                            command: [py3, '-m', 'numpy.f2py', '@INPUT@',
                            '--lower', '--build-dir', 'pyoptsparse/pyCONMIN']
                            )

py3.extension_module('conmin',
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
                install : true)

python_sources = [
    '__init__.py',
    'pyCONMIN.py',
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'pyoptsparse/pyCONMIN'
)