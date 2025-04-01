psqp_source = custom_target('psqpmodule.c',
                            input : ['source/f2py/psqp.pyf'],
                            output : ['psqpmodule.c', 'psqp-f2pywrappers.f'],
                            command: [py3, '-m', 'numpy.f2py', '@INPUT@',
                            '--lower', '--build-dir', 'pyoptsparse/pyPSQP', '--no-wrap-functions']
                            )

py3.extension_module('psqp',
                'source/closeunit.f',
                'source/mqsubs.f',
                'source/openunit.f',
                'source/pqsubs.f',
                'source/psqp.f',
                'source/psqp_wrap.f90',
                psqp_source,
                dependencies: [fortranobject_dep],
                subdir: 'pyoptsparse/pyPSQP',
                install: false)
