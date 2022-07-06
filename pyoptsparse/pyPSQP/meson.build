psqp_source = custom_target('psqpmodule.c',
                            input : ['source/f2py/psqp.pyf'],
                            output : ['psqpmodule.c', 'psqp-f2pywrappers.f'],
                            command: [py3_command, '-m', 'numpy.f2py', '@INPUT@',
                            '--lower', '--build-dir', 'pyoptsparse/pyPSQP', '--no-wrap-functions']
                            )

py3_target.extension_module('psqp',
                'source/closeunit.f',
                'source/mqsubs.f',
                'source/openunit.f',
                'source/pqsubs.f',
                'source/psqp.f',
                'source/psqp_wrap.f90',
                psqp_source,
                fortranobject_c,
                include_directories: [inc_np, inc_f2py],
                dependencies : py3_dep,
                subdir: 'pyoptsparse/pyPSQP',
                install : false)

#python_sources = [
#    '__init__.py',
#    'pyPSQP.py',
#    'LICENSE'
#]
#
#py3_target.install_sources(
#  python_sources,
#  pure: false,
#  subdir: 'pyoptsparse/pyPSQP'
#)
