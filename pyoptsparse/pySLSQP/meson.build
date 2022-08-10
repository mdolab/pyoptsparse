slsqp_source = custom_target('slsqpmodule.c',
                            input : ['source/f2py/slsqp.pyf'],
                            output : ['slsqpmodule.c'],
                            command: [py3_command, '-m', 'numpy.f2py', '@INPUT@',
                            '--lower', '--build-dir', 'pyoptsparse/pySLSQP']
                            )

py3_target.extension_module('slsqp',
                'source/closeunit.f',
                'source/daxpy.f',
                'source/dcopy.f',
                'source/drot.f',
                'source/drotg.f',
                'source/dscal.f',
                'source/h12.f',
                'source/hfti.f',
                'source/ldl.f',
                'source/ldp.f',
                'source/lsei.f',
                'source/lsi.f',
                'source/lsq.f',
                'source/nnls.f',
                'source/openunit.f',
                'source/slsqp.f',
                'source/slsqpb.f',
                slsqp_source,
                fortranobject_c,
                include_directories: [inc_np, inc_f2py],
                dependencies : py3_dep,
                subdir: 'pyoptsparse/pySLSQP',
                install : false)

#python_sources = [
#    '__init__.py',
#    'pySLSQP.py',
#    'LICENSE'
#]

#py3_target.install_sources(
#  python_sources,
#  pure: false,
#  subdir: 'pyoptsparse/pySLSQP'
#)
