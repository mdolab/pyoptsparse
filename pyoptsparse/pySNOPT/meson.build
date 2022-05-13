#snopt_source = custom_target('snoptmodule.c',
#                            input : ['source/f2py/snopt.pyf'],
#                            output : ['snoptmodule.c'],
#                            command: [py3, '-m', 'numpy.f2py', '@INPUT@',
#                            '--lower', '--build-dir', 'pyoptsparse/pySNOPT']
#                            )
                            #
#py3.extension_module('snopt',
#                'source/closeunit.f',
#                'source/openunit.f',
#                snopt_source,
#                fortranobject_c,
#                include_directories: [inc_np, inc_f2py],
#                dependencies : py3_dep,
#                subdir: 'pyoptsparse/pySNOPT',
#                install : true)

python_sources = [
    '__init__.py',
    'pySNOPT.py',
    'LICENSE'
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'pyoptsparse/pySNOPT'
)