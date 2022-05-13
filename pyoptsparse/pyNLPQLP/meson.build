#nlpqlp_source = custom_target('nlpqlpmodule.c',
#                            input : ['source/f2py/nlpqlp.pyf'],
#                            output : ['nlpqlpmodule.c'],
#                            command: [py3, '-m', 'numpy.f2py', '@INPUT@',
#                            '--lower', '--build-dir', 'pyoptsparse/pyNLPQLP']
#                            )
#
#py3.extension_module('nlpqlp',
#                'source/wrapper.F90',
#                nlpqlp_source,
#                fortranobject_c,
#                include_directories: [inc_np, inc_f2py],
#                dependencies : py3_dep,
#                subdir: 'pyoptsparse/pyNLPQLP',
#                install : true)

python_sources = [
    '__init__.py',
    'pyNLPQLP.py',
    'LICENSE'
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'pyoptsparse/pyNLPQLP'
)