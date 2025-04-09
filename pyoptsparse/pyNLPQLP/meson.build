fs = import('fs')
if fs.is_file('source' / 'NLPQLP.F')
  message('NLPQLP source code found!')
  HAS_NLPQLP = true
else
  message('NLPQLP not found!')
  HAS_NLPQLP = false
endif

if HAS_NLPQLP
  nlpqlp_source = custom_target('nlpqlpmodule.c',
    input : ['source/f2py/nlpqlp.pyf'],
    output : ['nlpqlpmodule.c'],
    command: [py3, '-m', 'numpy.f2py', '@INPUT@',
    '--lower', '--build-dir', 'pyoptsparse/pyNLPQLP']
  )

  py3.extension_module('nlpqlp',
    'source/wrapper.F90',
    'source/NLPQLP.F',
    'source/QL.F',
    nlpqlp_source,
    dependencies: [fortranobject_dep],
    subdir: 'pyoptsparse/pyNLPQLP',
    install: true
  )
endif
