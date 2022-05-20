fs = import('fs')
if fs.is_file('source' / 'NLPQLP.F')
  message('NLPQLP source code found!')
  HAS_SNOPT = true
else
  message('NLPQLP not found!')
  HAS_SNOPT = false
endif

if HAS_SNOPT
  nlpqlp_source = custom_target('nlpqlpmodule.c',
    input : ['source/f2py/nlpqlp.pyf'],
    output : ['nlpqlpmodule.c'],
    command: [py3_command, '-m', 'numpy.f2py', '@INPUT@',
    '--lower', '--build-dir', 'pyoptsparse/pyNLPQLP']
  )

  py3_target.extension_module('nlpqlp',
    'source/wrapper.F90',
    'source/NLPQLP.F',
    'source/QL.F',
    nlpqlp_source,
    '../fortranobject.c',
    include_directories: [inc_np, inc_f2py],
    dependencies : py3_dep,
    subdir: 'pyoptsparse/pyNLPQLP',
    install : false
  )
endif

#python_sources = [
#    '__init__.py',
#    'pyNLPQLP.py',
#    'LICENSE'
#]
#
#py3_target.install_sources(
#  python_sources,
#  pure: false,
#  subdir: 'pyoptsparse/pyNLPQLP'
#)
