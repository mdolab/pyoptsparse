fs = import('fs')
if fs.is_file('source' / 'snoptc.f')
  message('SNOPT source code found!')
  HAS_SNOPT = true
else
  message('SNOPT not found!')
  HAS_SNOPT = false
endif

if HAS_SNOPT
  c = run_command('source' / 'grab-all-fortran-files.py', check: true)
  snopt_source_files = c.stdout().strip().split('\n')

  snopt_source = custom_target('snoptmodule.c',
    input : ['source/f2py/snopt.pyf'],
    output : ['snoptmodule.c'],
    command: [py3_command, '-m', 'numpy.f2py', '@INPUT@',
    '--lower', '--build-dir', 'pyoptsparse/pySNOPT']
  )

  py3_target.extension_module('snopt',
    snopt_source,
    fortranobject_c,
    snopt_source_files,
    include_directories: [inc_np, inc_f2py],
    dependencies : py3_dep,
    subdir: 'pyoptsparse/pySNOPT',
    install : false,
    fortran_args: '-ffixed-line-length-80'
  )
endif

#python_sources = [
#    '__init__.py',
#    'pySNOPT.py',
#    'LICENSE'
#]

#py3_target.install_sources(
#  python_sources,
#  pure: false,
#  subdir: 'pyoptsparse/pySNOPT'
#)