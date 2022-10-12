
python_sources = [
  '__init__.py',
  'OptView.py',
  'opt_view_controller.py',
]

py3_target.install_sources(
  python_sources,
  pure: true,
  subdir: 'pyoptsparse/postprocessing'
)

asset_sources = [
    'assets/nicePlotsStyle',
    'assets/optViewIcon.gif',
    'assets/mdoLabLogo.png',
    'assets/pyOptSparse_logo.png',
]

py3_target.install_sources(
  asset_sources,
  pure: true,
  subdir: 'pyoptsparse/postprocessing/assets'
)

subdir('baseclasses')
subdir('sub_windows')
subdir('utils')