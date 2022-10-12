
python_sources = [
  '__init__.py',
  'OptView.py',
  'opt_view_controller.py',
  'baseclasses.py',
  'data_structures.py',
  'colors.py',
  'general_widgets.py',
  'utils.py',
]

py3_target.install_sources(
  python_sources,
  pure: true,
  subdir: 'pyoptsparse/postprocessing'
)

asset_sources = [
    'assets/nicePlotsStyle',
    'assets/OptViewIcon.png',
    'assets/pyOptSparse_logo.png',
]

py3_target.install_sources(
  asset_sources,
  pure: true,
  subdir: 'pyoptsparse/postprocessing/assets'
)

subdir('sub_windows')