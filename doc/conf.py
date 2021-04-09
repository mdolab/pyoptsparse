# Standard Python modules
import os
import sys

# External modules
from sphinx_mdolab_theme.config import *

# First party modules
import pyoptsparse

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.


sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("./"))  # to import custom Sphinx extension

# -- Project information -----------------------------------------------------

project = "pyOptSparse"
version = pyoptsparse.__version__

# -- General configuration ---------------------------------------------------

# mock import for autodoc
autodoc_mock_imports = ["baseclasses", "scipy", "numpy", "sqlitedict"]

# logo
html_logo = "_static/pyOptSparse_logo.svg"
html_theme_options["logo_only"] = True

# bibtex
bibtex_bibfiles = ["pyoptsparse.bib"]
