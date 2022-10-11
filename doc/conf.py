# Standard Python modules
import os
import re
import sys

# External modules
from sphinx_mdolab_theme.config import *

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.


sys.path.insert(0, os.path.abspath("../"))


# -- Project information -----------------------------------------------------

project = "pyOptSparse"
version = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("../pyoptsparse/__init__.py").read(),
)[0]

# -- General configuration ---------------------------------------------------

# mock import for autodoc
autodoc_mock_imports = ["baseclasses", "scipy", "numpy", "sqlitedict"]

# logo
html_logo = "_static/pyOptSparse_logo.svg"
html_theme_options["logo_only"] = True

# bibtex
bibtex_bibfiles = ["pyoptsparse.bib"]
