# Standard Python modules
import os
import sys

# External modules
from setuptools_scm import get_version
from sphinx_mdolab_theme.config import *

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.


sys.path.insert(0, os.path.abspath("../"))


# -- Project information -----------------------------------------------------

project = "pyOptSparse"
version = get_version(root="..", relative_to=__file__)

# -- General configuration ---------------------------------------------------

# mock import for autodoc
autodoc_mock_imports = ["baseclasses", "sqlitedict"]

# logo
html_logo = "_static/pyOptSparse_logo.svg"
html_theme_options["logo_only"] = True

# bibtex
bibtex_bibfiles = ["pyoptsparse.bib"]

# autolink
extensions.extend(["sphinx_codeautolink"])
codeautolink_concat_default = True
codeautolink_warn_on_missing_inventory = True
codeautolink_warn_on_failed_resolve = True
