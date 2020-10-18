from sphinx_mdolab_theme.config import *

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("./"))  # to import custom Sphinx extension

# -- Project information -----------------------------------------------------

project = "pyOptSparse"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions.extend(["numpydoc", "ext.optimizertable"])
numpydoc_show_class_members = False
