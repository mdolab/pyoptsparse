# --- Python 3.8 ---
"""
OptView - A GUI designed to assist viewing optimization problems with
pyOptSparse.  The app uses a MVC architecure to modularly handle
data management, functional control, and visualization.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import sys

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5 import QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================
from views.combo_box import ExtendedComboBox
from views.mpl_canvas import MplCanvas
