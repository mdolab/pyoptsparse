# --- Python 3.8 ---
"""
Main OptView controller.  Communicates between the main view and the
different plotting tabs created by the user.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5 import QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.sub_MVCs.tab_view import TabView


class OptViewController:
    def __init__(self, view):
        self._view = view

    def addTab(self):
        tab_name, ok_pressed = QtWidgets.QInputDialog.getText(
            self._view, "Enter Tab Name", "Tab Name:", QtWidgets.QLineEdit.Normal, ""
        )
        tab = TabView(self._view)
        self._view.tabs.addTab(tab, tab_name)

    def closeTab(self, current_index):
        self._view.tabs.removeTab(current_index)
