# --- Python 3.8 ---
"""
Custom widgets for the configure plot window
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5 import QtWidgets, QtCore

# ==============================================================================
# Extension modules
# ==============================================================================


class FileTreeWidgetItem(QtWidgets.QTreeWidgetItem):
    def __init__(self, *args, **kwargs):
        self.file = None
        super().__init__(*args, **kwargs)

    def setFile(self, file):
        self.file = file


class OptTreeWidgetItem(QtWidgets.QTreeWidgetItem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class OptTableWidgetItem(QtWidgets.QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class VarTreeWidgetItem(QtWidgets.QTreeWidgetItem):
    def __init__(self, *args, **kwargs):
        self.var = None
        super().__init__(*args, **kwargs)

    def setVar(self, var):
        self.var = var

    def setVarOpts(self):
        self.var.options["scaled"] = self.checkState(2) == QtCore.Qt.Checked
        self.var.options["bounds"] = self.checkState(3) == QtCore.Qt.Checked
        self.var.options["major_iter"] = self.checkState(4) == QtCore.Qt.Checked
