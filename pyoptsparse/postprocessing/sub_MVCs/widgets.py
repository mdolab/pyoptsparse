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


class PlotListWidget(QtWidgets.QWidget):
    def __init__(self, parent=None, controller=None, idx: int = 0):
        super(PlotListWidget, self).__init__(parent)

        self.controller = controller

        self.idx = idx

        self.title = QtWidgets.QLineEdit(f"Plot {idx}")

        self.configure_button = QtWidgets.QPushButton("Configure/Add Variables")
        self.configure_button.clicked.connect(self.configure)

        self.remove_button = QtWidgets.QPushButton("Remove Plot")
        self.remove_button.clicked.connect(self.remove)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.title, 1)
        layout.addWidget(self.configure_button, 1)
        layout.addWidget(self.remove_button, 1)

        self.setLayout(layout)

    def remove(self):
        self.controller.remove_plot(self.idx)

    def configure(self):
        self.controller.configure_view(self.idx, self.title.text())


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
