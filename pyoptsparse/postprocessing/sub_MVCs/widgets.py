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
from pyoptsparse.postprocessing.utils.button import Button
from pyoptsparse.postprocessing.utils.switch import Switch


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
        self.controller.configure_view(self.idx)


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


class VarTableWidget(QtWidgets.QTableWidget):
    def __init__(self, parent=None):
        super(VarTableWidget, self).__init__(parent)

    def addYvarRow(self, file_item, var_item):
        row = self.rowCount()
        self.setRowCount(row + 1)
        self.setItem(row, 0, file_item)
        self.setItem(row, 1, var_item)

        add_btn = Button("Add", self)
        self.setCellWidget(row, 2, add_btn)

        rem_btn = Button("Remove", self)
        self.setCellWidget(row, 3, rem_btn)

    def addXvarRow(self, file_item, var_item, row):
        row = self.rowCount()
        self.setRowCount(row + 1)
        self.setItem(row, 0, file_item)
        self.setItem(row, 1, var_item)


class FileTableWidgetItem(QtWidgets.QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        super(FileTableWidgetItem, self).__init__(*args, **kwargs)


class VarTableWidgetItem(QtWidgets.QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        super(VarTableWidgetItem, self).__init__(*args, **kwargs)
        self.var = None

    def setVar(self, var):
        self.var = var


class TableButtonWidget(
    QtWidgets.QPushButton,
):
    def __init__(self, row, *args, **kwargs):
        super(TableButtonWidget, self).__init__(*args, **kwargs)
        self.row = row


class IterSwitchWidget(Switch):
    def __init__(self, row: int, *args, **kwargs):
        super(IterSwitchWidget, self).__init__(*args, **kwargs)
        self.row = row
