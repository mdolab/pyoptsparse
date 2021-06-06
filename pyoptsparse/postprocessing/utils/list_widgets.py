# --- Python 3.8 ---
"""
Custom widget to display variable name, plot number, and close button
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


class VariableListWidget(QtWidgets.QWidget):
    def __init__(
        self,
        var_name: str = "Default Name",
        file_idx: int = 0,
        var_idx: int = 0,
        axis: str = "x",
        parent=None,
        controller=None,
    ):
        super(VariableListWidget, self).__init__(parent)

        self.controller = controller

        self.file_idx = file_idx
        self.var_idx = var_idx

        self.axis = axis

        self.label = QtWidgets.QLabel(var_name)

        self.remove_button = QtWidgets.QPushButton("X")
        self.remove_button.clicked.connect(self.remove)

        layout = QtWidgets.QHBoxLayout()

        layout.addWidget(self.label, 10)
        layout.addWidget(self.remove_button, 1)

        self.setLayout(layout)

    def remove(self):
        self.controller.remove_variable(self.var_idx, self.axis)


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
