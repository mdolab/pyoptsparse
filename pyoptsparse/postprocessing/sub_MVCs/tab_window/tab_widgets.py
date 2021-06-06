# --- Python 3.8 ---
"""
Custom widgets only used for the tab window
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
