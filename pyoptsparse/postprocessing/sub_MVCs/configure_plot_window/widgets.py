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
from PyQt5 import QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================


class VariableListWidget(QtWidgets.QWidget):
    def __init__(self, var_name: str = "Default Name", idx: int = 0, parent=None, controller=None):
        """
        Custom widget for displaying a "mini-view" inside of a stock
        PyQt5 list widget.  We have to slightly deviate from the MVC
        architecture to properly link each widget with the parent
        controller in a way that allows data passing.  Passing the
        index allows us to track which variable is associated with this
        custom widget.  It

        Parameters
        ----------
        var_name : str, optional
            The variable name to be displayed, by default "Default Name"
        idx : int, optional
            The global variable index from the plot model, by default 0
        parent : PyQt5.QtWidgets.QWidget, optional
            The parent widget, by default None
        controller : controller object, optional
            The parent controller that handles the operations
            on the underlying variable data and parent view, by default
            None
        """
        super(VariableListWidget, self).__init__(parent)

        self._controller = controller

        self.idx = idx

        self.label = QtWidgets.QLabel(var_name)

        self.remove_button = QtWidgets.QPushButton("X")
        self.remove_button.clicked.connect(self.remove)

        self.scaled_opt = QtWidgets.QCheckBox("Scaled")
        self.scaled_opt.clicked.connect(self.scaled_opt_selected)

        self.bounds_opt = QtWidgets.QCheckBox("Bounds")
        self.bounds_opt.clicked.connect(self.bounds_opt_selected)

        self.minor_iter_opt = QtWidgets.QCheckBox("Minor Iters")
        self.minor_iter_opt.clicked.connect(self.minor_iter_opt_selected)

        layout = QtWidgets.QHBoxLayout()

        layout.addWidget(self.label, 1)
        layout.addWidget(self.scaled_opt, 1)
        layout.addWidget(self.bounds_opt, 1)
        layout.addWidget(self.minor_iter_opt, 1)
        layout.addWidget(self.remove_button, 1)

        self.setLayout(layout)

    def remove(self):
        self._controller.remove_variable(self.idx)

    def scaled_opt_selected(self):
        self._controller.scale_opt_selected(self, self.idx)

    def bounds_opt_selected(self):
        self._controller.bounds_opt_selected(self, self.idx)

    def minor_iter_opt_selected(self):
        self._controller.minor_iter_opt_selected(self, self.idx)
