# --- Python 3.8 ---
"""
Controller for the main view.  Interacts with the data models and
handles all user input and response functionality.  Controller can
only update the view based on user input.  If a view state is changed
which requires a messagebox view, that view is created by the controller
but managed seperately.
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
from .model import HistoryFileModel


class SubWindowController:
    """
    Contains functionality for user input and software
    response for the main view.
    """

    def __init__(self, view):
        self._model = None
        self._view = view
        self._plot_controller = None

    def setPlotController(self, controller):
        self._plot_controller = controller

    def openFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(
            self._view, "Open History File", "", "History Files (*.sql)", options=options
        )

        if self._model is None:
            self._model = HistoryFileModel(file_name)
        else:
            self._model.changeFile(file_name)

        self._view.x_cbox.addItems(["test1", "test2", "test3"])
        self._view.y_cbox.addItems(["test1", "test2", "test3"])

    def newWindow(self):
        print("New window")

    def saveTecFile(self):
        print("Save Tec File")

    def refreshPlot(self):
        print("Refresh Plot")

    def clearPlot(self):
        if self._plot_controller is not None:
            self._plot_controller.clear()

    def addVarX(self):
        print("Add X Var")

    def addVarY(self):
        print("Add Y Var")

    def undoVarX(self):
        print("Undo X Var")

    def undoVarY(self):
        print("Undo Y Var")

    def clearAllX(self):
        print("Clear all X")

    def clearAllY(self):
        print("Clear all Y")

    def stackPlots(self):
        print("Stack Plots")

    def shareAxisY(self):
        print("Share Y axis")

    def absDeltaY(self):
        print("Absolute Delta Y axis")

    def minorIterX(self):
        print("Minor iters")

    def majorIterX(self):
        print("Major iters")

    def scaleFactor(self):
        print("Scale Factor")

    def minMaxPlot(self):
        print("Min/Max Plot")

    def autoRefresh(self):
        print("Auto Refresh")
