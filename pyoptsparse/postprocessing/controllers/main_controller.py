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


class MainController:
    """
    Contains functionality for user input and software
    response for the main view.
    """

    def __init__(self, view):
        self._plot_controller = None
        self._model = None
        self._view = view

    def openFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(
            self._view, "Open History File", "", "History Files (*.sql)", options=options
        )
        # TODO: Set model file name and load variable names

    def newWindow(self):
        print("New window")

    def saveTecFile(self):
        print("Save Tec File")

    def changePlotFontSize(self, value):
        print("Change Font")

    def refreshPlot(self):
        print("Refresh Plot")

    def clearPlot(self):
        print("Clear Plot")

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
