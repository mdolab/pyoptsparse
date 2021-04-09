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
from PyQt5 import QtWidgets, QtCore

# ==============================================================================
# Extension modules
# ==============================================================================
from .model import HistoryFileModel
from .state_controller import StateController


class SubWindowController:
    """
    Contains functionality for user input and software
    response for the main view.
    """

    def __init__(self, view):
        self._model = None
        self._view = view
        self._plot_controller = None
        self._state_controller = StateController(view)
        self._plot_options = 1

    def setInitialState(self):
        self._state_controller.setInitialState()

    def setPlotController(self, controller):
        self._plot_controller = controller

    def openFile(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(
            self._view, "Open History File", "", "History Files (*.sql)", options=options
        )
        return file_name

    def addFile(self):
        # If there is no model, then we need to load initial model and
        # data
        if self._model is None:
            file_name = self.openFile()
            self._model = HistoryFileModel(file_name)
            self._view.x_cbox.addItems(self._model.getNames())
            self._view.y_cbox.addItems(self._model.getNames())

        # If a model already exists, we prompt the user if they want
        # to clear all current data and load new data and new model
        else:
            buttonReply = QtWidgets.QMessageBox.question(
                self._view,
                "New File Warning",
                "Adding new file will lose old file data.\nDo you want to continue?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel,
                QtWidgets.QMessageBox.Cancel,
            )
            # If user clicks yes button, the view and model are reset
            if buttonReply == QtWidgets.QMessageBox.Yes:
                self.reset()
                file_name = self.openFile()
                self._model.changeFile(file_name)

                self._view.x_cbox.addItems(self._model.getNames())
                self._view.y_cbox.addItems(self._model.getNames())

        self._state_controller.setAddFileState()

    def reset(self):
        self._view.x_cbox.clear()
        self._view.y_cbox.clear()
        self._view.x_label.clear()
        self._view.y_label.clear()
        self._view.stack_plot_opt.setChecked(False)
        self._view.share_x_opt.setChecked(False)
        self._view.min_max_opt.setChecked(False)
        self._view.bound_opt.setChecked(False)
        self._view.minor_itr_opt.setChecked(False)
        self._view.major_itr_opt.setChecked(False)
        self._view.abs_delta_opt.setChecked(False)

        self._state_controller.setInitialState()

    def refreshFile(self):
        print("refresh file")

    def clearPlot(self):
        self._plot_options = 1
        self._plot_controller.clear()

    def addVarX(self):
        x_var = self._view.x_cbox.currentText()
        self._model.addX(x_var)
        x_names = "".join([str(i) + "\n" for i in self._model.x_vars.keys()])
        self._view.x_label.setText(x_names)

    def addVarY(self):
        y_var = self._view.y_cbox.currentText()
        self._model.addY(y_var)
        y_names = "".join([str(i) + "\n" for i in self._model.y_vars.keys()])
        self._view.y_label.setText(y_names)

    def undoVarX(self):
        self._model.undoX()
        x_names = "".join([str(i) + "\n" for i in self._model.x_vars.keys()])
        self._view.x_label.setText(x_names)

    def undoVarY(self):
        self._model.undoY()
        y_names = "".join([str(i) + "\n" for i in self._model.y_vars.keys()])
        self._view.y_label.setText(y_names)

    def clearAllX(self):
        self._model.clearX()
        self._view.x_label.setText("")

    def clearAllY(self):
        self._model.clearY()
        self._view.y_label.setText("")

    def scaleFactor(self):
        self._model.scaleY()

    def autoRefresh(self):
        print("Auto Refresh")

    def majorMinorIterX(self):
        # --- Only major iteration is checked ---
        if self._view.major_itr_opt.isChecked() and not self._view.minor_itr_opt.isChecked():
            # --- State control ---
            self._state_controller.setMajorMinorIterCheckedState()

            # --- clear x-data and set major iterations as x-data ---
            self._model.clearX()
            self._view.x_label.setText("Major Iterations")
            self._view.x_label.setAlignment(QtCore.Qt.AlignCenter)
            self._model.majorIterX()

        # --- Only minor iteration is checked ---
        elif not self._view.major_itr_opt.isChecked() and self._view.minor_itr_opt.isChecked():
            # --- State control ---
            self._state_controller.setMajorMinorIterCheckedState()

            # --- clear x-data and set major iterations as x-data ---
            self._model.clearX()
            self._view.x_label.setText("Minor Iterations")
            self._view.x_label.setAlignment(QtCore.Qt.AlignCenter)
            self._model.minorIterX()

        # --- Both iteration types are checked ---
        elif self._view.major_itr_opt.isChecked() and self._view.minor_itr_opt.isChecked():
            # --- State control ---
            self._state_controller.setMajorMinorIterCheckedState()

            # --- clear x-data and set major iterations as x-data ---
            self._model.clearX()
            self._view.x_label.setText("Major Iterations\nMinor Iterations")
            self._view.x_label.setAlignment(QtCore.Qt.AlignCenter)
            self._model.majorIterX()
            self._model.minorIterX()

        # --- No iteration types are checked ---
        else:
            # --- State control ---
            self._state_controller.setMajorMinorIterUncheckedState()

            # --- Unset x-data ---
            self._model.clearX()
            self._view.x_label.clear()

    def plot(self):
        pass
