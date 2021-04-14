# --- Python 3.8 ---
"""
Controller for the sub view.  Interacts with the data models and
handles all user input and response functionality.  Controller can
only update the view based on user input.  If a view state is changed
which requires a messagebox view, that view is created by the controller
but managed seperately.

State control is encapsulated within it's own controller class which
is specific to this sub view.
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
        self._plot_options = {"standard": False, "stacked": False}

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
        # --- Clear combobox and labels ---
        self._view.x_cbox.clear()
        self._view.y_cbox.clear()
        self._view.x_label.clear()
        self._view.y_label.clear()

        # --- Uncheck all options ---
        self._view.stack_plot_opt.setChecked(False)
        self._view.share_x_opt.setChecked(False)
        self._view.min_max_opt.setChecked(False)
        self._view.bound_opt.setChecked(False)
        self._view.minor_itr_opt.setChecked(False)
        self._view.major_itr_opt.setChecked(False)
        self._view.abs_delta_opt.setChecked(False)

        # --- State control ---
        self._state_controller.setInitialState()

    def refreshFile(self):
        print("refresh file")

    def clearPlot(self):
        # --- Clear plot using plot controller ---
        self._plot_controller.clear()

        # --- Reset plotting options ---
        # These will be set again if user chooses to plot using
        # same options checked
        for key in self._plot_options:
            self._plot_options[key] = False

    def addVarX(self):
        x_var = self._view.x_cbox.currentText()
        self._model.addX(x_var)
        x_names = "".join([str(i) + "\n" for i in self._model.x_vars.keys()])
        self._view.x_label.setText(x_names)

        # --- State Control ---
        if len(self._model.x_vars) > 0 and len(self._model.y_vars) > 0:
            self._state_controller.setPlotState(False)

        if len(self._model.y_vars) > 1 and len(self._model.x_vars) > 0:
            self._state_controller.setStackedPlotState(False)

        if len(self._model.y_vars) > 0 and len(self._model.x_vars) > 1:
            self._state_controller.setStackedPlotState(False)

    def addVarY(self):
        y_var = self._view.y_cbox.currentText()
        self._model.addY(y_var)
        y_names = "".join([str(i) + "\n" for i in self._model.y_vars.keys()])
        self._view.y_label.setText(y_names)

        # --- State control ---
        if len(self._model.x_vars) > 0 and len(self._model.y_vars) > 0:
            self._state_controller.setPlotState(False)

        if len(self._model.y_vars) > 0 and len(self._model.x_vars) > 1:
            self._state_controller.setStackedPlotState(False)

        if len(self._model.y_vars) > 1 and len(self._model.x_vars) > 0:
            self._state_controller.setStackedPlotState(False)

        if len(self._model.y_vars) > 0:
            self._state_controller.setAbsDeltaState(False)

    def undoVarX(self):
        self._model.undoX()
        x_names = "".join([str(i) + "\n" for i in self._model.x_vars.keys()])
        self._view.x_label.setText(x_names)

        # --- State control ---
        if len(self._model.x_vars) < 1:
            self._state_controller.setStackedPlotState(True)
            self._state_controller.setPlotState(True)

    def undoVarY(self):
        self._model.undoY()
        y_names = "".join([str(i) + "\n" for i in self._model.y_vars.keys()])
        self._view.y_label.setText(y_names)

        # --- State control ---
        if len(self._model.y_vars) < 1:
            self._state_controller.setStackedPlotState(True)

        if len(self._model.y_vars) < 1:
            self._state_controller.setAbsDeltaState(True)
            self._state_controller.setPlotState(True)

    def clearAllX(self):
        self._model.clearX()
        self._view.x_label.setText("")

        # --- State control ---
        self._state_controller.setClearVarState()

    def clearAllY(self):
        self._model.clearY()
        self._view.y_label.setText("")

        # --- State control ---
        self._state_controller.setClearVarState()

    def scaleVars(self):
        if self._view.scale_var_togg.isChecked():
            self._model.scaleY()
            self.plot()
        else:
            self._model.unscaleY()
            self.plot()

    def autoRefresh(self):
        print("Auto Refresh")

    def majorMinorIterX(self):
        # --- Only major iteration is checked ---
        if self._view.major_itr_opt.isChecked() and not self._view.minor_itr_opt.isChecked():
            # --- State control ---
            self._state_controller.setMajorMinorIterCheckedState()

            if len(self._model.y_vars) > 0:
                self._state_controller.setPlotState(False)

            # --- clear x-data and set major iterations as x-data ---
            self._model.clearX()
            self._view.x_label.setText("Major Iterations")
            self._view.x_label.setAlignment(QtCore.Qt.AlignCenter)
            self._model.majorIterX()

        # --- Only minor iteration is checked ---
        elif not self._view.major_itr_opt.isChecked() and self._view.minor_itr_opt.isChecked():
            # --- State control ---
            self._state_controller.setMajorMinorIterCheckedState()

            if len(self._model.y_vars) > 0:
                self._state_controller.setPlotState(False)

            # --- clear x-data and set major iterations as x-data ---
            self._model.clearX()
            self._view.x_label.setText("Minor Iterations")
            self._view.x_label.setAlignment(QtCore.Qt.AlignCenter)
            self._model.minorIterX()

        # --- Both iteration types are checked ---
        elif self._view.major_itr_opt.isChecked() and self._view.minor_itr_opt.isChecked():
            # --- State control ---
            self._state_controller.setMajorMinorIterCheckedState()

            if len(self._model.y_vars) > 0:
                self._state_controller.setPlotState(False)

            # --- clear x-data and set major iterations as x-data ---
            self._model.clearX()
            self._view.x_label.setText("Major Iterations\nMinor Iterations")
            self._view.x_label.setAlignment(QtCore.Qt.AlignCenter)
            self._model.majorIterX()
            self._model.minorIterX()

        # --- No iteration types are checked ---
        else:
            # --- Unset x-data ---
            self._model.clearX()
            self._view.x_label.clear()

            # --- State control ---
            self._state_controller.setMajorMinorIterUncheckedState()
            self._state_controller.setPlotState(True)

    def plot(self):
        if len(self._model.x_vars) == 1 and len(self._model.y_vars) > 0:
            self._plot_options["standard"] = True
        if len(self._model.x_vars) > 1 or len(self._model.y_vars) > 1 and self._view.stack_plot_opt.isChecked():
            self._plot_options["stacked"] = True
        if self._view.abs_delta_opt.isChecked():
            # TODO: Plot abs delta values of the y-variables
            pass
        if self._view.min_max_opt.isChecked():
            # TODO: Plot min and max of each variable/function
            pass
        if self._view.bound_opt.isChecked():
            # TODO: Plot the bounds for each y-variable
            pass
