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
from PyQt5 import QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================
from .model import Model


class SubWindowController:
    """
    Contains functionality for user input and software
    response for the main view.
    """

    def __init__(self, view):
        self._model = Model()
        self._view = view
        self._plot_controller = None
        self._options = {"standard": False, "stacked": False}

    def set_plot_controller(self, controller):
        self._plot_controller = controller

    def open_file(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        file_name, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self._view, "Open History File", "", "History Files (*.hst);; SQL File (*.sql)", options=options,
        )
        return file_name

    def reset_window(self):
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

    def refresh_file(self):
        print("refresh file")

    def clear_plot(self):
        # --- Clear plot using plot controller ---
        self._plot_controller.clear()

        # --- Reset plotting options ---
        print("Clear plotting options")

    def add_x_var(self):
        pass
        # var_name = self._view.x_cbox.currentText()  # get the text from the combobox
        # self._model.add_x_var(var_name)  # add the xvar to the model

    def add_y_var(self):
        pass
        # y_var = self._view.y_cbox.currentText()
        # self._model.add_y_var(y_var)

    def clear_x(self):
        self._model.clear_x_vars()

    def clear_y(self):
        self._model.clear_y_vars()

    def scale_vars(self):
        if self._view.scale_var_togg.isChecked():
            self._model.scaleY()
            self.plot()
        else:
            self._model.unscaleY()
            self.plot()

    def auto_refresh(self):
        print("Auto Refresh")

    def plot(self):
        pass
