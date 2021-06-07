# --- Python 3.8 ---
"""
Controller for the plot configuration and add variables view
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
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_view import VariableListWidget
from pyoptsparse.postprocessing.utils.data_structures import Variable
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_model import ConfigurePlotModel


class ConfigureController(object):
    def __init__(self, parent_model, plot_model):
        self._parent_model = parent_model
        self._plot_model = plot_model
        self._model = ConfigurePlotModel()
        self._view = None

    def set_view(self, view):
        self._view = view

    def file_selected(self):
        # --- Clear the comboboxes and variable lists ---
        self._view.x_cbox.clear()
        self._view.y_cbox.clear()
        self._view.x_list.clear()
        self._view.y_list.clear()

        idx = self._view.file_list.currentRow()
        file = self._parent_model.files[idx]
        var_names = file.get_all_variable_names()

        # --- Populate the combo boxes based on file selection ---
        for var in var_names:
            self._view.x_cbox.addItem(var)
            self._view.y_cbox.addItem(var)

        # --- Populate the variable lists based on file selection ---
        self.populate_x_list()
        self.populate_y_list()

    def populate_x_list(self):
        for i, var in enumerate(self._plot_model.x_vars):
            if var.file_idx == self._view.file_list.currentRow():
                x_item = QtWidgets.QListWidgetItem(self._view.x_list)
                x_item_widget = VariableListWidget(var.name, i, "x", self._view, self)
                x_item.setSizeHint(x_item_widget.sizeHint())
                self._view.x_list.addItem(x_item)
                self._view.x_list.setItemWidget(x_item, x_item_widget)

    def populate_y_list(self,):
        for var in self._plot_model.y_vars:
            if var.file_idx == self._view.file_list.currentRow():
                y_item = QtWidgets.QListWidgetItem(self._view.y_list)
                y_item_widget = VariableListWidget(var.name, var.var_idx, "y", self._view, self)
                y_item.setSizeHint(y_item_widget.sizeHint())
                self._view.y_list.addItem(y_item)
                self._view.y_list.setItemWidget(y_item, y_item_widget)

    def populate_files(self):
        for file in self._parent_model.files:
            self._view.file_list.addItem(file.name_short)

    def add_x_var(self):
        if len(self._plot_model.x_vars) < 1:
            # --- Get selected name from the combobox ---
            var_name = self._view.x_cbox.currentText()

            # --- Get the file index from the file list widget ---
            file_idx = self._view.file_list.currentRow()

            # --- Get the file at the specified index from from the parent model ---
            file = self._parent_model.files[file_idx]

            # --- Get a single variable from the file ---
            var_idx = len(self._plot_model.x_vars)  # absolute index of the variable in the plot model
            new_var = Variable(var_name, file_idx, var_idx)

            # --- Retrieve data from the file and store in variable ---
            file.get_single_variable(new_var)

            # --- Add variable to the plot model and local model ---
            self._plot_model.add_x_var(new_var)
            self._model.add_var(file_idx, "x")

            # --- Add the new variable to the view and set the identifier indices ---
            x_item = QtWidgets.QListWidgetItem(self._view.x_list)
            x_item_widget = VariableListWidget(var_name, var_idx, "x", self._view, self)
            x_item.setSizeHint(x_item_widget.sizeHint())
            self._view.x_list.addItem(x_item)
            self._view.x_list.setItemWidget(x_item, x_item_widget)
        else:
            QtWidgets.QMessageBox.warning(self._view, "X-Variable Warning", "OptView only supports 1 X-Variable")

    def add_y_var(self):
        # --- Get selected name from the combobox ---
        var_name = self._view.y_cbox.currentText()

        # --- Get the file index from the file list widget ---
        file_idx = self._view.file_list.currentRow()

        # --- Get the file at the specified index from from the parent model ---
        file = self._parent_model.files[file_idx]

        # --- Get a single variable from the file ---
        var_idx = len(self._plot_model.y_vars)  # absolute index of the variable in the plot model
        new_var = Variable(var_name, file_idx, var_idx)

        # --- Retrieve data from the file and store in variable ---
        file.get_single_variable(new_var)

        # --- Add variable to the plot model and local model ---
        self._plot_model.add_y_var(new_var)
        self._model.add_var(file_idx, "y")

        # --- Add the new variable to the view and set the identifier indices ---
        y_item = QtWidgets.QListWidgetItem(self._view.y_list)
        y_item_widget = VariableListWidget(var_name, var_idx, "y", self._view, self)
        y_item.setSizeHint(y_item_widget.sizeHint())
        self._view.y_list.addItem(y_item)
        self._view.y_list.setItemWidget(y_item, y_item_widget)

    def remove_variable(self, var_idx: int, axis: str):
        if axis == "x":
            # --- Remove variable from models ---
            self._plot_model.remove_x_var(var_idx)
            file_idx, rel_idx = self._model.remove_var(var_idx, "x")

            # --- Remove variable widget from the view ---
            self._view.x_list.takeItem(rel_idx)

            # --- Update the list widget view ---
            for i, val in enumerate(self._model.x_var_mapping):
                if val[0] == file_idx:
                    item = self._view.x_list.item(val[1])
                    widget = self._view.x_list.itemWidget(item)
                    widget.idx = i

        elif axis == "y":
            # --- Remove variable from models ---
            self._plot_model.remove_y_var(var_idx)
            file_idx, rel_idx = self._model.remove_var(var_idx, "y")

            # --- Remove variable widget from the view ---
            self._view.y_list.takeItem(rel_idx)

            # --- Update the list widget view ---
            for i, val in enumerate(self._model.y_var_mapping):
                if val[0] == file_idx:
                    item = self._view.y_list.item(val[1])
                    widget = self._view.y_list.itemWidget(item)
                    widget.idx = i

    def scale_opt_checked(self):
        if self._view.scale_opt.isChecked():
            self._plot_model.options["scaled"] = True
        else:
            self._plot_model.options["scaled"] = False

    def bounds_opt_checked(self):
        if self._view.bounds_opt.isChecked():
            self._plot_model.options["bounds"] = True
        else:
            self._plot_model.options["bounds"] = False

    def major_iter_toggled(self):
        if self._view.major_iter_togg.isChecked():
            print("True")
        else:
            print("False")

    def minor_iter_toggled(self):
        if self._view.minor_iter_togg.isChecked():
            print("True")
        else:
            print("False")

    def cancel(self):
        self._plot_model.x_vars = []
        self._plot_model.y_vars = []
        self._plot_model.options = {}
        self._view.close()
