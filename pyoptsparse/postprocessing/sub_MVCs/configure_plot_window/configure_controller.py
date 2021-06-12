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
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.widgets import VariableListWidget
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
        self._view.var_cbox.clear()
        self._view.var_list.clear()

        idx = self._view.file_list.currentRow()
        file = self._parent_model.files[idx]
        var_names = file.get_all_variable_names()

        # --- Populate the combo boxes based on file selection ---
        for var in var_names:
            self._view.var_cbox.addItem(var)

        # --- Populate the variable lists based on file selection ---
        self.populate_var_list()

    def populate_var_list(self,):
        for i, var in enumerate(self._plot_model.vars):
            if var.file_idx == self._view.file_list.currentRow():
                var_item = QtWidgets.QListWidgetItem(self._view.var_list)
                var_item_widget = VariableListWidget(var.name, i, self._view, self)
                var_item.setSizeHint(var_item_widget.sizeHint())
                self._view.var_list.addItem(var_item)
                self._view.var_list.setItemWidget(var_item, var_item_widget)

    def populate_files(self):
        for file in self._parent_model.files:
            self._view.file_list.addItem(file.name_short)

    def add_var(self):
        # --- Get selected name from the combobox ---
        var_name = self._view.var_cbox.currentText()

        # --- Get the file index from the file list widget ---
        file_idx = self._view.file_list.currentRow()

        # --- Get the file at the specified index from from the parent model ---
        file = self._parent_model.files[file_idx]

        # --- Add the variable to the overall list ---
        self._view.all_vars_list.addItem(f"{file.name_short} | {var_name}")

        # --- Get a single variable from the file ---
        var_idx = len(self._plot_model.vars)  # absolute index of the variable in the plot model
        new_var = Variable(var_name, file_idx, var_idx)

        # --- Retrieve data from the file and store in variable ---
        file.get_single_variable(new_var)

        # --- Add variable to the plot model and local model ---
        self._plot_model.add_var(new_var)
        self._model.add_var(file_idx)

        # --- Add the new variable to the view and set the identifier indices ---
        var_item = QtWidgets.QListWidgetItem(self._view.var_list)
        var_item_widget = VariableListWidget(var_name, var_idx, self._view, self)
        var_item.setSizeHint(var_item_widget.sizeHint())
        self._view.var_list.addItem(var_item)
        self._view.var_list.setItemWidget(var_item, var_item_widget)

    def remove_variable(self, var_idx: int):
        # --- Remove variable from overall list ---
        self._view.all_vars_list.takeItem(var_idx)

        # --- Remove variable from models ---
        self._plot_model.remove_var(var_idx)
        file_idx, rel_idx = self._model.remove_var(var_idx)

        # --- Remove variable widget from the view ---
        self._view.var_list.takeItem(rel_idx)

        # --- Update the list widget view ---
        for i, val in enumerate(self._model.var_mapping):
            if val[0] == file_idx:
                item = self._view.var_list.item(val[1])
                widget = self._view.var_list.itemWidget(item)
                widget.idx = i

    def scale_opt_selected(self, widget: VariableListWidget, var_idx: int):
        if widget.scaled_opt.isChecked():
            self._plot_model.vars[var_idx].options["scaled"] = True
            self._view.all_vars_list[var_idx].strip(" ").split("|")
        else:
            self._plot_model.vars[var_idx].options["scaled"] = False

    def bounds_opt_selected(self, widget: VariableListWidget, var_idx: int):
        if widget.bounds_opt.isChecked():
            self._plot_model.vars[var_idx].options["bounds"] = True
        else:
            self._plot_model.vars[var_idx].options["scaled"] = False

    def cancel(self):
        self._plot_model.clear_vars()
        self._plot_model.clear_options()
        self._view.close()
