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
from pyoptsparse.postprocessing.utils.list_widgets import VariableListWidget


class ConfigureController(object):
    def __init__(self, parent_model, plot_model):
        self._parent_model = parent_model
        self._plot_model = plot_model
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
        for var in self._plot_model.x_vars:
            if var.file_idx == self._view.file_list.currentRow():
                x_item = QtWidgets.QListWidgetItem(self._view.x_list)
                x_item_widget = VariableListWidget(var.name, var.file_idx, var.var_idx, "x", self._view, self)
                x_item.setSizeHint(x_item_widget.sizeHint())
                self._view.x_list.addItem(x_item)
                self._view.x_list.setItemWidget(x_item, x_item_widget)

    def populate_y_list(self,):
        for var in self._plot_model.y_vars:
            if var.file_idx == self._view.file_list.currentRow():
                y_item = QtWidgets.QListWidgetItem(self._view.y_list)
                y_item_widget = VariableListWidget(var.name, var.file_idx, var.var_idx, "y", self._view, self)
                y_item.setSizeHint(y_item_widget.sizeHint())
                self._view.y_list.addItem(y_item)
                self._view.y_list.setItemWidget(y_item, y_item_widget)

    def populate_files(self):
        for file in self._parent_model.files:
            self._view.file_list.addItem(file.name.split("/")[-1])

    def add_x_var(self):
        var_name = self._view.x_cbox.currentText()

        file_idx = self._view.file_list.currentRow()
        file = self._parent_model.files[file_idx]

        var = file.get_single_variable(var_name, file_idx)
        var.var_idx = var_idx = len(self._plot_model.x_vars)

        self._plot_model.add_x_var(var)

        x_item = QtWidgets.QListWidgetItem(self._view.x_list)
        x_item_widget = VariableListWidget(var_name, file_idx, var_idx, "x", self._view, self)
        x_item.setSizeHint(x_item_widget.sizeHint())
        self._view.x_list.addItem(x_item)
        self._view.x_list.setItemWidget(x_item, x_item_widget)

    def add_y_var(self):
        var_name = self._view.y_cbox.currentText()

        file_idx = self._view.file_list.currentRow()
        file = self._parent_model.files[file_idx]

        var = file.get_single_variable(var_name, file_idx)
        var.var_idx = var_idx = len(self._plot_model.y_vars)

        self._plot_model.add_y_var(var)

        y_item = QtWidgets.QListWidgetItem(self._view.y_list)
        y_item_widget = VariableListWidget(var_name, file_idx, var_idx, "y", self._view, self)
        y_item.setSizeHint(y_item_widget.sizeHint())
        self._view.y_list.addItem(y_item)
        self._view.y_list.setItemWidget(y_item, y_item_widget)

    def remove_variable(self, var_idx: int, axis: str):
        print(var_idx)
        if axis == "x":
            self._plot_model.remove_x_var(var_idx)
            self._view.x_list.takeItem(var_idx)

            # --- Loop over custom widgets and update the index and plot number ---
            for i in range(len(self._plot_model.x_vars)):
                item = self._view.x_list.item(i)
                widget = self._view.x_list.itemWidget(item)
                widget.var_idx = i
        elif axis == "y":
            self._plot_model.remove_y_var(var_idx)
            self._view.y_list.takeItem(var_idx)

            # --- Loop over custom widgets and update the index and plot number ---
            for i in range(len(self._plot_model.y_vars)):
                item = self._view.y_list.item(i)
                widget = self._view.y_list.itemWidget(item)
                widget.var_idx = i
