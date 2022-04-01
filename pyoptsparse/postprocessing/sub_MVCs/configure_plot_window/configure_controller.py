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
from PyQt5 import QtWidgets, QtCore

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.sub_MVCs.widgets import FileTreeWidgetItem, VarTableWidgetItem, FileTableWidgetItem
from pyoptsparse.postprocessing.utils.data_structures import Variable
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_model import ConfigureModel
from pyoptsparse.postprocessing.sub_MVCs.variables.y_controller import YController
from pyoptsparse.postprocessing.sub_MVCs.variables.x_controller import XController


class ConfigureController(object):
    def __init__(self, parent_model, plot_model):
        self._parent_model = parent_model
        self._plot_model = plot_model
        self._model = ConfigureModel()
        self._view = None
        self._current_file = None

        self._xtable_controller = None
        self._ytable_controller = None

    def setup_var_tables(self):
        self._xtable_controller = XController(self._view.x_table, self._plot_model, self._parent_model)
        self._ytable_controller = YController(self._view.y_table, self._plot_model, self._parent_model)

        self._view.x_table.setController(self._xtable_controller)
        self._view.y_table.setController(self._ytable_controller)

    def set_view(self, view):
        self._view = view

    def add_file(self):
        # --- Set file dialog options ---
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog

        # --- Open file dialog and get selected user files ---
        file_names, _ = QtWidgets.QFileDialog.getOpenFileNames(self._view, "Open History File", "", "", options=options)

        # --- Load files into the model ---
        self._parent_model.load_files(file_names)

        # --- Populate the files ---
        self.populate_files()

    def populate_files(self):
        for file in self._parent_model.files:
            file_item = FileTreeWidgetItem(self._view.file_tree)
            file_item.setFile(file)
            file_item.setText(0, file.name_short)
            self._view.file_tree.addTopLevelItem(file_item)

        if len(self._parent_model.files) > 0:
            self._current_file = self._parent_model.files[0]
            self.populate_vars()

    def file_selected(self, item, column):
        self._current_file = item.file
        self.populate_vars()

    def populate_vars(self):
        self.populate_x_var_checkbox()
        self._view.y_table.clear()
        self._view.x_table.clear()
        self._ytable_controller.populate_vars(self._current_file)

        if self._plot_model.x_var is None:
            new_var = Variable("iter")
            new_var.file = self._current_file
            self._plot_model.add_var(new_var, "x")

        self._xtable_controller.populate_vars()

    def populate_x_var_checkbox(self):
        self._view.x_cbox.clear()

        for var_name in self._current_file.get_all_x_var_names():
            self._view.x_cbox.addItem(var_name)

    def add_x_var(self):
        self._view.x_table.clear()
        var_name = self._view.x_cbox.currentText()

        # Create a new variable and add to the plot model
        var = Variable(var_name)
        var.file = self._current_file
        self._plot_model.add_var(var, "x")

        # Create a new variable widget item for the tree view
        file_item = FileTableWidgetItem(var.file.name_short)
        var_item = VarTableWidgetItem(var.name)
        var_item.var = var

        self._xtable_controller.add_row(file_item, var_item)

    def y_var_search(self, s):
        items = self._view.y_table.findItems(s, QtCore.Qt.MatchContains)
        if items:
            item = items[1]
            self._view.y_table.setCurrentItem(item)

    def cancel(self):
        self._plot_model.clear_vars()
        self._plot_model.clear_options()
        self._view.close()
