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
from pyoptsparse.postprocessing.sub_MVCs.widgets import FileTreeWidgetItem, VarTreeWidgetItem
from pyoptsparse.postprocessing.utils.data_structures import Variable
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_model import ConfigurePlotModel


class ConfigureController(object):
    def __init__(self, parent_model, plot_model):
        self._parent_model = parent_model
        self._plot_model = plot_model
        self._model = ConfigurePlotModel()
        self._view = None
        self._current_file = None

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
            self.populate_var_checkbox()

    def populate_vars(self):
        for var in self._plot_model.vars:
            var_item = VarTreeWidgetItem(self._view.var_tree)
            var_item.setVar(var)
            var_item.setText(0, var.file.name_short)
            var_item.setText(1, var.name)

            var_item.setCheckState(2, QtCore.Qt.Checked if var.options["scaled"] else QtCore.Qt.Unchecked)
            var_item.setCheckState(3, QtCore.Qt.Checked if var.options["bounds"] else QtCore.Qt.Unchecked)
            var_item.setCheckState(4, QtCore.Qt.Checked if var.options["major_iter"] else QtCore.Qt.Unchecked)
            self._view.var_tree.addTopLevelItem(var_item)

    def file_selected(self, item, column):
        self._current_file = item.file
        self.populate_var_checkbox()

    def populate_var_checkbox(self):
        self._view.var_cbox.clear()

        for var_name in self._current_file.get_all_var_names():
            self._view.var_cbox.addItem(var_name)

    def add_var(self):
        var_name = self._view.var_cbox.currentText()

        # Create a new variable and add to the plot model
        new_var = Variable(var_name)
        new_var.file = self._current_file
        self._plot_model.add_var(new_var)

        # Create a new variable widget item for the tree view
        new_var_item = VarTreeWidgetItem(self._view.var_tree)
        new_var_item.var = new_var

        new_var_item.setText(0, self._current_file.name_short)
        new_var_item.setText(1, var_name)
        new_var_item.setCheckState(2, QtCore.Qt.Unchecked)
        new_var_item.setCheckState(3, QtCore.Qt.Unchecked)
        new_var_item.setCheckState(4, QtCore.Qt.Unchecked)
        self._view.var_tree.addTopLevelItem(new_var_item)

    def remove_sel_var(self):
        sel_vars = self._view.var_tree.selectedItems()

        rem_idx = set()

        # Loop over all the selected variables and find them in the
        # plot model. We need to check both the filename and the
        # variable name in case the variable name exists in different
        # files.
        for i, plot_var in enumerate(self._plot_model.vars):
            for sel_var in sel_vars:
                if sel_var.var.name == plot_var.name and sel_var.var.file.name == plot_var.file.name:
                    rem_idx.add(i)

        # Remove variable from the plot
        self._plot_model.vars = [var for i, var in enumerate(self._plot_model.vars) if i not in rem_idx]

        # Remove variable from the tree view
        root = self._view.var_tree.invisibleRootItem()
        for item in sel_vars:
            root.removeChild(item)

    def cancel(self):
        self._plot_model.clear_vars()
        self._plot_model.clear_options()
        self._view.close()

    def plot(self):
        # We need to iterate over the variables in the view and set the
        # selected options in the plot variables.
        var_it = QtWidgets.QTreeWidgetItemIterator(self._view.var_tree)
        while var_it.value():
            item = var_it.value()
            item.setVarOpts()
            var_it += 1

        self._plot_model.plot()
        self._parent_model.canvas.draw()
        self._view.close()
