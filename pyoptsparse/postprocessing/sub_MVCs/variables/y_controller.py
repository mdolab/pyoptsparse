#!/usr/bin/env python
"""
@File    :   y_controller.py
@Time    :   2022/03/31
@Desc    :   None
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5 import QtCore, QtGui, QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.utils.data_structures import Variable
from pyoptsparse.postprocessing.sub_MVCs.widgets import TableButtonWidget as Button
from pyoptsparse.postprocessing.sub_MVCs.widgets import VarTableWidgetItem, FileTableWidgetItem


GREEN = QtGui.QColor(0, 255, 0, 20)
WHITE = QtGui.QColor(255, 255, 255)


class YController:
    def __init__(self, view, plot_model, parent_model):
        self._view = view
        self._parent_model = parent_model
        self._plot_model = plot_model

    def populate_vars(self, current_file):
        for name in current_file.get_all_y_var_names():
            var = Variable(name)
            var.file = current_file

            # Create a new variable widget item for the tree view
            var_item = VarTableWidgetItem(var.name)
            file_item = FileTableWidgetItem(var.file.name_short)
            var_item.var = var

            self.add_row(file_item, var_item)

        # Find all variables that are in the plot, highlight them in
        # green, and link the table variable to the plot variable
        for plot_var in self._plot_model.y_vars:
            for i in self._view.rowCount():
                var_item = self._view.item(i, 1)
                if var_item.var.name == plot_var.name and var_item.var.file.name == plot_var.file.name:
                    var_item.var = plot_var
                    self.setRowColor(i, GREEN)

    def add_row(self, file_item, var_item):
        row = self._view.rowCount()
        self._view.setRowCount(row + 1)
        self._view.setItem(row, 0, file_item)
        self._view.setItem(row, 1, var_item)

        scaled_opt_item = QtWidgets.QTableWidgetItem()
        scaled_opt_item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
        scaled_opt_item.setCheckState(QtCore.Qt.Unchecked)
        self._view.setItem(row, 2, scaled_opt_item)

        bounds_opt_item = QtWidgets.QTableWidgetItem()
        bounds_opt_item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
        bounds_opt_item.setCheckState(QtCore.Qt.Unchecked)
        self._view.setItem(row, 3, bounds_opt_item)

        add_btn = Button(4, "Add", self._view)
        add_btn.clicked.connect(self.add_var_to_plot)
        self._view.setCellWidget(row, 4, add_btn)

        rem_btn = Button(5, "Remove", self._view)
        rem_btn.clicked.connect(self.remove_var_from_plot)
        self._view.setCellWidget(row, 5, rem_btn)

        self._view.resizeColumnsToContents()

    def clear_vars(self):
        self._view.clear()

    def add_var_to_plot(self):
        sender = self._view.sender()
        selected_var = self._view.item(sender.row, 1)
        scaled_opt = self._view.item(sender.row, 2).checkState()
        bounds_opt = self._view.item(sender.row, 3).checkState()

        self._plot_model.add_var(selected_var.var, "y", scaled_opt=scaled_opt, bounds_opt=bounds_opt)

        self._plot_model.plot()
        self._parent_model.canvas.draw()

        self.setRowColor(sender.row, GREEN)

    def remove_var_from_plot(self):
        sender = self._view.sender()
        selected_var = self._view.item(sender.row, 1)
        rem_idx = None
        for i, plot_var in enumerate(self._plot_model.y_vars):
            if selected_var.var.name == plot_var.name and selected_var.var.file.name == plot_var.file.name:
                rem_idx = i

        # Remove variable from the plot
        if rem_idx is not None:
            self._plot_model.remove_var(rem_idx)

        # Update the plot
        self._plot_model.plot()
        self._parent_model.canvas.draw()

        self.setRowColor(sender.row, WHITE)

    def setRowColor(self, row, color):
        for j in range(self._view.columnCount() - 2):
            self._view.item(row, j).setBackground(color)
