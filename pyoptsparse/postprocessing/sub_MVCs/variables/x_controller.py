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

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.sub_MVCs.widgets import VarTableWidgetItem, FileTableWidgetItem, IterSwitchWidget


class XController:
    def __init__(self, view, plot_model, parent_model):
        self._view = view
        self._parent_model = parent_model
        self._plot_model = plot_model

    def populate_vars(self):
        var = self._plot_model.x_var
        # Create a new variable widget item for the tree view
        var_item = VarTableWidgetItem(var.name)
        file_item = FileTableWidgetItem(var.file.name_short)
        var_item.var = var

        self.add_row(file_item, var_item)

    def add_row(self, file_item, var_item):
        row = self._view.rowCount()
        self._view.setRowCount(row + 1)
        self._view.setItem(row, 0, file_item)
        self._view.setItem(row, 1, var_item)

        iter_switch = IterSwitchWidget(row, self._view)
        iter_switch.clicked.connect(self.iter_switch_togg)
        iter_switch.setToolTip("Turn on for minor iterations, off for major iterations")
        self._view.setCellWidget(row, 2, iter_switch)

        self._view.setHorizontalHeaderLabels(["File", "Name", "Minor/Major"])
        self._view.resizeColumnsToContents()
        self._view.resizeRowsToContents()

    def iter_switch_togg(self):
        sender = self._view.sender()
        x_var = self._plot_model.x_var
        if sender.isChecked():
            # Adjust the iteration option for the x-variable
            x_var.options["major"] = False
            x_var.file.get_var_data(x_var)
            for y_var in self._plot_model.y_vars:
                y_var.options["major"] = False
        else:
            x_var.options["major"] = True
            x_var.file.get_var_data(x_var)
            for y_var in self._plot_model.y_vars:
                y_var.options["major"] = True

        self._plot_model.plot()
        self._parent_model.canvas.draw()

    def clear_vars(self):
        self._view.clear()
