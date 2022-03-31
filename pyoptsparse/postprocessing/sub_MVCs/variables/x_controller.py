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
from pyoptsparse.postprocessing.sub_MVCs.widgets import VarTableWidgetItem, FileTableWidgetItem


class XController:
    def __init__(self, view, plot_model, parent_model):
        self._view = view
        self._parent_model = parent_model
        self._plot_model = plot_model

    def populate_vars(self):
        for var in self._plot_model.x_vars:
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

    def clear_vars(self):
        self._view.clear()
