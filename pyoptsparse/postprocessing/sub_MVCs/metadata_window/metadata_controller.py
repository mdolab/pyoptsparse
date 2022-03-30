#!/usr/bin/env python
"""
@File    :   metadata_controller.py
@Time    :   2022/03/30
@Desc    :   None
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5.QtCore import Qt

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.sub_MVCs.widgets import FileTreeWidgetItem, OptTreeWidgetItem, OptTableWidgetItem


class MetadataController:
    def __init__(self, model, files):
        self._view = None
        self._model = model
        self.files = files
        self._current_file = None

    def set_view(self, view):
        self._view = view

    def populate_files(self):
        for file in self.files:
            file_item = FileTreeWidgetItem(self._view.file_tree)
            file_item.setFile(file)
            file_item.setText(0, file.name_short)
            self._view.file_tree.addTopLevelItem(file_item)

        if len(self.files) > 0:
            self._current_file = self.files[0]
            self.populate_opts()

    def file_selected(self, item, column):
        self._current_file = item.file
        self.populate_opts()

    def search(self, s):
        items = self._view.opt_prob_table.findItems(s, Qt.MatchContains)
        if items:
            item = items[0]
            self._view.opt_prob_table.setCurrentItem(item)

    def populate_opts(self):
        self.clear_opts()

        metadata = self._current_file.get_metadata()
        for key, val in metadata.items():
            if key != "optOptions":
                item = OptTreeWidgetItem(self._view.opt_tree)
                item.setText(0, key)
                item.setText(1, f"{val}")
                self._view.opt_tree.addTopLevelItem(item)

        self._view.opt_prob_table.setRowCount(len(metadata["optOptions"].keys()))
        self._view.opt_prob_table.setColumnCount(2)
        for i, (key, val) in enumerate(metadata["optOptions"].items()):
            option = OptTableWidgetItem(key)
            value = OptTableWidgetItem(f"{val}")

            self._view.opt_prob_table.setItem(i, 0, option)
            self._view.opt_prob_table.setItem(i, 1, value)
        self._view.opt_prob_table.resizeColumnsToContents()
        self._view.opt_prob_table.setHorizontalHeaderLabels(["Option", "Value"])
        self._view.opt_prob_table.verticalHeader().setVisible(False)

    def clear_opts(self):
        # Clear everything
        self._view.opt_tree.clear()
        self._view.opt_prob_table.clear()
