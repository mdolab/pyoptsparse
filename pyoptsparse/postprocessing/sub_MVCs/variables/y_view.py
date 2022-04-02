#!/usr/bin/env python
"""
@File    :   y_view.py
@Time    :   2022/03/31
@Desc    :   None
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


class YTableWidget(QtWidgets.QTableWidget):
    def __init__(self, parent=None):
        super(YTableWidget, self).__init__(parent)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.setColumnCount(6)
        self.setHorizontalHeaderLabels(["File", "Name", "Scaled", "Bounds", "Add", "Remove"])
        self._controller = None

    def setController(self, controller):
        self._controller = controller
