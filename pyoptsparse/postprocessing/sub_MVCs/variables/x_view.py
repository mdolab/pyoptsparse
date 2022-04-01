#!/usr/bin/env python
"""
@File    :   x_view.py
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


class XTableWidget(QtWidgets.QTableWidget):
    def __init__(self, parent=None):
        super(XTableWidget, self).__init__(parent)
        self.setColumnCount(3)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self._controller = None

    def setController(self, controller):
        self._controller = controller
