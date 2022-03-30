#!/usr/bin/env python
"""
@File    :   metadata_view.py
@Time    :   2022/03/30
@Desc    :   View for meta
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


class MetadataView(QtWidgets.QDialog):
    def __init__(self, parent, controller, name):
        super(MetadataView, self).__init__(parent)
        self._center()
        self.setWindowTitle(name)
        self._controller = controller
        self._controller.set_view(self)
        self.resize(1000, 800)
        self._initView()

    def _initView(self):
        # --- Create layout ---
        layout = QtWidgets.QHBoxLayout()

        # ==============================================================================
        # File List - Left Layout
        # ==============================================================================
        # --- File list ---
        self.file_tree = QtWidgets.QTreeWidget(self)
        self.file_tree.setColumnCount(1)
        self.file_tree.setHeaderLabels(["File Name"])
        self.file_tree.itemClicked.connect(self._controller.file_selected)
        layout.addWidget(self.file_tree)

        # ==============================================================================
        # Options Table - Right Layout
        # ==============================================================================
        right_layout = QtWidgets.QVBoxLayout()
        layout.addLayout(right_layout, 3)

        self.opt_tree = QtWidgets.QTreeWidget(self)
        self.opt_tree.setColumnCount(2)
        self.opt_tree.setHeaderLabels(["Name", "Value"])
        right_layout.addWidget(self.opt_tree)

        self.query = QtWidgets.QLineEdit()
        self.query.setPlaceholderText("Search...")
        self.query.textChanged.connect(self._controller.search)
        right_layout.addWidget(self.query)

        self.opt_prob_table = QtWidgets.QTableWidget(self)
        right_layout.addWidget(self.opt_prob_table)

        self._controller.populate_files()

        self.setLayout(layout)

        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())
