# --- Python 3.8 ---
"""
Window view for each tab.
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
from pyoptsparse.postprocessing.utils.combo_box import ExtendedComboBox
from pyoptsparse.postprocessing.utils.button import Button


class ConfigurePlotView(QtWidgets.QDialog):
    def __init__(self, parent, controller, name: str):
        super(ConfigurePlotView, self).__init__(parent)
        self._center()  # centers the application in the middle of the screen
        self.setWindowTitle(name)  # sets the GUI title
        self._controller = controller
        self._controller.set_view(self)
        self.resize(1000, 800)
        self._initView()

    def _initView(self):
        # --- Create top layout
        layout = QtWidgets.QHBoxLayout()

        # --- Create sub layout for files ---
        left_layout = QtWidgets.QVBoxLayout()
        layout.addLayout(left_layout, 1)

        # --- Create sub layout for variables ---
        right_layout = QtWidgets.QVBoxLayout()
        layout.addLayout(right_layout, 2)

        # ==============================================================================
        # File Management - Top Left Layout
        # ==============================================================================
        # --- Add file(s) button ---
        self.add_file_btn = Button("Add file(s)", self)
        self.add_file_btn.clicked.connect(self._controller.add_file)
        left_layout.addWidget(self.add_file_btn)

        # --- File list ---
        self.file_tree = QtWidgets.QTreeWidget(self)
        self.file_tree.setColumnCount(1)
        self.file_tree.setHeaderLabels(["File Name"])
        self.file_tree.itemClicked.connect(self._controller.file_selected)
        left_layout.addWidget(self.file_tree)

        # ==============================================================================
        # Variables - Top Right Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.var_cbox = ExtendedComboBox(self)
        self.var_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.var_cbox.setToolTip("Type to search for variables")
        self.var_cbox.activated.connect(self._controller.add_var)
        right_layout.addWidget(self.var_cbox)

        # --- Add y-vars variable list ---
        self.var_tree = QtWidgets.QTreeWidget(self)
        self.var_tree.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.var_tree.setColumnCount(5)
        self.var_tree.setHeaderLabels(["File", "Name", "Scaled", "Bounds", "Major/Minor"])
        right_layout.addWidget(self.var_tree)

        # ==============================================================================
        # Button Layout - Mid right Layout
        # ==============================================================================
        self.plot_btn = Button("Plot", self)
        self.plot_btn.clicked.connect(self._controller.plot)
        right_layout.addWidget(self.plot_btn)

        self.remove_var_btn = Button("Remove Selected", self)
        self.remove_var_btn.clicked.connect(self._controller.remove_sel_var)
        right_layout.addWidget(self.remove_var_btn)

        self.cancel_btn = Button("Cancel", self)
        self.cancel_btn.clicked.connect(self._controller.cancel)
        right_layout.addWidget(self.cancel_btn)

        # --- Anything that needs to be done upon re-opening the window ---
        self._controller.populate_files()
        self._controller.populate_vars()

        # --- Set the main layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())
