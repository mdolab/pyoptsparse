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
        # Button Layout - Bottom Left Layout
        # ==============================================================================
        self.plot_btn = Button("Plot", self)
        self.plot_btn.clicked.connect(self._controller.plot)
        left_layout.addWidget(self.plot_btn)

        self.remove_var_btn = Button("Remove Selected", self)
        self.remove_var_btn.clicked.connect(self._controller.remove_sel_var)
        left_layout.addWidget(self.remove_var_btn)

        self.cancel_btn = Button("Cancel", self)
        self.cancel_btn.clicked.connect(self._controller.cancel)
        left_layout.addWidget(self.cancel_btn)

        # ==============================================================================
        # Y-Variables - Top Right Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.y_var_query = QtWidgets.QLineEdit(self)
        self.y_var_query.setPlaceholderText("Search...")
        self.y_var_query.textChanged.connect(self._controller.y_var_search)
        right_layout.addWidget(self.y_var_query)

        # --- Add y-vars variable list ---
        self.y_var_tree = QtWidgets.QTreeWidget(self)
        self.y_var_tree.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.y_var_tree.setColumnCount(4)
        self.y_var_tree.setHeaderLabels(["File", "Name", "Scaled", "Bounds"])
        right_layout.addWidget(self.y_var_tree)

        # ==============================================================================
        # X-Variables - Middle Right Layout
        # ==============================================================================
        self.x_var_cbox = ExtendedComboBox(self)
        self.x_var_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.x_var_cbox.setToolTip("Type to search for variables")
        self.x_var_cbox.activated.connect(self._controller.add_x_var)
        right_layout.addWidget(self.x_var_cbox)

        self.x_var_tree = QtWidgets.QTreeWidget(self)
        self.x_var_tree.setColumnCount(2)
        self.x_var_tree.setHeaderLabels(["File", "Name"])
        right_layout.addWidget(self.x_var_tree)

        # --- Anything that needs to be done upon re-opening the window ---
        self._controller.populate_files()

        # --- Set the main layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())
