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
from pyoptsparse.postprocessing.sub_MVCs.variables.y_view import YTableWidget
from pyoptsparse.postprocessing.sub_MVCs.variables.x_view import XTableWidget


class ConfigurePlotView(QtWidgets.QDialog):
    def __init__(self, parent, controller, name: str):
        super(ConfigurePlotView, self).__init__(parent)
        self._center()  # centers the application in the middle of the screen
        self.setWindowTitle(name)  # sets the GUI title
        self._controller = controller
        self._controller.set_view(self)
        self.resize(1000, 800)
        self._initView()

        # --- Anything that needs to be done upon re-opening the window ---
        self._controller.setup_var_tables()
        self._controller.populate_files()

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
        # Y-Variables - Top Right Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.y_query = QtWidgets.QLineEdit(self)
        self.y_query.setPlaceholderText("Search...")
        self.y_query.textChanged.connect(self._controller.y_var_search)
        right_layout.addWidget(self.y_query)

        # --- Add y-vars variable list ---
        self.y_table = YTableWidget(self)
        self.y_table.setColumnCount(6)
        self.y_table.setHorizontalHeaderLabels(["File", "Name", "Scaled", "Bounds", "Add", "Remove"])
        right_layout.addWidget(self.y_table)

        # ==============================================================================
        # X-Variables - Middle Right Layout
        # ==============================================================================
        self.x_cbox = ExtendedComboBox(self)
        self.x_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.x_cbox.setToolTip("Type to search for variables")
        self.x_cbox.activated.connect(self._controller.add_x_var)
        right_layout.addWidget(self.x_cbox)

        self.x_table = XTableWidget(self)
        right_layout.addWidget(self.x_table)

        # --- Set the main layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())
