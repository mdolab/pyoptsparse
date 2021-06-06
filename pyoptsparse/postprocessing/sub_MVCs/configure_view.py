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


class ConfigurePlotView(QtWidgets.QDialog):
    def __init__(self, parent, controller, name: str):
        super(ConfigurePlotView, self).__init__(parent)
        self._center()  # centers the application in the middle of the screen
        self.setWindowTitle(name)  # sets the GUI title
        self._controller = controller
        self._controller.set_view(self)
        self._initView()

    def _initView(self):
        # --- Create top level layout ---
        layout = QtWidgets.QVBoxLayout()

        # --- Create top layout
        top_layout = QtWidgets.QHBoxLayout()
        layout.addLayout(top_layout)

        # --- Create sub layout for files ---
        top_left_layout = QtWidgets.QVBoxLayout()
        top_layout.addLayout(top_left_layout, 1)

        # --- Create sub layout for x-variables ---
        top_center_layout = QtWidgets.QVBoxLayout()
        top_layout.addLayout(top_center_layout, 2)

        # --- Create sub layout for y-variables ---
        top_right_layout = QtWidgets.QVBoxLayout()
        top_layout.addLayout(top_right_layout, 2)

        # ==============================================================================
        # File Layout - Left most column of Sub Layout
        # ==============================================================================
        self.file_list = QtWidgets.QListWidget()
        self.file_list.clicked.connect(self._controller.file_selected)
        top_left_layout.addWidget(self.file_list)
        self._controller.populate_files()  # This must be called after file list is instantiated

        # ==============================================================================
        # X Variable Layout - Left Center column of Sub Layout
        # ==============================================================================
        # --- Add x-vars combobox ---
        self.x_cbox = ExtendedComboBox(self)
        self.x_cbox.setToolTip("Type to search for x-variables")
        self.x_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.x_cbox.activated.connect(self._controller.add_x_var)
        top_center_layout.addWidget(self.x_cbox)

        # --- Add x-vars variable list ---
        self.x_list = QtWidgets.QListWidget(self)
        top_center_layout.addWidget(self.x_list)

        # ==============================================================================
        # Y Variable Layout - Right Center column of Sub Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.y_cbox = ExtendedComboBox(self)
        self.y_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.y_cbox.setToolTip("Type to search for y-variables")
        self.y_cbox.activated.connect(self._controller.add_y_var)
        top_right_layout.addWidget(self.y_cbox)

        # --- Add y-vars variable list ---
        self.y_list = QtWidgets.QListWidget(self)
        top_right_layout.addWidget(self.y_list)

        # --- Set the main layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())
