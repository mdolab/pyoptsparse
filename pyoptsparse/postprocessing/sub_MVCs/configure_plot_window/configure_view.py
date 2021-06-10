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

        # --- Create sub layout for y-variables ---
        top_right_layout = QtWidgets.QVBoxLayout()
        top_layout.addLayout(top_right_layout, 2)

        mid_layout = QtWidgets.QHBoxLayout()
        layout.addLayout(mid_layout)

        mid_left_layout = QtWidgets.QVBoxLayout()
        mid_layout.addLayout(mid_left_layout)

        mid_right_layout = QtWidgets.QVBoxLayout()
        mid_layout.addLayout(mid_right_layout)

        # ==============================================================================
        # File List - Top Right Layout
        # ==============================================================================
        self.file_list = QtWidgets.QListWidget()
        self.file_list.clicked.connect(self._controller.file_selected)
        top_left_layout.addWidget(self.file_list)
        self._controller.populate_files()  # This must be called after file list is instantiated

        # ==============================================================================
        # Variables - Top Right Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.var_cbox = ExtendedComboBox(self)
        self.var_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.var_cbox.setToolTip("Type to search for y-variables")
        self.var_cbox.activated.connect(self._controller.add_var)
        top_right_layout.addWidget(self.var_cbox)

        # --- Add y-vars variable list ---
        self.var_list = QtWidgets.QListWidget(self)
        top_right_layout.addWidget(self.var_list)

        # ==============================================================================
        # Options Layout 1 - Mid Left Layout
        # ==============================================================================
        self.all_vars_list = QtWidgets.QListWidget(self)
        mid_left_layout.addWidget(self.all_vars_list)

        # ==============================================================================
        # Button Layout - Mid right Layout
        # ==============================================================================
        self.ok_btn = Button("Ok", self)
        mid_right_layout.addWidget(self.ok_btn)

        self.cancel_btn = Button("Cancel", self)
        self.cancel_btn.clicked.connect(self._controller.cancel)
        mid_right_layout.addWidget(self.cancel_btn)

        # --- Set the main layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())
