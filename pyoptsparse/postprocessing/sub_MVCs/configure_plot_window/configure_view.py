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
from PyQt5 import QtWidgets, QtCore

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.utils.combo_box import ExtendedComboBox
from pyoptsparse.postprocessing.utils.button import Button
from pyoptsparse.postprocessing.utils.switch import Switch


class VariableListWidget(QtWidgets.QWidget):
    def __init__(self, var_name: str = "Default Name", idx: int = 0, axis: str = "x", parent=None, controller=None):
        super(VariableListWidget, self).__init__(parent)

        self.controller = controller

        self.idx = idx

        self.axis = axis

        self.label = QtWidgets.QLabel(var_name)

        self.remove_button = QtWidgets.QPushButton("X")
        self.remove_button.clicked.connect(self.remove)

        layout = QtWidgets.QHBoxLayout()

        layout.addWidget(self.label, 10)
        layout.addWidget(self.remove_button, 1)

        self.setLayout(layout)

    def remove(self):
        self.controller.remove_variable(self.idx, self.axis)


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

        mid_layout = QtWidgets.QHBoxLayout()
        layout.addLayout(mid_layout)

        mid_left_layout = QtWidgets.QVBoxLayout()
        mid_layout.addLayout(mid_left_layout)

        mid_center_layout = QtWidgets.QVBoxLayout()
        mid_layout.addLayout(mid_center_layout)

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
        # X Variables - Top Center Layout
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
        # Y Variables - Top Right Layout
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

        # ==============================================================================
        # Options Layout 1 - Mid Left Layout
        # ==============================================================================
        self.bounds_opt = QtWidgets.QCheckBox("Y-Variable Bounds")
        self.bounds_opt.clicked.connect(self._controller.bounds_opt_checked)
        mid_left_layout.addWidget(self.bounds_opt)

        self.scale_opt = QtWidgets.QCheckBox("Scale Y-Variables")
        self.scale_opt.clicked.connect(self._controller.scale_opt_checked)
        mid_left_layout.addWidget(self.scale_opt)

        # ==============================================================================
        # Options Layout 2 -  Mid Center Layout
        # ==============================================================================
        major_iter_layout = QtWidgets.QHBoxLayout()
        mid_center_layout.addLayout(major_iter_layout)
        self.major_iter_togg = Switch(self)
        self.major_iter_togg.clicked.connect(self._controller.major_iter_toggled)
        self.major_iter_lbl = QtWidgets.QLabel("X-Axis as Major Iterations")
        self.major_iter_lbl.setBuddy(self.major_iter_togg)
        major_iter_layout.addWidget(self.major_iter_lbl)
        major_iter_layout.addWidget(self.major_iter_togg, alignment=QtCore.Qt.AlignRight)

        minor_iter_layout = QtWidgets.QHBoxLayout()
        mid_center_layout.addLayout(minor_iter_layout)
        self.minor_iter_togg = Switch(self)
        self.minor_iter_togg.clicked.connect(self._controller.minor_iter_toggled)
        self.minor_iter_lbl = QtWidgets.QLabel("X-Axis as Minor Iterations")
        self.minor_iter_lbl.setBuddy(self.minor_iter_togg)
        minor_iter_layout.addWidget(self.minor_iter_lbl)
        minor_iter_layout.addWidget(self.minor_iter_togg, alignment=QtCore.Qt.AlignRight)

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
