# --- Python 3.8 ---
"""
Main window for each OptView tab.
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
from .utils.combo_box import ExtendedComboBox
from .utils.button import Button
from .utils.switch import Switch
from .plot_view import PlotView
from .controller import SubWindowController
from .plot_controller import PlotController


class SubWindowView(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._controller = SubWindowController(self)
        self._initView()

    def _initView(self):
        # --- Create top level layout ---
        layout = QtWidgets.QVBoxLayout()

        # --- Create plot view and add to layout ---
        self.plot = PlotView(self)
        self._controller.set_plot_controller(PlotController(self.plot.canvas))  # Need to set a controller for our plot
        layout.addWidget(self.plot)

        # --- Create sublayout underneath the plot for buttons, forms, and options ---
        sub_layout = QtWidgets.QHBoxLayout()
        layout.addLayout(sub_layout)

        # --- Create sublayout for x-variables ---
        x_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(x_layout)

        # --- Create sublayout for y-variables ---
        y_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(y_layout)

        # --- Create sublayout for options ---
        opt1_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(opt1_layout)

        # --- Create sublayout for buttons ---
        button_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(button_layout)

        # ==============================================================================
        # X Variable Layout - Left Center column of Sub Layout
        # ==============================================================================
        # --- Add x-vars combobox ---
        self.x_cbox = ExtendedComboBox(self)
        self.x_cbox.setToolTip("Type to search for x-variables")
        self.x_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.x_cbox.resize(250, 30)
        self.x_cbox.activated.connect(self._controller.add_x_var)
        x_layout.addWidget(self.x_cbox)

        # --- Add x-vars variable list ---
        self.x_label = QtWidgets.QLabel(self)
        self.x_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self.x_label.resize(250, 100)
        x_layout.addWidget(self.x_label)

        # --- Add clear x-vars button ---
        self.x_clear_btn = Button("Clear x-var", self)
        self.x_clear_btn.setToolTip("Clear all x-variables")
        self.x_clear_btn.clicked.connect(self._controller.clear_x)
        x_layout.addWidget(self.x_clear_btn)

        # ==============================================================================
        # Y Variable Layout - Right Center column of Sub Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.y_cbox = ExtendedComboBox(self)
        self.y_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.y_cbox.setToolTip("Type to search for y-variables")
        self.y_cbox.resize(250, 30)
        self.y_cbox.activated.connect(self._controller.add_y_var)
        y_layout.addWidget(self.y_cbox)

        # --- Add y-vars variable list ---
        self.y_label = QtWidgets.QLabel(self)
        self.y_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self.y_label.resize(250, 100)
        y_layout.addWidget(self.y_label)

        # --- Add clear y-vars button ---
        self.y_clear_btn = Button("Clear y-var", self)
        self.y_clear_btn.setToolTip("Clear all y-variables")
        self.y_clear_btn.clicked.connect(self._controller.clear_y)
        y_layout.addWidget(self.y_clear_btn)

        # ==============================================================================
        # Options Layout - Sub-layout column for options
        # ==============================================================================
        # --- Stacked Plots ---
        self.stack_plot_opt = QtWidgets.QCheckBox("Stack plots")
        opt1_layout.addWidget(self.stack_plot_opt, QtCore.Qt.AlignLeft)

        # --- y-axis as absolute delta values ---
        self.abs_delta_opt = QtWidgets.QCheckBox("y-axis as absolute delta values")
        opt1_layout.addWidget(self.abs_delta_opt, QtCore.Qt.AlignLeft)

        # --- x-axis as minor iterations ---
        self.minor_itr_opt = QtWidgets.QCheckBox("x-axis as minor iterations")
        opt1_layout.addWidget(self.minor_itr_opt, QtCore.Qt.AlignLeft)

        # --- x-axis as major iterations ---
        self.major_itr_opt = QtWidgets.QCheckBox("x-axis as major iterations")
        opt1_layout.addWidget(self.major_itr_opt, QtCore.Qt.AlignLeft)

        # --- Min/Max arrays ---
        self.min_max_opt = QtWidgets.QCheckBox("Min/Max for arrays")
        opt1_layout.addWidget(self.min_max_opt, QtCore.Qt.AlignLeft)

        # --- Variable bounds and constraints ---
        self.bound_opt = QtWidgets.QCheckBox("Show variable/function bounds")
        opt1_layout.addWidget(self.bound_opt, QtCore.Qt.AlignLeft)

        # ==============================================================================
        # Button Layout - Sub-layout column for buttons
        # ==============================================================================
        # --- Add file ---
        self.add_file_btn = Button("Add file", self)
        self.add_file_btn.clicked.connect(self._controller.open_file)
        button_layout.addWidget(self.add_file_btn)

        # --- Refresh history file ---
        self.refresh_btn = Button("Refresh History File", self)
        self.refresh_btn.clicked.connect(self._controller.refresh_file)
        button_layout.addWidget(self.refresh_btn)

        # --- Plot ---
        self.plot_btn = Button("Plot", self)
        self.plot_btn.clicked.connect(self._controller.plot)
        button_layout.addWidget(self.plot_btn)

        # --- Clear Plot ---
        self.clear_plot_btn = Button("Clear Plot", self)
        self.clear_plot_btn.clicked.connect(self._controller.clear_plot)
        button_layout.addWidget(self.clear_plot_btn)

        # ==============================================================================
        # Switch Layout - Sub-layout rows of button layout for toggleable options
        # ==============================================================================
        # --- create the scale variables toggle ---
        scale_layout = QtWidgets.QHBoxLayout()
        button_layout.addLayout(scale_layout)
        self.scale_var_togg = Switch(self)
        self.scale_var_togg.clicked.connect(self._controller.scale_vars)
        self.scale_var_lbl = QtWidgets.QLabel("Apply Scaling Factor")
        self.scale_var_lbl.setBuddy(self.scale_var_togg)
        scale_layout.addWidget(self.scale_var_lbl)
        scale_layout.addWidget(self.scale_var_togg, alignment=QtCore.Qt.AlignRight)

        # --- create the auto-refresh toggle ---
        refresh_layout = QtWidgets.QHBoxLayout()
        button_layout.addLayout(refresh_layout)
        self.auto_refresh_togg = Switch(self)
        self.auto_refresh_lbl = QtWidgets.QLabel("Auto Refresh History File")
        self.auto_refresh_lbl.setBuddy(self.auto_refresh_togg)
        refresh_layout.addWidget(self.auto_refresh_lbl)
        refresh_layout.addWidget(self.auto_refresh_togg, alignment=QtCore.Qt.AlignRight)

        # --- Set the main layout ---
        self.setLayout(layout)
