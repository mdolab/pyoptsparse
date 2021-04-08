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
        self._controller.setPlotController(PlotController(self.plot.canvas))  # Need to set a controller for our plot
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

        # --- Create sublayout for LHS options ---
        opt1_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(opt1_layout)

        # --- Create sublayout for RHS options ---
        opt2_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(opt2_layout)

        # ==============================================================================
        # X Variable Layout - Left Center column of Sub Layout
        # ==============================================================================
        # --- Add x-vars combobox ---
        self.x_cbox = ExtendedComboBox(self)
        self.x_cbox.setToolTip("Type to search for x-variables")
        self.x_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.x_cbox.resize(250, 30)
        self.x_cbox.addItems(["test1", "test2", "test3"])
        self.x_cbox.currentIndexChanged.connect(self._controller.addVarX)
        x_layout.addWidget(self.x_cbox)

        # --- Add x-vars variable list ---
        self.x_label = QtWidgets.QLabel(self)
        self.x_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self.x_label.resize(250, 100)
        x_layout.addWidget(self.x_label)

        # --- Add undo x-vars button ---
        self.x_undo_btn = Button("Undo x-var", self)
        self.x_undo_btn.setToolTip("Undo add x-variable")
        self.x_undo_btn.clicked.connect(self._controller.undoVarX)
        x_layout.addWidget(self.x_undo_btn)

        # --- Add clear x-vars button ---
        self.x_clear_btn = Button("Clear x-var", self)
        self.x_clear_btn.setToolTip("Clear all x-variables")
        self.x_clear_btn.clicked.connect(self._controller.clearAllX)
        x_layout.addWidget(self.x_clear_btn)

        # ==============================================================================
        # Y Variable Layout - Right Center column of Sub Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.y_cbox = ExtendedComboBox(self)
        self.y_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.y_cbox.setToolTip("Type to search for y-variables")
        self.y_cbox.resize(250, 30)
        self.y_cbox.addItems(["test1", "test2", "test3"])
        self.y_cbox.currentIndexChanged.connect(self._controller.addVarY)
        y_layout.addWidget(self.y_cbox)

        # --- Add y-vars variable list ---
        self.y_label = QtWidgets.QLabel(self)
        self.y_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self.y_label.resize(250, 100)
        y_layout.addWidget(self.y_label)

        # --- Add undo y-vars button ---
        self.y_undo_btn = Button("Undo y-var", self)
        self.y_undo_btn.setToolTip("Undo add y-variable")
        self.y_undo_btn.clicked.connect(self._controller.undoVarY)
        y_layout.addWidget(self.y_undo_btn)

        # --- Add clear y-vars button ---
        self.y_clear_btn = Button("Clear y-var", self)
        self.y_clear_btn.setToolTip("Clear all y-variables")
        self.y_clear_btn.clicked.connect(self._controller.clearAllY)
        y_layout.addWidget(self.y_clear_btn)

        # ==============================================================================
        # Options Layout 1 - First sub-layout column for options
        # ==============================================================================
        # --- Stacked Plots ---
        self.stack_plot_opt = QtWidgets.QCheckBox("Stack plots")
        self.stack_plot_opt.clicked.connect(self._controller.stackPlots)
        opt1_layout.addWidget(self.stack_plot_opt)
        opt1_layout.setAlignment(self.stack_plot_opt, QtCore.Qt.AlignLeft)

        # --- Shared y-axis ---
        self.share_y_opt = QtWidgets.QCheckBox("Shared y-axis")
        self.share_y_opt.clicked.connect(self._controller.shareAxisY)
        opt1_layout.addWidget(self.share_y_opt)
        opt1_layout.setAlignment(self.share_y_opt, QtCore.Qt.AlignLeft)

        # --- y-axis as absolute delta values ---
        self.abs_delta_opt = QtWidgets.QCheckBox("y-axis as absolute delta values")
        self.abs_delta_opt.clicked.connect(self._controller.absDeltaY)
        opt1_layout.addWidget(self.abs_delta_opt)
        opt1_layout.setAlignment(self.abs_delta_opt, QtCore.Qt.AlignLeft)

        # --- x-axis as minor iterations ---
        self.minor_itr_opt = QtWidgets.QCheckBox("x-axis as minor iterations")
        self.minor_itr_opt.clicked.connect(self._controller.minorIterX)
        opt1_layout.addWidget(self.minor_itr_opt)
        opt1_layout.setAlignment(self.minor_itr_opt, QtCore.Qt.AlignLeft)

        # --- x-axis as major iterations ---
        self.major_itr_opt = QtWidgets.QCheckBox("x-axis as major iterations")
        self.major_itr_opt.clicked.connect(self._controller.majorIterX)
        opt1_layout.addWidget(self.major_itr_opt)
        opt1_layout.setAlignment(self.major_itr_opt, QtCore.Qt.AlignLeft)

        # --- Apply scaling factor ---
        self.scale_factor_opt = QtWidgets.QCheckBox("Apply Scaling Factor")
        self.scale_factor_opt.clicked.connect(self._controller.scaleFactor)
        opt1_layout.addWidget(self.scale_factor_opt)
        opt1_layout.setAlignment(self.scale_factor_opt, QtCore.Qt.AlignLeft)

        # ==============================================================================
        # Options Layout 2 - Second sub-layout column for options
        # ==============================================================================
        # --- Min/Max arrays ---
        self.min_max_opt = QtWidgets.QCheckBox("Min/Max for arrays")
        self.min_max_opt.clicked.connect(self._controller.minMaxPlot)
        opt2_layout.addWidget(self.min_max_opt)
        opt2_layout.setAlignment(self.min_max_opt, QtCore.Qt.AlignLeft)

        # --- Auto refresh data ---
        self.auto_refresh_opt = QtWidgets.QCheckBox("Automatically refresh history")
        self.auto_refresh_opt.clicked.connect(self._controller.autoRefresh)
        opt2_layout.addWidget(self.auto_refresh_opt)
        opt2_layout.setAlignment(self.auto_refresh_opt, QtCore.Qt.AlignLeft)

        # --- Set the main layout ---
        self.setLayout(layout)
