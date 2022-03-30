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
from pyoptsparse.postprocessing.utils.button import Button
from pyoptsparse.postprocessing.utils.switch import Switch
from pyoptsparse.postprocessing.sub_MVCs.plotting.plot_view import PlotView
from pyoptsparse.postprocessing.sub_MVCs.tab_window.tab_controller import TabViewController


class TabView(QtWidgets.QWidget):
    def __init__(self, parent):
        super(TabView, self).__init__(parent)
        self._controller = TabViewController(parent, self)
        self._initView()

    def _initView(self):
        # --- Create top level layout ---
        layout = QtWidgets.QVBoxLayout()

        # --- Create plot view and add to layout ---
        self.plot_view = PlotView(self)
        self._controller.set_model_canvas(self.plot_view.canvas)
        layout.addWidget(self.plot_view)

        # --- Create layout underneath the plot ---
        bottom_layout = QtWidgets.QHBoxLayout()
        layout.addLayout(bottom_layout)

        # --- Create sublayout for plot list ---
        plot_list_layout = QtWidgets.QVBoxLayout()
        bottom_layout.addLayout(plot_list_layout)

        # --- Create sublayout for buttons ---
        button_layout = QtWidgets.QVBoxLayout()
        bottom_layout.addLayout(button_layout)

        # ==============================================================================
        # File Layout - Left most column of Sub Layout
        # ==============================================================================
        self.plot_list = QtWidgets.QListWidget(self)
        plot_list_layout.addWidget(self.plot_list)

        # ==============================================================================
        # Button Layout - Sub-layout column for buttons
        # ==============================================================================
        # --- Add file ---
        self.add_file_btn = Button("Add file(s)", self)
        self.add_file_btn.clicked.connect(self._controller.open_files)
        button_layout.addWidget(self.add_file_btn)

        # --- Add Plot ---
        self.plot_btn = Button("Add Plot", self)
        self.plot_btn.clicked.connect(self._controller.add_plot)
        button_layout.addWidget(self.plot_btn)

        # --- Opt Problem Metadata ---
        self.meta_btn = Button("View Metadata", self)
        self.meta_btn.clicked.connect(self._controller.meta_view)
        button_layout.addWidget(self.meta_btn)

        # --- Manually refresh history file ---
        self.refresh_btn = Button("Refresh Files", self)
        self.refresh_btn.clicked.connect(self._controller.refresh)
        button_layout.addWidget(self.refresh_btn)

        # ==============================================================================
        # Switch Layout - Sub-layout rows of button layout for toggleable options
        # ==============================================================================

        # --- Auto refresh file Toggle ---
        # Need to add a sub layout for the toggle switch
        refresh_layout = QtWidgets.QHBoxLayout()
        button_layout.addLayout(refresh_layout)

        # Create and add the switch to the layout
        self.auto_refresh_togg = Switch(self)
        self.auto_refresh_togg.clicked.connect(self._controller.auto_refresh)
        self.auto_refresh_lbl = QtWidgets.QLabel("Auto Refresh Files")
        self.auto_refresh_lbl.setBuddy(self.auto_refresh_togg)  # This attaches the label to the switch widget

        # Still need to add and align each widget even though they are set as buddies
        refresh_layout.addWidget(self.auto_refresh_lbl)
        refresh_layout.addWidget(self.auto_refresh_togg, alignment=QtCore.Qt.AlignRight)

        # --- Set the main layout ---
        self.setLayout(layout)
