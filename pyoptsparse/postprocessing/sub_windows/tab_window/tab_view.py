# External modules
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QHBoxLayout, QLabel, QListWidget, QTreeWidget, QVBoxLayout, QWidget

# First party modules
from pyoptsparse.postprocessing.sub_windows.plotting.plot_view import PlotView
from pyoptsparse.postprocessing.utils.base_classes import Controller
from pyoptsparse.postprocessing.utils.button import Button
from pyoptsparse.postprocessing.utils.switch import Switch


class TabView(QWidget):
    def __init__(self, parent: QWidget = None, controller: Controller = None):
        """
        The view for new tabs.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget, optional
            The parent view, by default None
        controller : Controller, optional
            The tab controller, by default None
        """
        super(TabView, self).__init__(parent)
        self._controller = controller
        self._controller.set_view(self)
        self._initView()
        self._controller.populate_files()

    def _initView(self):
        """
        Initializes the tab view.
        """
        # --- Create top level layout ---
        layout = QVBoxLayout()

        # --- Create plot view and add to layout ---
        self.plot_view = PlotView(self)
        self._controller.set_model_canvas(self.plot_view.canvas)
        layout.addWidget(self.plot_view)

        # --- Create layout underneath the plot ---
        bottom_layout = QHBoxLayout()
        layout.addLayout(bottom_layout)

        # ==============================================================
        # Plot List - Left most column of Sub Layout
        # ==============================================================
        self.plot_list = QListWidget(self)
        bottom_layout.addWidget(self.plot_list, 3)

        # ==============================================================
        # File List - Middle column of Sub Layout
        # ==============================================================
        self.file_tree = QTreeWidget(self)
        self.file_tree.setColumnCount(1)
        self.file_tree.setHeaderLabels(["File Name"])
        bottom_layout.addWidget(self.file_tree, 1)

        # ==============================================================
        # Button Layout - Sub-layout column for buttons
        # ==============================================================
        # --- Create sublayout for buttons ---
        button_layout = QVBoxLayout()
        bottom_layout.addLayout(button_layout, 1)

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

        # --- Auto refresh file Toggle ---
        # Need to add a sub layout for the toggle switch
        refresh_layout = QHBoxLayout()
        button_layout.addLayout(refresh_layout)

        # Create and add the switch to the layout
        self.auto_refresh_togg = Switch(self)
        self.auto_refresh_togg.clicked.connect(self._controller.auto_refresh)
        self.auto_refresh_lbl = QLabel("Auto Refresh Files")
        self.auto_refresh_lbl.setBuddy(self.auto_refresh_togg)  # This attaches the label to the switch widget

        # Still need to add and align each widget even though they are set as buddies
        refresh_layout.addWidget(self.auto_refresh_lbl)
        refresh_layout.addWidget(self.auto_refresh_togg, alignment=Qt.AlignRight)

        # --- Set the main layout ---
        self.setLayout(layout)
