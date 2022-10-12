# External modules
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QKeySequence, QShortcut
from PyQt6.QtWidgets import QAbstractItemView, QHBoxLayout, QLabel, QTreeWidget, QVBoxLayout, QWidget

# First party modules
from pyoptsparse.postprocessing.baseclasses.controller import Controller
from pyoptsparse.postprocessing.baseclasses.view import View
from pyoptsparse.postprocessing.sub_windows.plotting.plot_view import PlotView
from pyoptsparse.postprocessing.utils.button import Button
from pyoptsparse.postprocessing.utils.switch import Switch
from pyoptsparse.postprocessing.utils.widgets import PlotList


class TabView(QWidget, View):
    def __init__(self, parent: QWidget = None, controller: Controller = None):
        """
        The view for new tabs.

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QWidget, optional
            The parent view, by default None
        controller : Controller, optional
            The tab controller, by default None
        """
        super(TabView, self).__init__(parent)
        self._controller = controller
        self._controller.view = self
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
        # Keyboard Shortcuts
        # ==============================================================
        self.add_file_action = QShortcut(QKeySequence.StandardKey.Open, self)
        self.add_plot_action = QShortcut(QKeySequence(QKeySequence.StandardKey.New), self)
        self.figure_options_action = QShortcut(QKeySequence.StandardKey.Find, self)
        self.tight_layout_action = QShortcut(QKeySequence("Ctrl+t"), self)
        self.save_figure_action = QShortcut(QKeySequence(QKeySequence.StandardKey.Save), self)
        self.figure_home_action = QShortcut(QKeySequence(QKeySequence.StandardKey.SelectStartOfDocument), self)

        self.add_file_action.activated.connect(self._controller.open_files)
        self.add_plot_action.activated.connect(self._controller.add_plot)
        self.figure_options_action.activated.connect(self._controller.mpl_figure_options)
        self.tight_layout_action.activated.connect(self._controller.mpl_tight_layout)
        self.save_figure_action.activated.connect(self._controller.mpl_save)
        self.figure_home_action.activated.connect(self._controller.mpl_home)

        # ==============================================================
        # Plot List - Left most column of Sub Layout
        # ==============================================================
        self.plot_list = PlotList(self, self._controller)
        self.plot_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.plot_list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
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
        refresh_layout.addWidget(self.auto_refresh_togg, alignment=Qt.AlignmentFlag.AlignRight)

        # --- Set the main layout ---
        self.setLayout(layout)
