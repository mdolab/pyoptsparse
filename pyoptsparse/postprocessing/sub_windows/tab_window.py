# Standard Python modules
from typing import List

# External modules
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QKeySequence, QShortcut
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QListWidgetItem,
    QMessageBox,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

# First party modules
from pyoptsparse.postprocessing.baseclasses import Controller, Model, View
from pyoptsparse.postprocessing.data_structures import File
from pyoptsparse.postprocessing.general_widgets import Button, FileTreeWidgetItem, PlotList, PlotListWidget, Switch
from pyoptsparse.postprocessing.sub_windows.configure_plot_window import ConfigureController, ConfigurePlotView
from pyoptsparse.postprocessing.sub_windows.metadata_window import MetadataController, MetadataModel, MetadataView
from pyoptsparse.postprocessing.sub_windows.plotting import PlotModel, PlotView


class TabModel(Model):
    def __init__(self, file_names: List = []):
        """
        The model for the tab view.

        Parameters
        ----------
        file_names : List, optional
            List of files to be pre-loaded, by default []
        """
        super(TabModel, self).__init__()
        self.canvas = None
        self.files = []
        self.plots = []
        self.sub_views = []
        self.timer = QTimer()

        if file_names:
            self.load_files(file_names)

    def load_files(self, file_names: List):
        """
        Loads files into the model.

        Parameters
        ----------
        file_names : List
            List of file names to be loaded.
        """
        curr_file_names = [file.name for file in self.files]
        for fp in file_names:
            if fp not in curr_file_names:
                file = File()
                file.load_file(fp)
                self.files.append(file)

    def add_plot(self, plot: Model, view: QDialog):
        """
        Adds a plot and the corresponding sub view to the model.

        Parameters
        ----------
        plot : Model
            The plot model being added.
        view : PyQt6.QtWidgets.QDialog
            The plot configuration sub view being added.
        """
        self.plots.append(plot)
        self.sub_views.append(view)
        self.canvas.draw()

    def remove_plot(self, idx: int):
        """
        Removes the plot and sub view at the given index.

        Parameters
        ----------
        idx : int
            The index of the plot being removed.

        Returns
        -------
        PyQt6.QtWidgets.QDialog
            The sub view associated with the plot being removed.
        """
        # --- Remove the plot object and clear the figure ---
        self.plots.pop(idx)
        view = self.sub_views.pop(idx)
        self.canvas.fig.clf()

        # --- Loop over existing plots and update the axes ---
        self.update_axes()

        self.canvas.draw()

        # --- If no plots exist then draw pyOptSparse logo ---
        if not self.plots:
            self.canvas.addImage()

        # --- Draw the canvas to show updates ---
        self.canvas.draw()

        return view

    def reorder(self, mapping):
        self.plots[:] = [self.plots[i] for i in mapping]
        self.sub_views = [self.sub_views[i] for i in mapping]

        self.update_axes()

    def swap(self, idx1, idx2):
        self.plots[idx1], self.plots[idx2] = self.plots[idx2], self.plots[idx1]
        self.sub_views[idx1], self.sub_views[idx2] = self.sub_views[idx2], self.sub_views[idx1]

        self.update_axes()

    def update_axes(self):
        num_plots = len(self.plots)
        self.canvas.fig.clf()
        for i, p in enumerate(self.plots):
            p.update_axis(self.canvas.fig.add_subplot(int(f"{num_plots}1{i+1}"), label=f"Plot {i}"))

    def draw_canvas(self):
        self.canvas.draw()


class TabController(Controller):
    def __init__(self, root: QWidget, file_names: List = []):
        """
        The controller for the tab view.

        Parameters
        ----------
        root : PyQt6.QtWidgets.QWidget
            The OptView main view
        file_names : List, optional
            Names of files to be pre-loaded in to the model,
            by default []
        """
        super(TabController, self).__init__()
        self._root = root
        self._model = TabModel(file_names=file_names)
        self._sub_views = []

    def open_files(self):
        """
        Opens a file dialog for the user to load files into the model.
        """
        # --- Open file dialog and get selected user files ---
        file_names, _ = QFileDialog.getOpenFileNames(
            self._view, "Open History File", "", "History File (*.hst)", options=QFileDialog.Option.DontUseNativeDialog
        )

        # --- Load files into the model ---
        self._model.load_files(file_names)

        self.populate_files()

    def populate_files(self):
        self._view.file_tree.clear()
        for file in self._model.files:
            file_item = FileTreeWidgetItem(self._view.file_tree)
            file_item.setFile(file)
            file_item.setText(0, file.short_name)
            self._view.file_tree.addTopLevelItem(file_item)

    def add_plot(self):
        """
        Adds a plot to the tab model.
        """
        # --- Get the number of the plot ---
        idx = len(self._model.plots)

        try:
            # --- Only allow 3 plots per tab ---
            if idx > 2:
                raise ValueError("Only 3 plots allowed per tab.")

            # --- Clear the plot to prepare for axis update ---
            self._model.canvas.fig.clf()

            # --- Update previous plots to reflect new number of axis ---
            for i, p in enumerate(self._model.plots):
                p.update_axis(self._model.canvas.fig.add_subplot(int(f"{idx+1}1{i+1}"), label=f"Plot {i}"))

            # --- Create a plot object and set its axis ---
            plot = PlotModel()
            label = f"Plot {idx}"
            plot.axis = self._model.canvas.fig.add_subplot(int(f"{idx+1}1{idx+1}"), label=label)

            # --- Create sub view and controller ---
            configure_plot_controller = ConfigureController(self._model, plot)
            sub_view = ConfigurePlotView(self._view, configure_plot_controller, label)
            self._model.add_plot(plot, sub_view)

            # --- Create socket for custom widget ---
            item = QListWidgetItem(self._view.plot_list)

            # --- Create custom plot list widget ---
            plot_list_widget = PlotListWidget(self._view, self, idx)

            # --- Size the list row to fit custom widget ---
            item.setSizeHint(plot_list_widget.sizeHint())
            item.setToolTip("Click and drag to re-order plots.")

            # --- Add the item and custom widget to the list ---
            self._view.plot_list.addItem(item)
            self._view.plot_list.setItemWidget(item, plot_list_widget)

            self.refresh_plots()

        except ValueError:
            # --- Show warning if more than 3 plots are added ---
            QMessageBox.warning(self._view, "Subplot Value Warning", "OptView can only handle 3 subplots")

    def remove_plot(self, idx: int):
        """
        Removes a plot at the given index from the model.

        Parameters
        ----------
        idx : int
            The index of the plot to be removed.
        """
        # --- Remove the plot from the model ---
        sub_view = self._model.remove_plot(idx)
        self._view.plot_list.takeItem(idx)

        if sub_view is not None:
            sub_view.close()
            sub_view.destroy()

        # --- Loop over custom widgets and update the index and plot number ---
        for i in range(len(self._model.plots)):
            item = self._view.plot_list.item(i)
            widget = self._view.plot_list.itemWidget(item)
            widget.idx = i
            widget.title.setText(f"Plot {i}")
            self._model.sub_views[i].setWindowTitle(f"Plot {i}")

        self.refresh_plots()

    def reorder_plots(self):
        mapping = []
        for i in range(self._view.plot_list.count()):
            item = self._view.plot_list.item(i)
            widget = self._view.plot_list.itemWidget(item)
            mapping.append(widget.idx)
            widget.idx = i

        self._model.reorder(mapping)
        self.refresh_plots()

    def configure_view(self, idx: int):
        """
        Opens the configuration view for the plot at the given index.

        Parameters
        ----------
        idx : int
            The index of the plot for which the configuration window
            is associated.
        """
        self._model.sub_views[idx].show()

    def auto_refresh(self):
        """
        Turns on auto refresh mode.  When activated, this function
        will refresh the history file and the plots every 5 seconds.
        """
        switch = self._view.auto_refresh_togg
        if switch.isChecked():
            time, ok = QInputDialog.getInt(self._view, "Refresh Time", "Enter refresh interval in seconds:")
            if ok:
                if time:
                    time = time * 1000
                else:
                    time = 5000

                self._model.timer.start(time)
                self._model.timer.timeout.connect(self.refresh)
            else:
                switch.setChecked(False)
        else:
            self._model.timer.stop()

    def refresh(self):
        """
        Performs a single refresh operation on the history file
        and the plots.
        """
        for file in self._model.files:
            file.refresh()

        for plot_model in self._model.plots:
            for y_var in plot_model.y_vars.values():
                y_var.file._set_data(y_var)

            x_var = plot_model.x_var
            if x_var is not None:
                x_var.file._set_data(x_var)

        self.refresh_plots()

    def refresh_plots(self):
        """
        Loops over all the plots in the model re-plots them with the
        new data from the refreshed history file.
        """
        for p in self._model.plots:
            p.plot()

        self._model.canvas.draw()

    def set_model_canvas(self, canvas: FigureCanvasQTAgg):
        """
        Sets the canvas for the model.

        Parameters
        ----------
        canvas : matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg
            The backend matplotlib canvas configured for qt5
        """
        self._model.canvas = canvas

    def meta_view(self):
        """
        Creates a meta data controller and spawns the meta data view.
        """
        meta_controller = MetadataController(MetadataModel(), self._model)
        MetadataView(self._root, meta_controller, "Metadata Viewer")

    def move_plot_up(self):
        item = self._view.plot_list.currentItem()
        row = self._view.plot_list.currentRow()
        if item is not None:
            widget = self._view.plot_list.itemWidget(item)
            if row > 0:
                self._model.swap(row, row - 1)

                # --- Create socket for custom widget ---
                new_item = item.clone()

                self._view.plot_list.insertItem(row - 1, new_item)
                self._view.plot_list.setItemWidget(new_item, widget)
                widget.idx = row - 1

                self._view.plot_list.takeItem(row + 1)
                self._view.plot_list.setCurrentRow(row - 1)

                self.refresh_plots()

    def move_plot_down(self):
        item = self._view.plot_list.currentItem()
        row = self._view.plot_list.currentRow()
        if item is not None:
            widget = self._view.plot_list.itemWidget(item)
            if row < self._view.plot_list.count() - 1:
                self._model.swap(row, row + 1)

                # --- Create socket for custom widget ---
                new_item = item.clone()

                self._view.plot_list.insertItem(row + 2, new_item)
                self._view.plot_list.setItemWidget(new_item, widget)
                widget.idx = row + 1

                self._view.plot_list.takeItem(row)
                self._view.plot_list.setCurrentRow(row + 1)

                self.refresh_plots()

    def mpl_figure_options(self):
        self._view.plot_view.toolbar.edit_parameters()

    def mpl_tight_layout(self):
        self._view.plot_view.canvas.fig.tight_layout()
        self.refresh_plots()

    def mpl_save(self):
        self._view.plot_view.toolbar.save_figure()

    def mpl_home(self):
        self._view.plot_view.toolbar.home()


class TabView(View):
    def __init__(self, parent: QWidget = None, controller: TabController = None):
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
