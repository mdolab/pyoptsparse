# Standard Python modules
from typing import List

# External modules
from PyQt5.QtWidgets import QFileDialog, QInputDialog, QListWidgetItem, QMessageBox, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

# First party modules
from pyoptsparse.postprocessing.sub_windows.configure_plot_window.configure_controller import ConfigureController
from pyoptsparse.postprocessing.sub_windows.configure_plot_window.configure_view import ConfigurePlotView
from pyoptsparse.postprocessing.sub_windows.metadata_window.metadata_controller import MetadataController
from pyoptsparse.postprocessing.sub_windows.metadata_window.metadata_model import MetadataModel
from pyoptsparse.postprocessing.sub_windows.metadata_window.metadata_view import MetadataView
from pyoptsparse.postprocessing.sub_windows.plotting.plot_model import PlotModel
from pyoptsparse.postprocessing.sub_windows.tab_window.tab_model import TabModel
from pyoptsparse.postprocessing.utils.base_classes import Controller
from pyoptsparse.postprocessing.utils.widgets import FileTreeWidgetItem, PlotListWidget


class TabController(Controller):
    def __init__(self, root: QWidget, file_names: List = []):
        """
        The controller for the tab view.

        Parameters
        ----------
        root : PyQt5.QtWidgets.QWidget
            The OptView main view
        file_names : List, optional
            Names of files to be pre-loaded in to the model,
            by default []
        """
        super(TabController, self).__init__()
        self._root = root
        self._model = TabModel(file_names=file_names)
        self._view = None
        self._sub_views = []

    def set_view(self, view: QWidget):
        """
        Sets the view of the controller.

        Parameters
        ----------
        view : PyQt5.QtWidgets.QWidget
            The view associated with this controller.
        """
        self._view = view

    def open_files(self):
        """
        Opens a file dialog for the user to load files into the model.
        """
        # --- Set file dialog options ---
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        # --- Open file dialog and get selected user files ---
        file_names, _ = QFileDialog.getOpenFileNames(
            self._view, "Open History File", "", "History File (*.hst)", options=options
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