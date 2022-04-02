# --- Python 3.8 ---
"""
Controller for the sub view.  Interacts with the data models and
handles all user input and response functionality.  Controller can
only update the view based on user input.  If a view state is changed
which requires a messagebox view, that view is created by the controller
but managed seperately.

State control is encapsulated within it's own controller class which
is specific to this sub view.
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
from pyoptsparse.postprocessing.sub_MVCs.tab_window.tab_model import TabModel
from pyoptsparse.postprocessing.sub_MVCs.plotting.plot_model import PlotModel
from pyoptsparse.postprocessing.sub_MVCs.widgets import PlotListWidget
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_view import ConfigurePlotView
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_controller import ConfigureController
from pyoptsparse.postprocessing.sub_MVCs.metadata_window.metadata_controller import MetadataController
from pyoptsparse.postprocessing.sub_MVCs.metadata_window.metadata_view import MetadataView
from pyoptsparse.postprocessing.sub_MVCs.metadata_window.metadata_model import MetadataModel


class TabViewController:
    """
    Contains functionality for user input and software
    response for the main view.
    """

    def __init__(self, root, view, files=[]):
        self._root = root
        self._model = TabModel(files)
        self._view = view
        self._sub_views = []

    def open_files(self):
        # --- Set file dialog options ---
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog

        # --- Open file dialog and get selected user files ---
        file_names, _ = QtWidgets.QFileDialog.getOpenFileNames(self._view, "Open History File", "", "", options=options)

        # --- Load files into the model ---
        self._model.load_files(file_names)

    def add_plot(self):
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
            item = QtWidgets.QListWidgetItem(self._view.plot_list)

            # --- Create custom plot list widget ---
            plot_list_widget = PlotListWidget(self._view, self, idx)

            # --- Size the list row to fit custom widget ---
            item.setSizeHint(plot_list_widget.sizeHint())

            # --- Add the item and custom widget to the list ---
            self._view.plot_list.addItem(item)
            self._view.plot_list.setItemWidget(item, plot_list_widget)

            self.refresh_plots()

        except ValueError:
            # --- Show warning if more than 3 plots are added ---
            QtWidgets.QMessageBox.warning(self._view, "Subplot Value Warning", "OptView can only handle 3 subplots")

    def remove_plot(self, idx):
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

    def configure_view(self, idx: int):
        self._model.sub_views[idx].show()

    def auto_refresh(self):
        if self._view.auto_refresh_togg.isChecked():
            self._view.refresh_btn.setEnabled(False)
            self._model.timer.start(5000)
            self._model.timer.timeout.connect(self.refresh)
        else:
            self._model.timer.stop()
            self._view.refresh_btn.setEnabled(True)

    def refresh(self):
        for file in self._model.files:
            file.refresh()

        for plot_model in self._model.plots:
            for var in plot_model.vars:
                var.file.get_var_data(var)

        self.refresh_plots()

    def refresh_plots(self):
        for p in self._model.plots:
            p.plot()

        self._model.canvas.draw()

    def set_model_canvas(self, canvas):
        self._model.canvas = canvas

    def meta_view(self):
        meta_controller = MetadataController(MetadataModel(), self._model.files)
        MetadataView(self._root, meta_controller, "Metadata Viewer")
