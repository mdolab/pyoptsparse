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
from pyoptsparse.postprocessing.sub_MVCs.tab_window.tab_widgets import PlotListWidget
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_view import ConfigurePlotView
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_controller import ConfigureController


class TabViewController:
    """
    Contains functionality for user input and software
    response for the main view.
    """

    def __init__(self, root, view):
        self._root = root
        self._model = TabModel()
        self._view = view

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
                raise ValueError

            # --- Clear the plot to prepare for axis update ---
            self.clear_plot()

            # --- Update previous plots to reflect new number of axis ---
            for i, p in enumerate(self._model.plots):
                p.update_axis(self._model.canvas.fig.add_subplot(int(f"{idx+1}1{i+1}")))

            # --- Create a plot object and set its axis ---
            plot = PlotModel()
            plot.axis = self._model.canvas.fig.add_subplot(int(f"{idx+1}1{idx+1}"))
            self._model.add_plot(plot)

            # --- Create socket for custom widget ---
            item = QtWidgets.QListWidgetItem(self._view.plot_list)

            # --- Create custom widget ---
            plot_list_widget = PlotListWidget(self._view, self, idx)

            # --- Size the list row to fit custom widget ---
            item.setSizeHint(plot_list_widget.sizeHint())

            # --- Add the item and custom widget to the list ---
            self._view.plot_list.addItem(item)
            self._view.plot_list.setItemWidget(item, plot_list_widget)

            # TODO: Redraw all plots after updating axis
        except ValueError:
            # --- Show warning if more than 3 plots are added ---
            QtWidgets.QMessageBox.warning(self._view, "Subplot Value Warning", "OptView can only handle 3 subplots")

    def remove_plot(self, idx):
        # --- Remove the plot from the model ---
        self._model.remove_plot(idx)
        self._view.plot_list.takeItem(idx)

        # --- Loop over custom widgets and update the index and plot number ---
        for i in range(len(self._model.plots)):
            item = self._view.plot_list.item(i)
            widget = self._view.plot_list.itemWidget(item)
            widget.idx = i
            widget.title.setText(f"Plot {i}")

    def clear_plot(self):
        self._model.canvas.fig.clf()

    def configure_view(self, idx: int, name: str):
        configure_plot_controller = ConfigureController(self._model, self._model.plots[idx])
        ConfigurePlotView(self._root, configure_plot_controller, name)

    def auto_refresh(self):
        print("Auto Refresh")

    def set_model_canvas(self, canvas):
        self._model.canvas = canvas
