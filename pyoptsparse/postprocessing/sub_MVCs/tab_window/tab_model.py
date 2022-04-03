# Standard Python modules
from typing import List

# External modules
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QDialog

# Local modules
from pyoptsparse.postprocessing.utils.data_structures import File
from pyoptsparse.postprocessing.utils.base_classes import Model


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
        view : PyQt5.QtWidgets.QDialog
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
        PyQt5.QtWidgets.QDialog
            The sub view associated with the plot being removed.
        """
        # --- Remove the plot object and clear the figure ---
        self.plots.pop(idx)
        view = self.sub_views.pop(idx)
        self.canvas.fig.clf()

        # --- Loop over existing plots and update the axes ---
        n_plots = len(self.plots)
        for i, p in enumerate(self.plots):
            p.update_axis(self.canvas.fig.add_subplot(int(f"{n_plots}1{i+1}")))

        self.canvas.draw()

        # --- If no plots exist then draw pyOptSparse logo ---
        if not self.plots:
            self.canvas.addImage()

        # --- Draw the canvas to show updates ---
        self.canvas.draw()

        return view