# --- Python 3.8 ---
"""
Data structure and management class for history files.  The controller
only has access to top level data for plotting.  Data manipulation
only occurs here and not in the controller or views.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5 import QtCore

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.utils.data_structures import File


class TabModel(object):
    """Manages top-level data for the controller"""

    def __init__(self, file_names=[]):
        self.canvas = None
        self.files = []
        self.plots = []
        self.sub_views = []
        self.timer = QtCore.QTimer()

        if file_names:
            self.load_files(file_names)

    def refresh(self):
        """Refresh all files and variable lists"""
        pass

    def load_files(self, file_names: list):
        for fp in file_names:
            file = File()
            file.load_file(fp)
            self.files.append(file)

    def add_plot(self, plot, view):
        self.plots.append(plot)
        self.sub_views.append(view)
        self.canvas.draw()

    def remove_plot(self, idx):
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
