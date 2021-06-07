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

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.utils.data_structures import File


class TabModel(object):
    """Manages top-level data for the controller"""

    def __init__(self):
        self.canvas = None
        self.files = []
        self.plots = []

    def refresh(self):
        """Refresh all files and variable lists"""
        pass

    def load_files(self, file_names: list):
        for i, fp in enumerate(file_names):
            self.files.append(File(fp, i))

    def add_plot(self, plot):
        self.plots.append(plot)
        self.canvas.draw()

    def remove_plot(self, idx):
        # --- Remove the plot object and clear the figure ---
        self.plots.pop(idx)
        self.canvas.fig.clf()

        # --- Loop over existing plots and update the axes ---
        n_plots = len(self.plots)
        for i, p in enumerate(self.plots):
            p.update_axis(self.canvas.fig.add_subplot(int(f"{n_plots}1{i+1}")))

        # TODO: Replot existing plots

        # --- If no plots exist then draw pyOptSparse logo ---
        if not self.plots:
            self.canvas.addImage()

        # --- Draw the canvas to show updates ---
        self.canvas.draw()
