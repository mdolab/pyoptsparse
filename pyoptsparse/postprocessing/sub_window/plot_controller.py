# --- Python 3.8 ---
"""
Controller for all matplotlib canvas related input and response
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


class PlotController:
    """Handles all matplotlib input and updates the view accordingly."""

    def __init__(self, canvas):
        self.canvas = canvas

    def plot(self, x_data=[], y_data=[]):
        """
        Plot function for updating the Canvas

        Parameters
        ----------
        x_data : list, optional
            List of x data to be plotted, by default []
        y_data : list, optional
            List of y data to be plotted, by default []
        """

        self.canvas.axes.plot(x_data, y_data)
        self.canvas.draw()  # draw updates the plot

    def stackedPlot(self, x_data=[], y_data_1=[], y_data_2=[]):
        pass

    def clear(self):
        """Clears the matplotlib canvas"""

        self.canvas.axes.cla()
        self.canvas.draw()  # draw updates the plot
