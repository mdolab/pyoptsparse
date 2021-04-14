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

    def plot(self, x_data=[], y_data=[], options={}):
        if options["standard"]:  # plot normal x-y
            self.canvas.fig.clf()
            ax = self.canvas.fig.add_subplot(111)
            ax.plot(x_data, y_data)
            self.canvas.draw()  # draw updates the plot

        elif options["stacked"]:  # stacked plots
            self.canvas.fig.clf()
            ax1 = self.canvas.fig.add_subplot(211)
            x1 = [i for i in range(100)]
            y1 = [i ** 0.5 for i in x1]
            ax1.set(title="Plot 1")
            ax1.plot(x1, y1, "b.-")

            ax2 = self.canvas.fig.add_subplot(212)
            x2 = [i for i in range(100)]
            y2 = [i for i in x2]
            ax2.set(title="Plot 2")
            ax2.plot(x2, y2, "b.-")
            self.canvas.draw()

        elif options == 3:  # shared X-axis plot
            self.canvas.fig.clf()
            ax1 = self.canvas.fig.add_subplot(111)
            ax1.plot(x_data, y_data)

            ax2 = ax1.twinx()
            ax2.plot(x_data, [y ** 2 for y in y_data])

            self.canvas.draw()

    def clear(self):
        """Clears the matplotlib canvas"""
        self.canvas.fig.clf()
        self.canvas.addImage()
        self.canvas.draw()  # draw updates the plot
