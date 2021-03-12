# --- Python 3.8 ---
"""
View class for the matplotlib plotting canvas
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt5 import QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================

# --- Set matplotlib package settings ---
matplotlib.use("Qt5Agg")


class MplCanvas(FigureCanvasQTAgg):
    """Configures a canvas using the matplotlib backend for PyQT5"""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        """
        Initializer for the matplotlib canvas

        Parameters
        ----------
        parent : PyQt5 view, optional
            The view to embed the canvas, by default None
        width : int, optional
            Width of the plot canvas, by default 5
        height : int, optional
            Height of the plot canvas, by default 4
        dpi : int, optional
            Display resolution for the canvas, by default 100
        """

        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

        self.setParent(parent)

        FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)
        self.plot()

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

        self.axes.plot(x_data, y_data)
        self.draw()

    def clear(self):
        """Clears the matplotlib canvas"""

        self.axes.cla()
        self.draw()
