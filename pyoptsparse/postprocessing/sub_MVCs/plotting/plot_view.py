# --- Python 3.8 ---
"""
View class for the matplotlib plotting canvas
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from PIL import Image
import os

# ==============================================================================
# External Python modules
# ==============================================================================
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5 import QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================

# --- Set matplotlib backend settings to use Qt5 ---
matplotlib.use("Qt5Agg")
dir_path = os.path.dirname(os.path.realpath(__file__))
plt.style.use(os.path.join(dir_path, "nicePlotsStyle"))


class MplCanvas(FigureCanvasQTAgg):
    """Configures a canvas using the matplotlib backend for PyQT5"""

    def __init__(self, parent=None):
        """
        Initializer for the matplotlib canvas

        Parameters
        ----------
        parent : Figure, optional
            Matplotlib figure used to set the canvas
        """
        super(MplCanvas, self).__init__(parent)
        self.fig = parent

        self.addImage()

        FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    def addImage(self):
        """
        Adds the pyoptsparse logo to the canvas as an axis.
        """
        self.img = Image.open("assets/pyOptSparse_logo.png")
        axes = self.fig.add_subplot(111)
        axes.imshow(self.img, alpha=0.5)
        axes.axis("off")


class PlotView(QtWidgets.QWidget):
    """Creates the plot view with a canvas and matplotlib toolbar"""

    def __init__(self, parent=None, width: int = 10, height: int = 5, dpi: int = 100):
        """
        Constructor

        Parameters
        ----------
        parent : QtWidgets.QtWidget, optional
            Parent Window, by default None
        width : int, optional
            Figure Width, by default 10
        height : int, optional
            Figure Height, by default 5
        dpi : int, optional
            Figure dpi, by default 100
        """
        super(PlotView, self).__init__(parent)

        # Create three "plot" QPushButton widgets

        # Create a maptlotlib FigureCanvas object,
        # which defines a single set of axes as its axes attribute
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.canvas = MplCanvas(parent=fig)

        # Create toolbar for the figure:
        # * First argument: the canvas that the toolbar must control
        # * Second argument: the toolbar's parent (self, the PlotterWidget)
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Define and apply widget layout
        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
