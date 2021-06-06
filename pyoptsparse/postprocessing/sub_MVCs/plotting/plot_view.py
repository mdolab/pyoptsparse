# --- Python 3.8 ---
"""
View class for the matplotlib plotting canvas
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from PIL import Image

# ==============================================================================
# External Python modules
# ==============================================================================
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5 import QtWidgets

# ==============================================================================
# Extension modules
# ==============================================================================

# --- Set matplotlib backend settings to use Qt5 ---
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
        self.parent = parent
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(MplCanvas, self).__init__(self.fig)

        self.addImage()

        FigureCanvasQTAgg.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)
        self.setParent(parent)

    def addImage(self):
        self.img = Image.open("assets/pyOptSparse_logo.png")
        axes = self.fig.add_subplot(111)
        axes.imshow(self.img, alpha=0.5)
        axes.axis("off")


class PlotView(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(PlotView, self).__init__(parent)

        # Create three "plot" QPushButton widgets

        # Create a maptlotlib FigureCanvas object,
        # which defines a single set of axes as its axes attribute
        self.canvas = MplCanvas(self, width=10, height=5, dpi=100)

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
