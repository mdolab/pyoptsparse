# Standard Python modules
import os

# External modules
from PIL import Image
from PyQt5.QtWidgets import QSizePolicy, QVBoxLayout, QWidget
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# ======================================================================
# Set matplotlib backend and plt style
# ======================================================================
matplotlib.use(backend="Qt5Agg")
plt.style.use(os.path.join("pyoptsparse", "postprocessing", "assets", "nicePlotsStyle"))


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent: Figure = None):
        """
        Matplotlib canvas using the QTAgg backend.

        Parameters
        ----------
        parent : Figure, optional
            Matplotlib figure used to set the canvas
        """
        super(MplCanvas, self).__init__(parent)
        self.fig = parent

        self.addImage()  # Add the pyoptparse image

        # Set the size policy so the plot can be resized with the parent
        # view
        FigureCanvasQTAgg.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    def addImage(self):
        """
        Adds the pyoptsparse logo to the canvas as an axis.
        """
        self.img = Image.open(os.path.join("pyoptsparse", "postprocessing", "assets", "pyOptSparse_logo.png"))
        axes = self.fig.add_subplot(111)
        axes.imshow(self.img, alpha=0.5)
        axes.axis("off")


class PlotView(QWidget):
    def __init__(self, parent: QWidget = None, width: int = 10, height: int = 5, dpi: int = 100):
        """
        Constructor

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget, optional
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
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
