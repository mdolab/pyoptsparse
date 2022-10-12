# Standard Python modules
import os

# External modules
from PIL import Image
from PyQt6.QtWidgets import QSizePolicy, QVBoxLayout, QWidget
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.qt_editor import figureoptions
from matplotlib.figure import Figure
import matplotlib.patheffects as patheffects
import matplotlib.pyplot as plt
import numpy as np

# First party modules
from pyoptsparse.postprocessing.baseclasses import Model, View
from pyoptsparse.postprocessing.data_structures import Variable
from pyoptsparse.postprocessing.sub_windows.mpl_figureoptions import figure_edit
from pyoptsparse.postprocessing.utils import ASSET_PATH

# ======================================================================
# Set matplotlib backend and plt style
# ======================================================================
matplotlib.use(backend="Qt5Agg")
plt.style.use(os.path.join(ASSET_PATH, "nicePlotsStyle"))
figureoptions.figure_edit = figure_edit


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
        FigureCanvasQTAgg.setSizePolicy(self, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    def addImage(self):
        """
        Adds the pyoptsparse logo to the canvas as an axis.
        """
        self.img = Image.open(os.path.join(ASSET_PATH, "pyOptSparse_logo.png"))
        axes = self.fig.add_subplot(111)
        axes.imshow(self.img, alpha=0.5)
        axes.axis("off")


class PlotView(View):
    def __init__(self, parent: QWidget = None, width: int = 10, height: int = 5, dpi: int = 100):
        """
        Constructor

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QWidget, optional
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


class PlotModel(Model):
    def __init__(self):
        """
        Model for each plot in the tab view.
        """
        super(PlotModel, self).__init__()
        self.y_vars = {}
        self.x_var = None
        self.axis = None

    def add_var(self, var: Variable, axis: str):
        """
        Adds a variable to the data model

        Parameters
        ----------
        var: Variable object
            The variable object to be added
        axis: str
            Either "x" or "y" depending on the axis to which the
            variable will be added.
        """
        if axis == "x":
            self.x_var = var
        elif axis == "y":
            if var.full_name not in self.y_vars:
                self.y_vars[var.full_name] = var

    def remove_var(self, var: Variable, axis: str):
        """
        Removes a variable from the data model

        Parameters
        ----------
        selected_var: Variable
            The variable to be removed from the model.
        idx: int
            The index of the variable to be removed
        """
        if axis == "x":
            self.x_var = None
        elif axis == "y":
            if var.full_name in self.y_vars:
                del self.y_vars[var.full_name]

    def clear_vars(self):
        """
        Clears the x and y variable data
        """
        self.x_var = None
        self.y_vars = {}

    def update_axis(self, axis):
        """
        Updates the matplotlib axis

        Parameters
        ----------
        axis : matplotlib.axis.Axis
            The subplot axis to be updated.
        """
        self.axis = axis

    def clear_axis(self):
        """
        Clears the current axis.
        """
        for artist in self.axis.lines + self.axis.collections:
            artist.remove()

    def plot(self):
        """
        Plots the x and y variables currently stored in the model.
        """
        self.clear_axis()

        for y_var in self.y_vars.values():
            if self.x_var.name == "iter":
                if self.x_var.options["major"]:
                    x_data = self.x_var.data_major
                else:
                    x_data = np.arange(0, len(y_var.data_minor), 1)
            else:
                if self.x_var.options["major"]:
                    x_data = self.x_var.data_major
                else:
                    x_data = self.x_var.data_minor

            if y_var.options["major"]:
                y_data = y_var.data_major
            else:
                y_data = y_var.data_minor

            if y_var.options["scale"]:
                y_data = y_data * y_var.scale

            self.axis.plot(
                x_data,
                y_data,
                marker=".",
                label=y_var.full_name if y_var.label is None else y_var.label,
            )

            if y_var.options["bounds"]:
                label = y_var.full_name + ":upper_bound" if y_var.label is None else y_var.label + ":upper_bound"
                if y_var.bounds.upper is not None:
                    if y_var.options["scale"]:
                        self.axis.axhline(
                            y=y_var.scaled_bounds.upper, path_effects=[patheffects.withTickedStroke()], label=label
                        )
                    else:
                        self.axis.axhline(
                            y=y_var.bounds.upper, path_effects=[patheffects.withTickedStroke()], label=label
                        )

                if y_var.bounds.lower is not None:
                    label = y_var.full_name + ":lower_bound" if y_var.label is None else y_var.label + ":lower_bound"
                    if y_var.options["scale"]:
                        self.axis.axhline(
                            y=y_var.scaled_bounds.lower,
                            path_effects=[patheffects.withTickedStroke(angle=-135)],
                            label=label,
                        )
                    else:
                        self.axis.axhline(
                            y=y_var.bounds.lower,
                            path_effects=[patheffects.withTickedStroke(angle=-135)],
                            label=label,
                        )

        if self.axis.legend_ is not None:
            self.axis.legend()

        self.axis.relim()
        self.axis.autoscale_view()
