# External modules
import matplotlib.patheffects as patheffects
import numpy as np

# First party modules
from pyoptsparse.postprocessing.baseclasses.model import Model
from pyoptsparse.postprocessing.utils.data_structures import Variable


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
        self.axis.clear()

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

        self.axis.relim()
        self.axis.autoscale_view()
