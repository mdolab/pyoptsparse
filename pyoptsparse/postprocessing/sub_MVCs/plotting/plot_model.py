# External modules
import matplotlib.patheffects as patheffects

# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Model
from pyoptsparse.postprocessing.utils.data_structures import Variable


class PlotModel(Model):
    def __init__(self):
        """
        Model for each plot in the tab view.
        """
        super(PlotModel, self).__init__()
        self.y_vars = []
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
            names = [y_var.name for y_var in self.y_vars]
            if var.name not in names:
                self.y_vars.append(var)

    def remove_var(self, selected_var: Variable, axis: str):
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
            rem_idx = None
            for _i, plot_var in enumerate(self.y_vars):
                if selected_var.name == plot_var.name and selected_var.file.name == plot_var.file.name:
                    rem_idx = _i

            if rem_idx is not None:
                self.y_vars.pop(rem_idx)

    def clear_vars(self):
        """
        Clears the x and y variable data
        """
        self.x_var = None
        self.y_vars = []

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
        Steps include:
            1) Clear the axis
            2) Get the x-variable data from the history file
            3) Loop over the y-variables, get their data from the history
            file, and then add them to the plot.
            4) Check bounds option and plot bounds if True
            5) Reset the axis limits
            6) Autoscale the axis
        """
        self.clear_axis()
        if self.x_var is not None:
            self.x_var.file.get_var_data(self.x_var)
        for y_var in self.y_vars:
            # Set some default plotting options for the variable
            y_var.setPlotOptions(marker=".")
            y_var.file.get_var_data(y_var)
            if len(self.x_var.values) != len(y_var.values):
                return False
            self.axis.plot(
                self.x_var.values,
                y_var.values,
                # color=y_var.plot_options["color"],
                marker=y_var.plot_options["marker"],
                label=y_var.name,
            )

            if y_var.options["bounds"]:
                if y_var.bounds["upper"] is not None:
                    for ub in y_var.bounds["upper"]:
                        if ub is not None:
                            self.axis.axhline(
                                y=ub,
                                path_effects=[patheffects.withTickedStroke()],
                            )
                if y_var.bounds["lower"] is not None:
                    for lb in y_var.bounds["lower"]:
                        if lb is not None:
                            self.axis.axhline(
                                y=lb,
                                path_effects=[patheffects.withTickedStroke(angle=-135)],
                            )
        self.axis.relim()
        self.axis.autoscale_view()

        return True
