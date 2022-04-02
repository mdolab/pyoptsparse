# --- Python 3.8 ---
"""
Model which controls individual subplots.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
import matplotlib.patheffects as patheffects

# ==============================================================================
# Extension modules
# ==============================================================================


class PlotModel(object):
    """Manages top-level data for the controller"""

    def __init__(self):
        self.y_vars = []
        self.x_var = None
        self.axis = None

    def add_var(self, var, axis, scaled_opt=False, bounds_opt=False):
        """
        Adds a y-variable to the data model

        Parameters
        ----------
        var: Variable object
            The variable object to be added
        """
        if axis == "x":
            self.x_var = var
        elif axis == "y":
            var.options["bounds"] = bounds_opt
            var.options["scaled"] = scaled_opt
            self.y_vars.append(var)

    def remove_var(self, idx: int, axis: str):
        """
        Removes a variable from the data model

        Parameters
        ----------
        idx: int
            The index of the variable to be removed
        """
        if axis == "x":
            self.x_var = None
        elif axis == "y":
            self.y_vars.pop(idx)

    def clear_vars(self):
        """Resets the variables to an empty list"""
        self.x_var = None
        self.y_vars = []

    def update_axis(self, axis):
        self.axis = axis

    def clear_axis(self):
        for artist in self.axis.lines + self.axis.collections:
            artist.remove()
        self.axis.clear()

    def plot(self):
        self.clear_axis()
        if self.x_var is not None:
            self.x_var.file.get_var_data(self.x_var)
        for i, y_var in enumerate(self.y_vars):
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
