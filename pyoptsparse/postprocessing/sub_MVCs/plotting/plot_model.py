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
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================


class PlotModel(object):
    """Manages top-level data for the controller"""

    def __init__(self):
        self.vars = []
        self.axis = None
        self.colors = [
            ("Yellow", "#e29400ff"),
            ("Blue", "#1E90FF"),
            ("Red", "#E21A1A"),
            ("Green", "#00a650ff"),
            ("Maroon", "#800000ff"),
            ("Orange", "#ff8f00"),
            ("Purple", "#800080ff"),
            ("Cyan", "#00A6D6"),
            ("Black", "#000000ff"),
            ("Grey", "#5a5758ff"),
        ]

    def add_var(self, var):
        """
        Adds a y-variable to the data model

        Parameters
        ----------
        var: Variable object
            The variable object to be added
        """
        self.vars.append(var)

    def remove_var(self, idx: int):
        """
        Removes a variable from the data model

        Parameters
        ----------
        idx: int
            The index of the variable to be removed
        """
        self.vars.pop(idx)

    def clear_vars(self):
        """Resets the variables to an empty list"""
        self.vars = []

    def clear_options(self):
        self.options = {"scaled": False, "bounds": False, "major_iter": True}

    def update_axis(self, axis):
        self.axis = axis

    def clear_axis(self):
        for artist in self.axis.lines + self.axis.collections:
            artist.remove()

    def plot(self):
        self.clear_axis()
        for i, var in enumerate(self.vars):
            # Set some default plotting options for the variable
            var.setPlotOptions(color=self.colors[i][1], marker=".")
            var.file.get_var_data(var)
            # We need to check for multi-dimensional input
            # If multi-dim, the iters will be the number of rows in the array, a.k.a the first dimension
            iters = np.arange(0, var.values.size, 1) if var.values.ndim == 1 else np.arange(0, var.values.shape[0], 1)
            self.axis.plot(
                iters, var.values, color=var.plot_options["color"], marker=var.plot_options["marker"], label=var.name
            )

            if var.options["bounds"]:
                if var.bounds["upper"] is not None:
                    for ub in var.bounds["upper"]:
                        if ub is not None:
                            self.axis.axhline(y=ub, linestyle="--", color=var.plot_options["color"])
                if var.bounds["lower"] is not None:
                    for lb in var.bounds["lower"]:
                        if lb is not None:
                            self.axis.axhline(y=lb, linestyle="--", color=var.plot_options["color"])
