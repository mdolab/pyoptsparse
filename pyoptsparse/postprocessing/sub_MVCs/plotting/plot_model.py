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

    def add_var(self, var):
        """
        Adds a y-variable to the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        self.vars.append(var)

    def remove_var(self, idx: int):
        """
        Removes a y-variable from the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        self.vars.pop(idx)

    def clear_vars(self):
        """Resets the y-variables"""
        self.vars = []

    def clear_options(self):
        self.options = {"scaled": False, "bounds": False, "major_iter": True, "minor_iter": False}

    def update_axis(self, axis):
        self.axis.cla()
        self.axis = axis

    def plot(self):
        for var in self.vars:
            # Parse variable options for plot data
            if var.options["scaled"] and var.options["major_iter"]:
                key1, key2 = "major_iter", "scaled"
            elif var.options["scaled"] and var.options["minor_iter"]:
                key1, key2 = "minor_iter", "scaled"
            elif not var.options["scaled"] and var.options["major_iter"]:
                key1, key2 = "major_iter", "unscaled"
            elif not var.options["scaled"] and var.options["minor_iter"]:
                key1, key2 = "minor_iter", "unscaled"

            print(np.arange(0, len(var.data[key1][key2]), 1))
            print(var.data[key1][key2])

            # Plot variable data
            self.axis.plot(np.arange(0, len(var.data[key1][key2]), 1), var.data[key1][key2])

            # Parse variable options for bounds data
            if var.options["scaled"] and var.options["bounds"]:
                pass  # plot scalaed bounds
            elif not var.options["scaled"] and var.options["bounds"]:
                pass  # plot unscaled bounds

    def refresh(self):
        """Refresh all files and variable lists"""
        pass
