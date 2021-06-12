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
        if self.options["major_iter"]:
            if self.options["scaled"]:
                for y_var in self.y_vars:
                    self.axis.plot(self.x_vars[0].data["major_iter"]["scaled"],)

    def refresh(self):
        """Refresh all files and variable lists"""
        pass
