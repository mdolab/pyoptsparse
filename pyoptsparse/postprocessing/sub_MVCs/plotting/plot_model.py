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
        self.options = {}
        self.x_vars = []
        self.y_vars = []
        self.axis = None

    def add_x_var(self, var):
        """
        Adds an x-variable to the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        self.x_vars.append(var)

    def add_y_var(self, var):
        """
        Adds a y-variable to the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        self.y_vars.append(var)

    def remove_x_var(self, idx: int):
        """
        Removes an x-variable from the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        self.x_vars.pop(idx)

    def remove_y_var(self, idx: int):
        """
        Removes a y-variable from the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        self.y_vars.pop(idx)

    def clear_x_vars(self):
        """Resets the x-variables"""
        self.x_vars = []

    def clear_y_vars(self):
        """Resets the y-variables"""
        self.y_vars = []

    def update_axis(self, axis):
        self.axis.cla()
        self.axis = axis

    def refresh(self):
        """Refresh all files and variable lists"""
        pass
