# --- Python 3.8 ---
"""
Data structure and management class for history files.  The controller
only has access to top level data for plotting.  Data manipulation
only occurs here and not in the controller or views.
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
from .data_structures import Variable, File


class Model(object):
    """Manages top-level data for the controller"""

    def __init__(self):
        self.x_vars = []
        self.y_vars = []
        self.files = []

    def add_x_var(self, var_name: str):
        """
        Adds an x-variable to the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        var = Variable()
        self.x_vars.append(var)

    def add_y_var(self, var_name: str):
        """
        Adds a y-variable to the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        var = Variable()
        self.y_vars.append(var)

    def remove_x_var(self, var_name: str):
        """
        Removes an x-variable from the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        for i, var in enumerate(self.x_vars):
            if var.name == var_name:
                self.x_vars.pop(i)

    def remove_y_var(self, var_name: str):
        """
        Removes a y-variable from the data model

        Parameters
        ----------
        var_name : str
            Name of the variable
        """
        for i, var in enumerate(self.y_vars):
            if var.name == var_name:
                self.y_vars.pop(i)

    def clear_x_vars(self):
        """Resets the x-variables"""
        self.x_vars = []

    def clear_y_vars(self):
        """Resets the y-variables"""
        self.y_vars = []

    def refresh(self):
        """Refresh all files and variable lists"""
        pass

    def load_files(self, file_names: list):
        for fp in file_names:
            self.files.append(File(fp))
