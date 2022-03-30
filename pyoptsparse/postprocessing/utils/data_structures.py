# --- Python 3.8 ---
"""
This module contains the custom data structures for storing variables
and files.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from pathlib import PurePath

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.pyOpt_history import History


class Variable:
    """Data structure for storing data related to each variable"""

    def __init__(self, var_name: str):
        self.name = var_name
        self.options = {"scaled": False, "bounds": False, "major_iter": False}
        self.bounds = {"upper": None, "lower": None}
        self.values = None
        self.file = None
        self.plot_options = {}

    def setPlotOptions(self, **kwargs):
        self.plot_options = kwargs


class File:
    """
    Data structure for holding files and setting the backend API for
    handling the data.
    - File is only accessed here to encapsulate the functionality


    Parameters
    ----------
    file_name : str
        Name of the file
    """

    def __init__(self):
        self.name = None
        self.short_name = None
        self.file_reader = None
        self.dv_names = []
        self.con_names = []
        self.obj_names = []

    def load_file(self, file_name: str):
        self.file_reader = History(file_name)
        self.name = file_name
        self.name_short = PurePath(file_name).name
        self.dv_names = self.file_reader.getDVNames()
        self.obj_names = self.file_reader.getObjNames()
        self.con_names = self.file_reader.getConNames()
        self.all_names = [k for k in self.file_reader.getValues().keys()]

    def refresh(self):
        self.load_file(self.name)

    def get_var_data(self, var: Variable):
        """
        Gets all iterations (minor and major), bounds, and scaling
        of a single variable from either the OpenMDAO case reader or the
        pyOptSparse history API and stores them in a variable object.

        Parameters
        ----------
        var: Variable
            A Variable object which stores all data required for
            plotting.
        """

        major_opt = var.options["major_iter"]
        bound_opt = var.options["bounds"]
        scale_opt = var.options["scaled"]

        if var.name in self.dv_names:
            var_info = self.file_reader.getDVInfo(var.name)
        elif var.name in self.con_names:
            var_info = self.file_reader.getConInfo(var.name)
        elif var.name in self.obj_names:
            var_info = self.file_reader.getObjInfo(var.name)
        else:
            var_info = None
            raise ValueError("Variable is not part of the optimization problem.")

        # If the scale option is selected, use the scale factor.  Otherwise, use all ones.
        if isinstance(var_info["scale"], float):
            scale = np.array([var_info["scale"]]) if scale_opt else np.ones(1)
        else:
            scale = np.array(var_info["scale"]) if scale_opt else np.ones(len(var_info["scale"]))

        # The data is returned from the hist api as a dictionary where the key is
        # the variable name and the values are the
        data = self.file_reader.getValues(names=var.name, major=major_opt, scale=scale_opt)
        var.values = np.array(list(data.values()))[0]

        # Only get the bounds if the option is True and the requested
        # variable is not an objective function.
        if bound_opt and var.name not in self.obj_names:

            upper = np.zeros(len(var_info["upper"]))  # Initialize the upper bounds
            lower = np.zeros(len(var_info["lower"]))  # Initialize the lower bounds

            var.bounds["upper"] = np.where(upper != 1e30, upper * scale, None)
            var.bounds["lower"] = np.where(lower != 1e30, lower * scale, None)

    def get_all_x_var_names(self):
        x_name_filter = ["nMajor", "nMinor"] + self.dv_names
        return [name for name in self.all_names if name not in x_name_filter]

    def get_all_y_var_names(self):
        y_name_filter = (
            ["step", "optimality", "feasibility", "merit", "penalty"] + self.dv_names + self.con_names + self.obj_names
        )
        return [name for name in self.all_names if name not in y_name_filter]

    def get_metadata(self):
        return self.file_reader.getMetadata()
