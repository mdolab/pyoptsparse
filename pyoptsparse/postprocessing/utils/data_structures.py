# Standard Python modules
from pathlib import PurePath

# First party modules
from pyoptsparse.pyOpt_history import History


class Variable:
    def __init__(self, var_name: str):
        """
        Data structure for storing variables.

        Parameters
        ----------
        var_name : str
            The name of the variable.
        """
        self.name = var_name
        self.idx = 0
        self.full_name = None
        self.options = {"scale": False, "bounds": False, "major": True}
        self.bounds = {"upper": None, "lower": None}
        self.bounds_scaled = {"upper": None, "lower": None}
        self.data_major = None
        self.data_minor = None
        self.file = None
        self.scale = 1.0
        self.plot_options = {}
        self.label = None

    def __eq__(self, other):
        return self.name == other.name and self.idx == other.idx

    def set_plot_options(self, **kwargs):
        """
        Stores matplotlib options for this variable.
        """
        self.plot_options = kwargs

    def set_scaled_bounds(self):
        if self.bounds["upper"] is not None:
            self.bounds_scaled["upper"] = self.bounds["upper"] * self.scale

        if self.bounds["lower"] is not None:
            self.bounds_scaled["lower"] = self.bounds["lower"] * self.scale

    def set_full_name(self):
        self.full_name = f"{self.name}_{self.idx}"

    def set_label(self, label: str):
        self.label = label


class File:
    def __init__(self):
        """
        Data structure for storing and accessing pyoptparse history
        files.
        """
        self.name = None
        self.short_name = None
        self.reader = None
        self.dv_names = []
        self.con_names = []
        self.obj_names = []
        self.func_names = []
        self.y_vars = []
        self.x_vars = []
        self.y_var_groups = set()

    def load_file(self, file_name: str):
        """
        Loads a file using the pyoptsparse history api.

        Parameters
        ----------
        file_name : str
            The name of the file.
        """
        self.reader = History(file_name)
        self.name = file_name
        self.name_short = PurePath(file_name).name
        self.dv_names = self.reader.getDVNames()
        self.obj_names = self.reader.getObjNames()
        self.con_names = self.reader.getConNames()
        self.func_names = self.reader.getExtraFuncsNames()
        self.all_names = [k for k in self.reader.getValues().keys()]

    def load_variables(self):
        self.y_vars = self.get_y_variables()
        self.x_vars = self.get_x_variables()

    def refresh(self):
        """
        Calls load file to refresh the data.
        """
        self.load_file(self.name)
        self.load_variables()

    def get_y_variables(self):
        y_var_names = self.get_all_y_var_names()
        variables = []
        for name in y_var_names:
            self.y_var_groups.add(name)
            data_major = self.reader.getValues(name, major=True, scale=False)[name]
            data_minor = self.reader.getValues(name, major=False, scale=False)[name]
            if data_major.ndim > 1:
                for i, cols in enumerate(zip(data_major.T, data_minor.T)):
                    var = Variable(name)
                    var.idx = i
                    var.data_major = cols[0]
                    var.data_minor = cols[1] if len(cols[0]) != len(cols[1]) else None
                    var.file = self
                    self.set_bounds(var)
                    self.set_scale(var)
                    var.set_scaled_bounds()
                    var.set_full_name()
                    variables.append(var)
            else:
                var = Variable(name)
                var.data_major = cols[0]
                var.data_minor = cols[1] if len(cols[0]) != len(cols[1]) else None
                var.file = self
                self.set_bounds(var)
                self.set_scale(var)
                var.set_scaled_bounds()
                var.set_full_name()
                variables.append(var)

        return variables

    def get_x_variables(self):
        x_var_names = self.get_all_x_var_names()
        variables = []
        for name in x_var_names:
            data_major = self.reader.getValues(name, major=True, scale=False)[name]
            data_minor = self.reader.getValues(name, major=False, scale=False)[name]
            if data_major.ndim > 1:
                for i, cols in enumerate(zip(data_major.T, data_minor.T)):
                    var = Variable(name)
                    var.idx = i
                    var.data_major = cols[0]
                    var.data_minor = cols[1] if len(cols[0]) != len(cols[1]) else None
                    var.file = self
                    var.set_full_name()
                    variables.append(var)
            else:
                var = Variable(name)
                var.data_major = data_major
                var.data_minor = data_minor if len(data_minor) != len(data_major) else None
                var.file = self
                var.set_full_name()
                variables.append(var)

        return variables

    def set_data(self, var: Variable):
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

        data_major = self.reader.getValues(var.name, major=True, scale=False)[var.name]
        data_minor = self.reader.getValues(var.name, major=False, scale=False)[var.name]

        if data_major.ndim > 1:
            for i, cols in enumerate(zip(data_major.T, data_minor.T)):
                if var.idx == i:
                    var.data_major = cols[0]
                    var.data_minor = cols[1] if len(cols[0]) != len(cols[1]) else None
        else:
            var.data_major = data_major
            var.data_minor = data_minor if len(data_minor) != len(data_major) else None

    def set_bounds(self, var: Variable):
        bounded_var_names = self.dv_names + self.con_names
        var_info = None
        if var.name in bounded_var_names:
            if var.name in self.dv_names:
                var_info = self.reader.getDVInfo(var.name)
            else:
                var_info = self.reader.getConInfo(var.name)

        if var_info is not None:
            if "lower" in var_info:
                lower_bound = var_info["lower"][var.idx]
                var.bounds["lower"] = lower_bound
            if "upper" in var_info:
                upper_bound = var_info["upper"][var.idx]
                var.bounds["upper"] = upper_bound

    def set_scale(self, var: Variable):
        scaled_var_names = self.dv_names + self.con_names + self.obj_names
        if var.name in scaled_var_names:
            if var.name in self.obj_names:
                var.scale = self.reader.getObjInfo(var.name)["scale"]
            elif var.name in self.con_names:
                var.scale = self.reader.getConInfo(var.name)["scale"]
            elif var.name in self.dv_names:
                var.scale = self.reader.getDVInfo(var.name)["scale"][var.idx]

    def get_all_x_var_names(self):
        """
        Returns all the possible x-variable names listed in the filter.
        """
        x_name_filter = ["iter", "time"] + self.dv_names
        return [name for name in self.all_names if name in x_name_filter]

    def get_all_y_var_names(self):
        """
        Returns all the possible y-variable names listed in the filter
        """
        y_name_filter = (
            ["step", "time", "optimality", "feasibility", "merit", "penalty"]
            + self.dv_names
            + self.con_names
            + self.obj_names
        )
        return [name for name in self.all_names if name in y_name_filter]

    def get_metadata(self):
        """
        Returns the metadata from the history file.
        """
        return self.reader.getMetadata()
