# Standard Python modules
from pathlib import PurePath
from typing import Dict, List

# External modules
import numpy as np

# First party modules
from pyoptsparse.pyOpt_history import History


class Bounds:
    def __init__(self):
        self._lower = None
        self._upper = None

    @property
    def lower(self):
        return self._lower

    @property
    def upper(self):
        return self._upper

    @lower.setter
    def lower(self, lower: float):
        self._lower = lower

    @upper.setter
    def upper(self, upper: float):
        self._upper = upper


class Variable:
    def __init__(self, var_name: str):
        """
        Data structure for storing variables.

        Parameters
        ----------
        var_name : str
            The name of the variable.
        """
        self._name = var_name
        self._idx = 0
        self._full_name = None
        self._options = {"scale": False, "bounds": False, "major": True}
        self._bounds = Bounds()
        self._scaled_bounds = Bounds()
        self._data_major = None
        self._data_minor = None
        self._file = None
        self._scale = 1.0
        self._label = None

    def __eq__(self, other):
        return self.name == other.name and self.idx == other.idx

    @property
    def name(self):
        return self._name

    @property
    def idx(self):
        return self._idx

    @property
    def full_name(self):
        return self._full_name

    @property
    def options(self):
        return self._options

    @property
    def bounds(self):
        return self._bounds

    @property
    def scaled_bounds(self):
        return self._scaled_bounds

    @property
    def data_major(self):
        return self._data_major

    @property
    def data_minor(self):
        return self._data_minor

    @property
    def file(self):
        return self._file

    @property
    def scale(self):
        return self._scale

    @property
    def label(self):
        return self._label

    @name.setter
    def name(self, name: str):
        self._name = name

    @idx.setter
    def idx(self, idx: int):
        self._idx = idx

    @full_name.setter
    def full_name(self, full_name: str):
        self._full_name = full_name

    @options.setter
    def options(self, options: Dict):
        self._options = options

    @bounds.setter
    def bounds(self, bounds: Bounds):
        self._bounds = bounds

    @scaled_bounds.setter
    def scaled_bounds(self, scaled_bounds: Bounds):
        self._scaled_bounds = scaled_bounds

    @data_major.setter
    def data_major(self, data_major: np.ndarray):
        self._data_major = data_major

    @data_minor.setter
    def data_minor(self, data_minor: np.ndarray):
        self._data_minor = data_minor

    @file.setter
    def file(self, file):
        self._file = file

    @scale.setter
    def scale(self, scale: float):
        self._scale = scale

    @label.setter
    def label(self, label: str):
        self._label = label

    def compute_scaled_bounds(self):
        if self._bounds.upper is not None:
            self._scaled_bounds.upper = self._bounds.upper * self.scale

        if self._bounds.lower is not None:
            self._scaled_bounds.lower = self._bounds.lower * self.scale


class File:
    def __init__(self):
        """
        Data structure for storing and accessing pyoptparse history
        files.
        """
        self._name = None
        self._short_name = None
        self._reader = None
        self._dv_names = []
        self._con_names = []
        self._obj_names = []
        self._func_names = []
        self._all_names = []
        self._y_vars = {}
        self._x_vars = {}
        self._metadata = None

    @property
    def name(self):
        return self._name

    @property
    def short_name(self):
        return self._short_name

    @property
    def reader(self):
        return self._reader

    @property
    def dv_names(self):
        return self._dv_names

    @property
    def con_names(self):
        return self._con_names

    @property
    def obj_names(self):
        return self._obj_names

    @property
    def func_names(self):
        return self._func_names

    @property
    def all_names(self):
        return self._all_names

    @property
    def y_vars(self):
        return self._y_vars

    @property
    def x_vars(self):
        return self._x_vars

    @property
    def metadata(self):
        return self._metadata

    @name.setter
    def name(self, name: str):
        self._name = name
        self._short_name = PurePath(name).name

    @reader.setter
    def reader(self, reader: History):
        self._reader = reader

    @dv_names.setter
    def dv_names(self, dv_names: List[str]):
        self._dv_names = dv_names

    @con_names.setter
    def con_names(self, con_names: List[str]):
        self._con_names = con_names

    @obj_names.setter
    def obj_names(self, obj_names: List[str]):
        self._obj_names = obj_names

    @func_names.setter
    def func_names(self, func_names: List[str]):
        self._func_names = func_names

    @all_names.setter
    def all_names(self, all_names: List[str]):
        self._all_names = all_names

    @y_vars.setter
    def y_vars(self, y_vars: Dict):
        self._y_vars = y_vars

    @x_vars.setter
    def x_vars(self, x_vars: Dict):
        self._x_vars = x_vars

    @metadata.setter
    def metadata(self, metadata: Dict):
        self._metadata = metadata

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
        self.dv_names = self.reader.getDVNames()
        self.obj_names = self.reader.getObjNames()
        self.con_names = self.reader.getConNames()
        self.func_names = self.reader.getExtraFuncsNames()
        self.all_names = [k for k in self.reader.getValues().keys()]
        self.metadata = self.reader.getMetadata()
        self.x_vars = self._get_x_variables()
        self.y_vars = self._get_y_variables()

    def refresh(self):
        """
        Calls load file to refresh the data.
        """
        self.load_file(self.name)

    def _get_y_variables(self):
        y_var_names = self._get_all_y_var_names()
        variables = {}
        for name in y_var_names:
            data_major = self.reader.getValues(name, major=True, scale=False)[name]
            data_minor = self.reader.getValues(name, major=False, scale=False)[name]

            # If the major iteration data has array dimensions greater
            # than 1, the variable is vector valued and we need to
            # create a variable object for each column of the array.
            if data_major.shape[1] > 1:
                for i, cols in enumerate(zip(data_major.T, data_minor.T)):
                    var = Variable(name)
                    var.idx = i
                    var.data_major = cols[0]
                    var.data_minor = cols[1] if len(cols[0]) != len(cols[1]) else None
                    var.file = self
                    self._set_bounds(var)
                    self._set_scale(var)
                    var.compute_scaled_bounds()
                    var.full_name = f"{var.name}_{var.idx}"
                    variables[var.full_name] = var

            # If the major iteration data is one dimensional we can
            # simply create a single variable and set the data directly
            else:
                var = Variable(name)
                var.data_major = data_major
                var.data_minor = data_minor if len(data_minor) != len(data_major) else None
                var.file = self
                self._set_bounds(var)
                self._set_scale(var)
                var.compute_scaled_bounds()
                var.full_name = f"{var.name}"
                variables[var.full_name] = var

        return variables

    def _get_x_variables(self):
        x_var_names = self._get_all_x_var_names()
        variables = {}
        for name in x_var_names:
            data_major = self.reader.getValues(name, major=True, scale=False)[name]
            data_minor = self.reader.getValues(name, major=False, scale=False)[name]
            if data_major.shape[1] > 1:
                for i, cols in enumerate(zip(data_major.T, data_minor.T)):
                    var = Variable(name)
                    var.idx = i
                    var.data_major = cols[0]
                    var.data_minor = cols[1] if len(cols[0]) != len(cols[1]) else None
                    var.file = self
                    var.full_name = f"{var.name}_{var.idx}"
                    variables[var.full_name] = var
            else:
                var = Variable(name)
                var.data_major = data_major
                var.data_minor = data_minor if len(data_minor) != len(data_major) else None
                var.file = self
                var.full_name = f"{var.name}"
                variables[var.full_name] = var

        return variables

    def _set_data(self, var: Variable):
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

    def _set_bounds(self, var: Variable):
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
                var.bounds.lower = lower_bound
            if "upper" in var_info:
                upper_bound = var_info["upper"][var.idx]
                var.bounds.upper = upper_bound

    def _set_scale(self, var: Variable):
        scaled_var_names = self.dv_names + self.con_names + self.obj_names
        if var.name in scaled_var_names:
            if var.name in self.obj_names:
                var.scale = self.reader.getObjInfo(var.name)["scale"]
            elif var.name in self.con_names:
                var.scale = self.reader.getConInfo(var.name)["scale"]
            elif var.name in self.dv_names:
                var.scale = self.reader.getDVInfo(var.name)["scale"][var.idx]

    def _get_all_x_var_names(self):
        """
        Returns all the possible x-variable names listed in the filter.
        """
        x_name_filter = ["iter", "time"] + self.dv_names
        return [name for name in self.all_names if name in x_name_filter]

    def _get_all_y_var_names(self):
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
