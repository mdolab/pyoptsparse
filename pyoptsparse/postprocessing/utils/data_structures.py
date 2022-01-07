# --- Python 3.8 ---
"""
This module contains the custom data structures for storing variables
and files with different backends.

Current backend options are: 1) pyOptSparse history file, 2) OpenMDAO
case recorder file.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
import openmdao.api as om
from sqlite3 import OperationalError

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.pyOpt_history import History


class FileError(Exception):
    pass


class Variable(object):
    """Data structure for storing data related to each variable"""

    def __init__(self, var_name: str, file_idx: int, var_idx: int = 0):
        self.name = var_name
        self.file_idx = file_idx
        self.vectorized = False  # Is the variable an array of values
        self.options = {"scaled": False, "bounds": False, "major_iter": True, "minor_iter": False}
        self.bounds = {"unscaled": {"upper": None, "lower": None}, "scaled": {"upper": None, "lower": None}}
        self.data = {"major_iter": {"scaled": [], "unscaled": []}, "minor_iter": {"scaled": [], "unscaled": []}}


class File(object):
    """
    Data structure for holding files and setting the backend API for
    handling the data.
    - File is only accessed here to encapsulate the functionality


    Parameters
    ----------
    file_name : str
        Name of the file
    """

    def __init__(self, file_name: str, file_idx: int):
        self.name = file_name
        self.name_short = file_name.split("/")[-1]
        self.file_index = file_idx
        self.backend = None
        self.file_reader = None

        # Try to use OpenMDAO CaseReader to open, if Operational error
        # then use pyOptSparse History class.
        # If pyOptSparse History fails, raise file format error
        try:
            self.file_reader = om.CaseReader(self.name)
            self.backend = "OpenMDAO"
        except OperationalError:
            self.file_reader = History(self.name)
            self.backend = "pyOptSparse"
        except Exception as e:
            message = (
                f"File '{self.name.split('/')[-1]}' could not be opened by the OpenMDAO case reader or the "
                + "pyOptSparse History API"
            )
            n = len(message)
            print("+", "-" * (n), "+")
            print("|", " " * n, "|")
            print("|", message, "|")
            print("|", " " * n, "|")
            print("+", "-" * (n), "+")
            raise e

    def refresh(self):
        pass

    def get_single_variable(self, var: Variable):
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

        # --- Use the case reader for OpenMDAO ---
        if self.backend == "OpenMDAO":
            # --- Need to sort through metadata to get bounds and scaling ---
            try:  # Check for upper bound
                upper_bound = self.file_reader.problem_metadata["variables"][var.name]["upper"]
            except KeyError:
                upper_bound = None

            try:  # Check for lower bound
                lower_bound = self.file_reader.problem_metadata["variables"][var.name]["lower"]
            except KeyError:
                lower_bound = None

            try:  # Check for scaler
                scaler = self.file_reader.problem_metadata["variables"][var.name]["scaler"]
                if not scaler:  # Ensure scaler is not None
                    scaler = 1.0
            except KeyError:
                scaler = 1.0

            try:  # Check for adder
                adder = self.file_reader.problem_metadata["variables"][var.name]["adder"]
                if not adder:  # Ensure adder is not None
                    adder = 0.0
            except KeyError:
                adder = 0.0

            # --- Checks for bounds and scaling ---
            # See OpenMDAO docs for scaling formula: scaled_bound = scaler * (bound + adder)
            if upper_bound is None and lower_bound is None:  # No bounds exist
                var.bounds = None

            elif upper_bound is None and lower_bound is not None:  # Only a lower bound
                var.bounds["unscaled"]["upper"] = None
                var.bounds["scaled"]["upper"] = None

                var.bounds["unscaled"]["lower"] = lower_bound
                var.bounds["scaled"]["lower"] = scaler * (lower_bound + adder)

            elif lower_bound is None and upper_bound is not None:  # Only an upper bound
                var.bounds["unscaled"]["lower"] = None
                var.bounds["scaled"]["lower"] = None

                var.bounds["unscaled"]["upper"] = upper_bound
                var.bounds["scaled"]["upper"] = scaler * (upper_bound + adder)

            else:  # All bounds exist
                var.bounds["unscaled"]["upper"] = upper_bound
                var.bounds["scaled"]["upper"] = scaler * (upper_bound + adder)

                var.bounds["unscaled"]["lower"] = lower_bound
                var.bounds["scaled"]["lower"] = scaler * (lower_bound + adder)

            cases = self.file_reader.get_cases("driver")
            unscaled_list = []
            scaled_list = []
            for case in cases:
                unscaled_val = case.get_val(var.name)
                scaled_val = scaler * (unscaled_val + adder)
                unscaled_list.append(unscaled_val)
                scaled_list.append(scaled_val)

            var.data["minor_iter"]["unscaled"] = unscaled_list
            var.data["minor_iter"]["scaledd"] = scaled_list

        # --- Use the History API for pyOptSparse ---
        elif self.backend == "pyOptSparse":
            major_iter_vals_unscaled = self.file_reader.getValues(names=var.name, major=True)
            major_iter_vals_scaled = self.file_reader.getValues(names=var.name, major=True, scale=True)
            minor_iter_vals_unscaled = self.file_reader.getValues(names=var.name, major=False, scale=False)
            minor_iter_vals_scaled = self.file_reader.getValues(names=var.name, major=False, scale=True)

            var.data["major_iter"]["unscaled"] = major_iter_vals_unscaled[var.name].flatten()
            var.data["major_iter"]["scaled"] = major_iter_vals_scaled[var.name].flatten()
            var.data["minor_iter"]["unscaled"] = minor_iter_vals_unscaled[var.name].flatten()
            var.data["minor_iter"]["scaled"] = minor_iter_vals_scaled[var.name].flatten()

            try:
                info = self.file_reader.getDVInfo(var.name)
                try:
                    scale = float(info["scale"][0])
                except TypeError:
                    scale = 1.0
                var.bounds["unscaled"]["upper"] = info["upper"][0]
                var.bounds["unscaled"]["lower"] = info["lower"][0]
                var.bounds["scaled"]["upper"] = info["upper"][0] * scale
                var.bounds["scaled"]["lower"] = info["lower"][0] * scale
            except KeyError:
                try:
                    info = self.file_reader.getConInfo(var.name)
                    try:
                        scale = float(info["scale"][0])
                    except TypeError:
                        scale = 1.0
                    var.bounds["unscaled"]["upper"] = info["upper"][0]
                    var.bounds["unscaled"]["lower"] = info["lower"][0]
                    var.bounds["scaled"]["upper"] = info["upper"][0] * scale
                    var.bounds["scaled"]["lower"] = info["lower"][0] * scale
                except KeyError:
                    info = self.file_reader.getObjInfo(var.name)
                    var.bounds = None

        return var

    def get_all_variable_names(self):
        # --- Use the case reader for OpenMDAO ---
        if self.backend == "OpenMDAO":
            cases = self.file_reader.get_cases()
            obj_dict = cases[0].get_objectives()
            con_dict = cases[0].get_constraints()
            dv_dict = cases[0].get_design_vars()

            obj_names = [i for i in obj_dict.keys()]
            con_names = [i for i in con_dict.keys()]
            dv_names = [i for i in dv_dict.keys()]
            return obj_names + con_names + dv_names

        # --- Use the History API for pyOptSparse ---
        elif self.backend == "pyOptSparse":
            con_names = self.file_reader.getConNames()
            dv_names = self.file_reader.getDVNames()
            obj_names = self.file_reader.getObjNames()
            return con_names + dv_names + obj_names


if __name__ == "__main__":
    file_py = File("/home/lamkina/Packages/pyoptsparse/pyoptsparse/postprocessing/HISTORY.hst", 0)
    file_om = File("/home/lamkina/Packages/pyoptsparse/pyoptsparse/postprocessing/RECORDER.sql", 1)

    print(file_om.get_all_variable_names())
    print(file_py.get_all_variable_names())
    print(file_om.get_single_variable("dinc1", 1))
