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

    def __init__(self):
        self.name = None
        self.vectorized = False  # Is the variable an array of values
        self.bounds = {"upper": [], "lower": []}
        self.plot_number = 0
        self.data = {"major_iter": {"scaled": [], "unscaled": []}, "minor_iter": {"scaled": [], "unscaled": []}}


class File(object):
    """
    Data structure for holding files and setting the backend API for
    handling the data.

    Parameters
    ----------
    file_name : str
        Name of the file
    """

    def __init__(self, file_name: str):
        self.file_name = file_name
        self.backend_name = None
        self.file_reader = None

        # Try to use OpenMDAO CaseReader to open, if Operational error
        # then use pyOptSparse History class.
        # If pyOptSparse History fails, raise file format error
        try:
            self.file_reader = om.CaseReader(file_name)
            self.backend = "OpenMDAO"
        except OperationalError:
            self.file_reader = History(self.file_name)
            self.backend = "pyOptSparse"
        except Exception as e:
            message = (
                f"File '{self.file_name.split('/')[-1]}' could not be opened by the OpenMDAO case reader or the "
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

    def get_variable(self, var_name):
        if self.backend == "OpenMDAO":
            pass
        elif self.backend == "pyOptSparse":
            pass


if __name__ == "__main__":
    File("/home/lamkina/Packages/pyoptsparse/pyoptsparse/postprocessing/HISTORY.hst")
