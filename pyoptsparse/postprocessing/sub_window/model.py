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
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse import History


class HistoryFileModel:
    """Manages top-level data for the controller"""

    def __init__(self, fp: str):
        self._file = History(fileName=fp)

        self._names = self._file.getConNames() + self._file.getObjNames() + self._file.getDVNames()

        self.x_vars = {}
        self.y_vars = {}

        self.last_x: list = []
        self.last_y: list = []

    def getNames(self):
        return self._names

    def changeFile(self, fp: str):
        self._file = History(fileName=fp)
        self.clearX()
        self.clearY()

    def addX(self, name: str):
        if name not in self.x_vars:
            self.x_vars = {**self.x_vars, **self._file.getValues(names=name)}
            self.last_x.append(name)

    def addY(self, name: str):
        if name not in self.y_vars:
            self.y_vars = {**self.y_vars, **self._file.getValues(names=name)}
            self.last_y.append(name)

    def undoX(self):
        if len(self.last_x) != 0:
            self.x_vars.pop(self.last_x[-1])
            self.last_x.pop(-1)

    def undoY(self):
        if len(self.last_y) != 0:
            self.y_vars.pop(self.last_y[-1])
            self.last_y.pop(-1)

    def clearX(self):
        self.x_vars = {}
        self.last_x = []

    def clearY(self):
        self.y_vars = {}
        self.last_x = []

    def scaleX(self):
        for key in self.x_vars.keys():
            self.x_vars[key] = self._file.getValues(names=str(key), scale=True)

    def unscaleX(self):
        for key in self.x_vars.keys():
            self.x_vars[key] = self._file.getValues(names=str(key), scale=False)

    def scaleY(self):
        for key in self.y_vars.keys():
            self.y_vars[key] = self._file.getValues(names=str(key), scale=True)

    def unscaleY(self):
        for key in self.y_vars.keys():
            self.y_vars[key] = self._file.getValues(names=str(key), scale=False)

    def majorIterX(self):
        self.x_vars["major_iterations"] = self._file.getValues(names="nMajor")["nMajor"].flatten()

    def minorIterX(self):
        self.x_vars["minor_iterations"] = np.arange(
            0, len(self._file.getValues(names=self._names[0], major=False)[self._names[0]]), 1
        )
        print(self._file.getValues(names="nMinor")["nMinor"].flatten())


if __name__ == "__main__":
    model = HistoryFileModel("test.sql")
    model.minorIterX()
    print(model.x_vars)
    model.majorIterX()
    print(model.x_vars)
