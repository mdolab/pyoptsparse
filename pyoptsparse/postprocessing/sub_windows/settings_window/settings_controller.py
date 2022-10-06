# Standard Python modules
import os

# External modules
import pkg_resources

# First party modules
from pyoptsparse.postprocessing.baseclasses.controller import Controller

ASSET_PATH = pkg_resources.resource_filename("pyoptsparse", os.path.join("postprocessing", "assets"))


class SettingsController(Controller):
    def __init__(self):
        """
        The controller for the tab view.

        Parameters
        ----------
        root : PyQt5.QtWidgets.QWidget
            The OptView main view
        file_names : List, optional
            Names of files to be pre-loaded in to the model,
            by default []
        """
        super(SettingsController, self).__init__()
        self._appearance_view = None
        self._keyshort_view = None
        self._rc_param_view = None

    @property
    def appearance_view(self):
        return self._appearance_view

    @property
    def keyshort_view(self):
        return self._keyshort_view

    @property
    def rc_param_view(self):
        return self._rc_param_view

    @appearance_view.setter
    def appearance_view(self, view):
        self._appearance_view = view

    @keyshort_view.setter
    def keyshort_view(self, view):
        self._keyshort_view = view

    @rc_param_view.setter
    def rc_param_view(self, view):
        self._rc_param_view = view

    def populate_rc_params(self):
        with open(os.path.join(ASSET_PATH, "nicePlotsStyle"), "r") as file:
            for line in file:
                if line != "\n":
                    self._rc_param_view.param_list.addItem(line)

    def populate_shortcuts(self):
        print(self._settings_model.shortcuts)
        for key, val in self._settings_model.shortcuts.items():
            self._keyshort_view.shortcut_list.addItem(f"{key}: {val}")

    def toggle_dark_mode(self):
        pass
