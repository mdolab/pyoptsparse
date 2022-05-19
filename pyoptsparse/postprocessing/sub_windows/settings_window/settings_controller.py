# First party modules
from pyoptsparse.postprocessing.sub_windows.settings_window.settings_model import SettingsModel
from pyoptsparse.postprocessing.utils.base_classes import Controller


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
        self._model = SettingsModel()
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

    def populate(self):
        for task, command in self._model._shortcuts():
            pass

    def toggle_dark_mode(self):
        pass
