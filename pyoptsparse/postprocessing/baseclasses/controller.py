# External modules
from PyQt5.QtWidgets import QShortcut, QWidget

# First party modules
from pyoptsparse.postprocessing.sub_windows.settings_window.settings_model import SettingsModel


class Controller(object):
    def __init__(self, *args, **kwargs):
        """
        Base class for controller objects.
        """
        self._view = None

        # We add the settings model to the base class because all
        # controllers need access to the same settings for the entire
        # program.
        self._settings_model = SettingsModel()

    @property
    def view(self):
        return self._view

    @view.setter
    def view(self, view: QWidget):
        self._view = view

    def add_shortcut(self, shortcut: QShortcut, desc: str):
        self._settings_model.add_shortcut(shortcut, desc)
