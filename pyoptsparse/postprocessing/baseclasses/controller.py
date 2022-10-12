# External modules
from PyQt6.QtGui import QShortcut
from PyQt6.QtWidgets import QWidget

# First party modules
from pyoptsparse.postprocessing.sub_windows.settings_window.settings_model import SettingsModel


class Controller(object):
    def __init__(self, *args, **kwargs):
        """
        Base class for controller objects.
        """
        self._view = None

    @property
    def view(self):
        return self._view

    @view.setter
    def view(self, view: QWidget):
        self._view = view
