# External modules
from PyQt5.QtWidgets import QWidget


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


class Model(object):
    def __init__(self, *args, **kwargs):
        """
        Base class for model objects.
        """
        pass
