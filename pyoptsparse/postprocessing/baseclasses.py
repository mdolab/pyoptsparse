# External modules
from PyQt6.QtWidgets import QWidget


class View(QWidget):
    def __init__(self, *args, **kwargs) -> None:
        super(View, self).__init__(*args, **kwargs)


class Controller(object):
    def __init__(self, *args, **kwargs):
        """
        Base class for controller objects.
        """
        self._view = None

    @property
    def view(self) -> View:
        return self._view

    @view.setter
    def view(self, view: View):
        self._view = view


class Model(object):
    def __init__(self, *args, **kwargs) -> None:
        """
        Base class for model objects.
        """
        pass
