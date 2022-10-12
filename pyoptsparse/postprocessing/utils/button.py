# External modules
from PyQt6.QtWidgets import QPushButton


class Button(QPushButton):
    def __init__(self, *args, **kwargs):
        """
        Inherits the PyQt6 push button class and implements a custom
        button format.
        """
        super().__init__(*args, **kwargs)
        self.resize(150, 30)
