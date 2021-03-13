# --- Python 3.8 ---
"""
Inherits the PyQt5 push button class and implements a custom button
format.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5.QtWidgets import QPushButton

# ==============================================================================
# Extension modules
# ==============================================================================


class Button(QPushButton):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.resize(150, 30)
        # self.setStyleSheet("background-color: Grey; color: white;")
