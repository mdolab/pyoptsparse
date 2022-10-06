# Standard Python modules
from typing import Dict

# External modules
from PyQt5.QtWidgets import QShortcut

# First party modules
from pyoptsparse.postprocessing.baseclasses.model import Model


class SettingsModel(Model):
    def __init__(self, *args, **kwargs):
        super(SettingsModel, self).__init__(*args, **kwargs)
        self._shortcuts = {}

    @property
    def shortcuts(self) -> Dict:
        return self._shortcuts

    def add_shortcut(self, shortcut: QShortcut, desc=str):
        self._shortcuts[shortcut.key().toString()] = desc
