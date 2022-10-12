# Standard Python modules
from typing import Dict

# First party modules
from pyoptsparse.postprocessing.baseclasses.model import Model


class SettingsModel(Model):
    def __init__(self, *args, **kwargs):
        super(SettingsModel, self).__init__(*args, **kwargs)
        self._shortcuts = {}

    @property
    def shortcuts(self) -> Dict:
        return self._shortcuts

    def add_shortcut(self, shortcut: Dict):
        self._shortcuts = {**self._shortcuts, **shortcut}
