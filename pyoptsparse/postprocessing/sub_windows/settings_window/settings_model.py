# Standard Python modules
from typing import Dict

# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Model


class SettingsModel(Model):
    def __init__(self, *args, **kwargs):
        super(SettingsModel, self).__init__(*args, **kwargs)
        self._shortcuts = {}

    @property
    def shortcuts(self) -> Dict:
        return self._shortcuts()

    def edit_shortcut(self, key, val):
        self.verify_key_command(val)
        if key in self._shortcuts:
            self._shortcuts[key] = val

    def verify_key_command(self, command):
        pass
