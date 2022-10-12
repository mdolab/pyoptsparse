# Standard Python modules
import os

# External modules
from PyQt6.QtGui import QKeySequence
from PyQt6.QtWidgets import QTableWidgetItem
import pkg_resources

# First party modules
from pyoptsparse.postprocessing.baseclasses.controller import Controller

# Local modules
from .settings_model import SettingsModel

ASSET_PATH = pkg_resources.resource_filename("pyoptsparse", os.path.join("postprocessing", "assets"))


class SettingsController(Controller):
    def __init__(self):
        """
        The controller for the tab view.

        Parameters
        ----------
        root : PyQt6.QtWidgets.QWidget
            The OptView main view
        file_names : List, optional
            Names of files to be pre-loaded in to the model,
            by default []
        """
        super(SettingsController, self).__init__()
        self._appearance_view = None
        self._keyshort_view = None
        self._rc_param_view = None
        self._model = SettingsModel()

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

    def populate_rc_params(self):
        with open(os.path.join(ASSET_PATH, "nicePlotsStyle"), "r") as file:
            for line in file:
                if line != "\n":
                    current_row_count = self._rc_param_view.rc_table.rowCount()
                    self._rc_param_view.rc_table.insertRow(current_row_count)
                    param_item = QTableWidgetItem(line.split(":")[0])
                    val_item = QTableWidgetItem(line.split(":")[1].strip("\n"))
                    self._rc_param_view.rc_table.setItem(current_row_count, 0, param_item)
                    self._rc_param_view.rc_table.setItem(current_row_count, 1, val_item)

    def populate_shortcuts(self):
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.Open).toString(), "Opens the add file menu.")
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.New).toString(), "Add a subplot to the figure.")
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.Find).toString(), "Opens the figure options menu.")
        self._add_shortcut("Ctrl+T", "Applies the tight layout format to the figure.")
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.Save).toString(), "Opens the save figure menu.")
        self._add_shortcut(
            QKeySequence(QKeySequence.StandardKey.SelectStartOfDocument).toString(),
            "Resets the figure to the default home view.",
        )
        self._add_shortcut("Ctrl+Up", "Moves a subplot up.")
        self._add_shortcut("Ctrl+Down", "Moves a subplot down.")
        self._add_shortcut(QKeySequence("Alt+h").toString(), "Opens the help/settings menu")

    def _add_shortcut(self, key, val):
        self._model.add_shortcut({key: val})
        current_row_count = self._keyshort_view.shortcut_table.rowCount()
        self._keyshort_view.shortcut_table.insertRow(current_row_count)
        key_item = QTableWidgetItem(key)
        val_item = QTableWidgetItem(val)
        self._keyshort_view.shortcut_table.setItem(current_row_count, 0, key_item)
        self._keyshort_view.shortcut_table.setItem(current_row_count, 1, val_item)
