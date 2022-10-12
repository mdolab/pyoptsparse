# Standard Python modules
from typing import List

# External modules
from PyQt6.QtWidgets import QInputDialog, QLineEdit

# First party modules
from pyoptsparse.postprocessing.baseclasses import Controller
from pyoptsparse.postprocessing.sub_windows.settings_window import SettingsController, SettingsView
from pyoptsparse.postprocessing.sub_windows.tab_window import TabController, TabView


class OptViewController(Controller):
    def __init__(self):
        """
        Main OptView controller.  Communicates between the main view and
        the different plotting tabs created by the user.

        Parameters
        ----------
        view : PyQt6.QtWidgets.QWidget
            MainView instance.
        """
        super(OptViewController, self).__init__()
        self._dark_mode = False
        self._settings_controller = SettingsController()

    def addTab(self, interactive: bool = True, tab_name: str = "Home", file_names: List = []):
        """
        Adds a tab to the main view.
        """
        if interactive:
            tab_name, ok_pressed = QInputDialog.getText(
                self.view, "Enter Tab Name", "Tab Name:", QLineEdit.EchoMode.Normal, ""
            )
        tab_controller = TabController(self.view, file_names=file_names)
        tab_view = TabView(self.view, tab_controller)
        self.view.tabs.addTab(tab_view, tab_name)

    def closeTab(self, current_index: int):
        """
        Closes a tab in the main view.

        Parameters
        ----------
        current_index : int
            The index of the current tab.
        """
        self.view.tabs.removeTab(current_index)

    def configure_settings(self):
        SettingsView(parent=self.view, controller=self._settings_controller)
        self._settings_controller.populate_rc_params()
        self._settings_controller.populate_shortcuts()
