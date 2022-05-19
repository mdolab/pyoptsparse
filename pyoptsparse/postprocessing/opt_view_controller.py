# Standard Python modules
from typing import List

# External modules
from PyQt5.QtWidgets import QApplication, QInputDialog, QLineEdit
import qdarkstyle

# First party modules
from pyoptsparse.postprocessing.sub_windows.settings_window.settings_controller import SettingsController
from pyoptsparse.postprocessing.sub_windows.settings_window.settings_view import SettingsView
from pyoptsparse.postprocessing.sub_windows.tab_window.tab_controller import TabController
from pyoptsparse.postprocessing.sub_windows.tab_window.tab_view import TabView
from pyoptsparse.postprocessing.utils.base_classes import Controller


class OptViewController(Controller):
    def __init__(self):
        """
        Main OptView controller.  Communicates between the main view and
        the different plotting tabs created by the user.

        Parameters
        ----------
        view : PyQt5.QtWidgets.QWidget
            MainView instance.
        """
        super(OptViewController, self).__init__()
        self._dark_mode = False

    def addTab(self, interactive: bool = True, tab_name: str = "Home", file_names: List = []):
        """
        Adds a tab to the main view.
        """
        if interactive:
            tab_name, ok_pressed = QInputDialog.getText(self.view, "Enter Tab Name", "Tab Name:", QLineEdit.Normal, "")
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

    def toggleDarkMode(self):
        app = QApplication.instance()
        if not self._dark_mode:
            self._dark_mode = True
            app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
        else:
            app.setStyleSheet("")

    def configure_settings(self):
        SettingsView(parent=self.view, controller=SettingsController())