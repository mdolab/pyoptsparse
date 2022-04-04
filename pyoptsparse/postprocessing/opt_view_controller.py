# Standard Python modules
from typing import List

# External modules
from PyQt5.QtWidgets import QInputDialog, QLineEdit, QWidget

# First party modules
from pyoptsparse.postprocessing.sub_MVCs.tab_window.tab_controller import TabController
from pyoptsparse.postprocessing.sub_MVCs.tab_window.tab_view import TabView


class OptViewController:
    def __init__(self, view: QWidget):
        """
        Main OptView controller.  Communicates between the main view and
        the different plotting tabs created by the user.

        Parameters
        ----------
        view : PyQt5.QtWidgets.QWidget
            MainView instance.
        """
        self._view = view

    def addTab(self, interactive: bool = True, tab_name: str = "Home", file_names: List = []):
        """
        Adds a tab to the main view.
        """
        if interactive:
            tab_name, ok_pressed = QInputDialog.getText(self._view, "Enter Tab Name", "Tab Name:", QLineEdit.Normal, "")
        tab_controller = TabController(self._view, file_names=file_names)
        tab_view = TabView(self._view, tab_controller)
        self._view.tabs.addTab(tab_view, tab_name)

    def closeTab(self, current_index: int):
        """
        Closes a tab in the main view.

        Parameters
        ----------
        current_index : int
            The index of the current tab.
        """
        self._view.tabs.removeTab(current_index)
