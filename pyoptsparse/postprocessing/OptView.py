# Standard Python modules
import argparse
import os
import sys
from typing import List

# External modules
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction, QApplication, QMenuBar, QTabWidget, QVBoxLayout, QWidget, qApp

# First party modules
from pyoptsparse.postprocessing.baseclasses.view import View
from pyoptsparse.postprocessing.opt_view_controller import OptViewController

QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)


class MainView(QWidget, View):
    def __init__(self, *args, file_names: List = [], **kwargs):
        """
        OptView - A GUI designed to assist viewing optimization problems
        with pyOptSparse.  The app uses a MVC architecure to modularly
        handle data management, functional control, and visualization.

        Parameters
        ----------
        file_names : List, optional
            List of file names to load on startup, by default []
        """
        super(MainView, self).__init__(*args, **kwargs)
        # Set the controller for the main OptView window
        self._controller = OptViewController()
        self._controller.view = self
        self.file_names = file_names

        self._center()  # centers the application in the middle of the screen
        self.setWindowTitle("OptView")  # sets the GUI title
        self.setWindowIcon(QIcon("assets/OptViewIcon.gif"))  # sets the OptView logo
        self.resize(1200, 800)
        self._initUI()  # Initialize the UI

    def _initUI(self):
        """
        Initializes the user inteface for the MainView widget.
        """
        # --- Set the main layout ---
        self.layout = QVBoxLayout()

        # ==============================================================
        # Menu Bar - First item added to top-level layout
        # ==============================================================
        # --- Add the menu bar object ---
        menu_bar = QMenuBar(self)

        # --- Add file sub-directory with sub-actions ---
        file_menu = menu_bar.addMenu("File")

        self.new_tab_action = QAction("New Tab...", self)
        self.new_tab_action.triggered.connect(lambda: self._controller.addTab(interactive=True))
        file_menu.addAction(self.new_tab_action)

        self.dark_mode_action = QAction("Dark Mode", self, checkable=True)
        self.dark_mode_action.setStatusTip("Toggle Dark/Light Mode")
        self.dark_mode_action.triggered.connect(self._controller.toggleDarkMode)
        file_menu.addAction(self.dark_mode_action)

        self.settings_action = QAction("OptView Settings...", self)
        self.settings_action.setStatusTip("View/Configure OptView settings")
        self.settings_action.triggered.connect(self._controller.configure_settings)
        file_menu.addAction(self.settings_action)

        self.exit_action = QAction("Close", self)
        self.exit_action.triggered.connect(qApp.quit)
        file_menu.addAction(self.exit_action)

        self.layout.addWidget(menu_bar)

        # ==============================================================
        # TabView - The central widget that contains the tabs
        # ==============================================================
        # --- Create the tab control framework ---
        self.tabs = QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self._controller.closeTab)

        # --- Add the first tab to the view ---
        self._controller.addTab(interactive=False, tab_name="Home", file_names=self.file_names)
        self.layout.addWidget(self.tabs)

        # --- Set the layout ---
        self.setLayout(self.layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        """
        Centers the view on the screen.
        """
        qr = self.frameGeometry()
        self.move(qr.center())


def main():
    # Create command line arguments here
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*", type=str, default=[], help="File(s) to load into OptView on startup.")
    args = parser.parse_args()

    # Check to make sure the files exist and are history files
    if args.files:
        for file in args.files:
            if not os.path.exists(file):
                raise FileNotFoundError(f"History file: {file} doesn't exist.")
            elif not file.endswith(".hst"):
                raise NameError(f"File: {file} is not a readable history file.")

    # Setup the app and the UI
    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    MainView(file_names=args.files)

    # This launches the app
    sys.exit(app.exec_())
