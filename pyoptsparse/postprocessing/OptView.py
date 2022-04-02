# --- Python 3.8 ---
"""
OptView - A GUI designed to assist viewing optimization problems with
pyOptSparse.  The app uses a MVC architecure to modularly handle
data management, functional control, and visualization.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import sys
import os
import argparse

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5 import QtWidgets, QtGui, QtCore
from click import FileError

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse.postprocessing.sub_MVCs.tab_window.tab_view import TabView
from pyoptsparse.postprocessing.opt_view_controller import OptViewController

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)


class MainView(QtWidgets.QWidget):
    def __init__(self, *args, files=[], **kwargs):
        super().__init__(*args, **kwargs)
        self._center()  # centers the application in the middle of the screen
        self.setWindowTitle("OptView")  # sets the GUI title
        self.setWindowIcon(QtGui.QIcon("assets/OptViewIcon.gif"))  # sets the OptView logo
        self._controller = OptViewController(self)
        self.files = files
        self._initUI()

    def _initUI(self):
        # --- Set the main layout ---
        self.layout = QtWidgets.QVBoxLayout()

        # ==============================================================================
        # Menu Bar - First item added to top-level layout
        # ==============================================================================
        # --- Add the menu bar object ---
        menu_bar = QtWidgets.QMenuBar(self)

        # --- Add file sub-directory with sub-actions ---
        file_menu = menu_bar.addMenu("File")

        new_window_action = QtWidgets.QAction("New Tab", self)
        new_window_action.triggered.connect(self._controller.addTab)
        file_menu.addAction(new_window_action)

        exit_action = QtWidgets.QAction("Exit", self)
        exit_action.triggered.connect(QtWidgets.qApp.quit)
        file_menu.addAction(exit_action)

        self.layout.addWidget(menu_bar)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self._controller.closeTab)
        tab1 = TabView(files=self.files, parent=self)

        self.tabs.addTab(tab1, "Home")

        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", nargs="+", type=str, help="File(s) to load into OptView on startup.")
    args = parser.parse_args()
    for file in args.files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"History file: {file}")
        elif not file.endswith(".hst"):
            raise NameError(f"File: {file} is not a readable history file.")
    app = QtWidgets.QApplication(sys.argv)
    w = MainView(files=args.files)
    app.exec_()
