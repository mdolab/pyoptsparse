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

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5 import QtWidgets, QtCore
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

# ==============================================================================
# Extension modules
# ==============================================================================
from views.utils.combo_box import ExtendedComboBox
from views.utils.button import Button
from views.mpl_canvas import MplCanvas
from controllers.main_controller import MainController


class MainView(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.left = 10
        self.top = 10
        self.width = 1000
        self.height = 750
        self.title = "OptView"
        self.layout = QtWidgets.QGridLayout()
        self._controller = MainController(self)
        self._initUI()

    def _initUI(self):
        # --- Set the window title and geometry ---
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        # self.setStyleSheet("background-color: white; foreground-color: ;white")

        # --- Add plot to the main view ---
        self._controller.plot_canvas = sc = MplCanvas(self, width=10, height=5, dpi=100)

        # --- Add matplotlib toolbar to the top of the view ---
        toolbar = NavigationToolbar(sc, self)
        toolbar.setMovable(False)
        self.addToolBar(QtCore.Qt.TopToolBarArea, toolbar)

        # --- Add "Add Files" button to the main view ---
        add_file_btn = Button("Add Files", self)
        add_file_btn.setToolTip("Add files to current application")
        self.layout.addWidget(add_file_btn, 4, 0)

        # # --- Add x-vars combobox ---
        # x_cbox = ExtendedComboBox(self)
        # x_cbox.setToolTip("Type to search for x-variables")
        # x_cbox.resize(250, 30)
        # x_cbox.move(col2, row1)

        # # --- Add x-vars combobox ---
        # y_cbox = ExtendedComboBox(self)
        # y_cbox.setToolTip("Type to search for y-variables")
        # y_cbox.resize(250, 30)
        # y_cbox.move(col3, row1)

        # --- Show the view ---
        self.show()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = MainView()
    app.exec_()
