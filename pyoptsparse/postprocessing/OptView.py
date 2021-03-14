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
from views.plot_view import PlotView
from controllers.main_controller import MainController

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)


class MainView(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._controller = MainController(self)
        self._initUI()

    def _initUI(self):
        self._center()  # centers the application in the middle of the screen
        self.setWindowTitle("OptView")

        # --- Create top level layout ---
        layout = QtWidgets.QVBoxLayout()

        # --- Create plot view and add to layout ---
        plot = PlotView(self)
        layout.addWidget(plot)

        # --- Create sublayout underneath the plot for buttons, forms, and options ---
        sub_layout = QtWidgets.QHBoxLayout()
        layout.addLayout(sub_layout)

        # --- Create sublayout for right hand side buttons ---
        button_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(button_layout)

        # --- Create layout for x-variables ---
        x_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(x_layout)

        # --- Create layout for y-variables ---
        y_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(y_layout)

        # --- Create layout for options ---
        opt_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(opt_layout)

        # ==============================================================================
        # Button Layout - Left hand column of Sub Layout
        # ==============================================================================
        # --- Add "Add Files" button to the main view ---
        add_file_btn = Button("Add Files", self)
        add_file_btn.setToolTip("Add files to current application")
        button_layout.addWidget(add_file_btn)

        # ==============================================================================
        # X Variable Layout - Left Center column of Sub Layout
        # ==============================================================================
        # --- Add x-vars combobox ---
        x_cbox = ExtendedComboBox(self)
        x_cbox.setToolTip("Type to search for x-variables")
        x_cbox.resize(250, 30)
        x_layout.addWidget(x_cbox)

        # ==============================================================================
        # Y Variable Layout - Right Center column of Sub Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        y_cbox = ExtendedComboBox(self)
        y_cbox.setToolTip("Type to search for y-variables")
        y_cbox.resize(250, 30)
        y_layout.addWidget(y_cbox)

        # ==============================================================================
        # Options Layout - Right hand column of Sub Layout
        # ==============================================================================
        # --- Add options checkboxes ---
        stack_plot_opt = QtWidgets.QCheckBox("Stack plots")
        opt_layout.addWidget(stack_plot_opt)
        opt_layout.setAlignment(stack_plot_opt, QtCore.Qt.AlignCenter)

        # --- Set the main layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        cp = QtWidgets.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = MainView()
    app.exec_()
