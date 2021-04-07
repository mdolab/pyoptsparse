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
from PyQt5 import QtWidgets, QtCore, QtGui
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
        self.setWindowTitle("OptView")  # sets the GUI title
        self.setWindowIcon(QtGui.QIcon("assets/OptViewIcon.gif"))  # sets the OptView logo

        # --- Create top level layout ---
        layout = QtWidgets.QVBoxLayout()

        # ==============================================================================
        # Menu Bar - First item added to top-level layout
        # ==============================================================================
        # --- Add the menu bar object ---
        menu_bar = QtWidgets.QMenuBar(self)

        # --- Add file sub-directory with sub-actions ---
        file_menu = menu_bar.addMenu("File")

        new_window_action = QtWidgets.QAction("New Window", self)
        new_window_action.triggered.connect(self._controller.newWindow)
        file_menu.addAction(new_window_action)

        load_action = QtWidgets.QAction("Load File...", self)
        load_action.triggered.connect(self._controller.openFile)
        file_menu.addAction(load_action)

        save_tec_action = QtWidgets.QAction("Save As Tec File", self)
        save_tec_action.triggered.connect(self._controller.saveTecFile)
        file_menu.addAction(save_tec_action)

        exit_action = QtWidgets.QAction("Exit", self)
        exit_action.triggered.connect(QtWidgets.qApp.quit)
        file_menu.addAction(exit_action)

        # --- Add format sub-directory with sub-actions ---
        format_menu = menu_bar.addMenu("Format")

        font_action = QtWidgets.QAction("Font size", self)
        font_action.triggered.connect(self._controller.changePlotFontSize)
        format_menu.addAction(font_action)

        refresh_plot_action = QtWidgets.QAction("Refresh Plot", self)
        refresh_plot_action.triggered.connect(self._controller.refreshPlot)
        format_menu.addAction(refresh_plot_action)

        clear_plot_action = QtWidgets.QAction("Clear Plot", self)
        clear_plot_action.triggered.connect(self._controller.clearPlot)
        format_menu.addAction(clear_plot_action)

        layout.addWidget(menu_bar)

        # --- Create plot view and add to layout ---
        plot = PlotView(self)
        layout.addWidget(plot)

        # --- Create sublayout underneath the plot for buttons, forms, and options ---
        sub_layout = QtWidgets.QHBoxLayout()
        layout.addLayout(sub_layout)

        # --- Create sublayout for x-variables ---
        x_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(x_layout)

        # --- Create sublayout for y-variables ---
        y_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(y_layout)

        # --- Create sublayout for LHS options ---
        opt1_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(opt1_layout)

        # --- Create sublayout for RHS options ---
        opt2_layout = QtWidgets.QVBoxLayout()
        sub_layout.addLayout(opt2_layout)

        # ==============================================================================
        # X Variable Layout - Left Center column of Sub Layout
        # ==============================================================================
        # --- Add x-vars combobox ---
        self.x_cbox = ExtendedComboBox(self)
        self.x_cbox.setToolTip("Type to search for x-variables")
        self.x_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.x_cbox.resize(250, 30)
        self.x_cbox.addItems(["test1", "test2", "test3"])
        self.x_cbox.currentIndexChanged.connect(self._controller.addVarX)
        x_layout.addWidget(x_cbox)

        # --- Add x-vars variable list ---
        self.x_label = QtWidgets.QLabel(self)
        self.x_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self.x_label.resize(250, 100)
        x_layout.addWidget(x_label)

        # --- Add undo x-vars button ---
        self.x_undo_btn = Button("Undo x-var", self)
        self.x_undo_btn.setToolTip("Undo add x-variable")
        self.x_undo_btn.clicked.connect(self._controller.undoVarX)
        x_layout.addWidget(x_undo_btn)

        # --- Add clear x-vars button ---
        self.x_clear_btn = Button("Clear x-var", self)
        self.x_clear_btn.setToolTip("Clear all x-variables")
        self.x_clear_btn.clicked.connect(self._controller.clearAllX)
        x_layout.addWidget(x_clear_btn)

        # ==============================================================================
        # Y Variable Layout - Right Center column of Sub Layout
        # ==============================================================================
        # --- Add y-vars combobox ---
        self.y_cbox = ExtendedComboBox(self)
        self.y_cbox.setInsertPolicy(QtWidgets.QComboBox.NoInsert)
        self.y_cbox.setToolTip("Type to search for y-variables")
        self.y_cbox.resize(250, 30)
        self.y_cbox.addItems(["test1", "test2", "test3"])
        self.y_cbox.currentIndexChanged.connect(self._controller.addVarY)
        y_layout.addWidget(y_cbox)

        # --- Add y-vars variable list ---
        self.y_label = QtWidgets.QLabel(self)
        self.y_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self.y_label.resize(250, 100)
        y_layout.addWidget(y_label)

        # --- Add undo y-vars button ---
        self.y_undo_btn = Button("Undo y-var", self)
        self.y_undo_btn.setToolTip("Undo add y-variable")
        self.y_undo_btn.clicked.connect(self._controller.undoVarY)
        y_layout.addWidget(y_undo_btn)

        # --- Add clear y-vars button ---
        self.y_clear_btn = Button("Clear y-var", self)
        self.y_clear_btn.setToolTip("Clear all y-variables")
        self.y_clear_btn.clicked.connect(self._controller.clearAllY)
        y_layout.addWidget(y_clear_btn)

        # ==============================================================================
        # Options Layout 1 - First sub-layout column for options
        # ==============================================================================
        # --- Stacked Plots ---
        self.stack_plot_opt = QtWidgets.QCheckBox("Stack plots")
        self.stack_plot_opt.clicked.connect(self._controller.stackPlots)
        opt1_layout.addWidget(self.stack_plot_opt)
        opt1_layout.setAlignment(self.stack_plot_opt, QtCore.Qt.AlignLeft)

        # --- Shared y-axis ---
        self.share_y_opt = QtWidgets.QCheckBox("Shared y-axis")
        self.share_y_opt.clicked.connect(self._controller.shareAxisY)
        opt1_layout.addWidget(self.share_y_opt)
        opt1_layout.setAlignment(self.share_y_opt, QtCore.Qt.AlignLeft)

        # --- y-axis as absolute delta values ---
        self.abs_delta_opt = QtWidgets.QCheckBox("y-axis as absolute delta values")
        self.abs_delta_opt.clicked.connect(self._controller.absDeltaY)
        opt1_layout.addWidget(self.abs_delta_opt)
        opt1_layout.setAlignment(self.abs_delta_opt, QtCore.Qt.AlignLeft)

        # --- x-axis as minor iterations ---
        self.minor_itr_opt = QtWidgets.QCheckBox("x-axis as minor iterations")
        self.minor_itr_opt.clicked.connect(self._controller.minorIterX)
        opt1_layout.addWidget(self.minor_itr_opt)
        opt1_layout.setAlignment(self.minor_itr_opt, QtCore.Qt.AlignLeft)

        # --- x-axis as major iterations ---
        self.major_itr_opt = QtWidgets.QCheckBox("x-axis as major iterations")
        self.major_itr_opt.clicked.connect(self._controller.majorIterX)
        opt1_layout.addWidget(self.major_itr_opt)
        opt1_layout.setAlignment(self.major_itr_opt, QtCore.Qt.AlignLeft)

        # --- Apply scaling factor ---
        self.scale_factor_opt = QtWidgets.QCheckBox("Apply Scaling Factor")
        self.scale_factor_opt.clicked.connect(self._controller.scaleFactor)
        opt1_layout.addWidget(self.scale_factor_opt)
        opt1_layout.setAlignment(self.scale_factor_opt, QtCore.Qt.AlignLeft)

        # ==============================================================================
        # Options Layout 2 - Second sub-layout column for options
        # ==============================================================================
        # --- Min/Max arrays ---
        self.min_max_opt = QtWidgets.QCheckBox("Min/Max for arrays")
        self.min_max_opt.clicked.connect(self._controller.minMaxPlot)
        opt2_layout.addWidget(self.min_max_opt)
        opt2_layout.setAlignment(self.min_max_opt, QtCore.Qt.AlignLeft)

        # --- Auto refresh data ---
        self.auto_refresh_opt = QtWidgets.QCheckBox("Automatically refresh history")
        self.auto_refresh_opt.clicked.connect(self._controller.autoRefresh)
        opt2_layout.addWidget(self.auto_refresh_opt)
        opt2_layout.setAlignment(self.auto_refresh_opt, QtCore.Qt.AlignLeft)

        # --- Set the main layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()

    def _center(self):
        qr = self.frameGeometry()
        self.move(qr.center())


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = MainView()
    app.exec_()
