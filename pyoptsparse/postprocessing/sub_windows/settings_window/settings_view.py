# External modules
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import QDialog, QDialogButtonBox, QTableWidget, QTabWidget, QVBoxLayout, QWidget, QHeaderView

# First party modules
from pyoptsparse.postprocessing.baseclasses.controller import Controller
from pyoptsparse.postprocessing.baseclasses.view import View


class SettingsView(QDialog, View):
    def __init__(self, parent: QWidget = None, controller: Controller = None):
        super(SettingsView, self).__init__(parent)
        self._controller = controller
        self._controller.view = self
        self.resize(800, 800)

        self._initUI()

    def _initUI(self):
        self.setWindowTitle("OptView Settings")
        layout = QVBoxLayout()

        # --- Create the tab control framework ---
        self.tabs = QTabWidget()
        self.tabs.setTabsClosable(False)

        # --- Add the first tab to the view ---
        self.tabs.addTab(KeyboardShortuctsView(parent=self, controller=self._controller), "Keyboard Shortcuts")
        self.tabs.addTab(MplRcParametersView(parent=self, controller=self._controller), "Matplotlib RC Parameters")
        layout.addWidget(self.tabs, 2)

        btnBox = QDialogButtonBox()
        btnBox.setStandardButtons(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btnBox.accepted.connect(self.accept)
        btnBox.rejected.connect(self.reject)
        layout.addWidget(btnBox)

        # --- Set the layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()


class KeyboardShortuctsView(QWidget, View):
    def __init__(self, parent: QDialog = None, controller: Controller = None):
        super(KeyboardShortuctsView, self).__init__(parent)
        self._controller = controller
        self._controller.keyshort_view = self

        self._initUI()

    def _initUI(self):
        layout = QVBoxLayout()
        self.shortcut_table = QTableWidget(self)
        font = QFont()
        font.setPointSize(16)
        self.shortcut_table.setFont(font)
        self.shortcut_table.setShowGrid(False)
        self.shortcut_table.verticalHeader().setVisible(False)
        self.shortcut_table.setColumnCount(2)
        self.shortcut_table.setHorizontalHeaderLabels(["Sequence", "Description"])
        self.shortcut_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.shortcut_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self.shortcut_table)

        self.setLayout(layout)


class MplRcParametersView(QWidget, View):
    def __init__(self, parent: QDialog = None, controller: Controller = None):
        super(MplRcParametersView, self).__init__(parent)

        self._controller = controller
        self._controller.rc_param_view = self

        self._initUI()

    def _initUI(self):
        layout = QVBoxLayout()
        self.rc_table = QTableWidget(self)
        font = QFont()
        font.setPointSize(16)
        self.rc_table.setFont(font)
        self.rc_table.setShowGrid(False)
        self.rc_table.verticalHeader().setVisible(False)
        self.rc_table.setColumnCount(2)
        self.rc_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self.rc_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.rc_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self.rc_table)

        self.setLayout(layout)
