# External modules
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QDialogButtonBox, QListWidget, QListWidgetItem, QTabWidget, QVBoxLayout, QWidget

# First party modules
from pyoptsparse.postprocessing.baseclasses.controller import Controller
from pyoptsparse.postprocessing.baseclasses.view import View
from pyoptsparse.postprocessing.utils.widgets import DarkModeListWidget


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
        self.tabs.addTab(AppearanceSettingsView(parent=self, controller=self._controller), "Appearance")
        self.tabs.addTab(KeyboardShortuctsView(parent=self, controller=self._controller), "Keyboard Shortcuts")
        self.tabs.addTab(MplRcParametersView(parent=self, controller=self._controller), "Matplotlib RC Parameters")
        layout.addWidget(self.tabs, 3)

        btnBox = QDialogButtonBox()
        btnBox.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btnBox.accepted.connect(self.accept)
        btnBox.rejected.connect(self.reject)
        layout.addWidget(btnBox)

        # --- Set the layout ---
        self.setLayout(layout)

        # --- Show the view ---
        self.show()


class AppearanceSettingsView(QWidget, View):
    def __init__(self, parent: QDialog = None, controller: Controller = None):
        super(AppearanceSettingsView, self).__init__(parent)
        self._controller = controller
        self._controller.appearance_view = self

        self._initUI()

    def _initUI(self):
        layout = QVBoxLayout()

        self.widget_list = QListWidget(self)
        layout.addWidget(self.widget_list)

        dark_mode_item = QListWidgetItem(self.widget_list)
        dark_mode_item.setFlags(Qt.NoItemFlags)
        dark_mode_widget = DarkModeListWidget(self.widget_list, self._controller)
        dark_mode_item.setSizeHint(dark_mode_widget.sizeHint())
        self.widget_list.addItem(dark_mode_item)
        self.widget_list.setItemWidget(dark_mode_item, dark_mode_widget)

        self.setLayout(layout)


class KeyboardShortuctsView(QWidget, View):
    def __init__(self, parent: QDialog = None, controller: Controller = None):
        super(KeyboardShortuctsView, self).__init__(parent)
        self._controller = controller
        self._controller.keyshort_view = self

        self._initUI()

    def _initUI(self):
        layout = QVBoxLayout()
        self.shortcut_list = QListWidget(self)
        layout.addWidget(self.shortcut_list)

        self.setLayout(layout)


class MplRcParametersView(QWidget, View):
    def __init__(self, parent: QDialog = None, controller: Controller = None):
        super(MplRcParametersView, self).__init__(parent)

        self._controller = controller
        self._controller.rc_param_view = self

        self._initUI()

    def _initUI(self):
        layout = QVBoxLayout()
        self.param_list = QListWidget(self)
        layout.addWidget(self.param_list)

        self.setLayout(layout)
