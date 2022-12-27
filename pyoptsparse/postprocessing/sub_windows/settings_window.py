# Standard Python modules
import os
from typing import Dict

# External modules
from PyQt6.QtGui import QFont, QKeySequence
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHeaderView,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

# First party modules
from pyoptsparse.postprocessing.baseclasses import Controller, Model, View
from pyoptsparse.postprocessing.utils import ASSET_PATH


class SettingsModel(Model):
    def __init__(self, *args, **kwargs):
        super(SettingsModel, self).__init__(*args, **kwargs)
        self._shortcuts = {}

    @property
    def shortcuts(self) -> Dict:
        return self._shortcuts

    def add_shortcut(self, shortcut: Dict):
        self._shortcuts = {**self._shortcuts, **shortcut}


class SettingsController(Controller):
    def __init__(self):
        """
        The controller for the tab view.

        Parameters
        ----------
        root : PyQt6.QtWidgets.QWidget
            The OptView main view
        file_names : List, optional
            Names of files to be pre-loaded in to the model,
            by default []
        """
        super(SettingsController, self).__init__()
        self._appearance_view = None
        self._keyshort_view = None
        self._rc_param_view = None
        self._model = SettingsModel()

    @property
    def appearance_view(self):
        return self._appearance_view

    @property
    def keyshort_view(self):
        return self._keyshort_view

    @property
    def rc_param_view(self):
        return self._rc_param_view

    @appearance_view.setter
    def appearance_view(self, view):
        self._appearance_view = view

    @keyshort_view.setter
    def keyshort_view(self, view):
        self._keyshort_view = view

    @rc_param_view.setter
    def rc_param_view(self, view):
        self._rc_param_view = view

    def populate_rc_params(self):
        with open(os.path.join(ASSET_PATH, "nicePlotsStyle"), "r") as file:
            for line in file:
                if line != "\n":
                    current_row_count = self._rc_param_view.rc_table.rowCount()
                    self._rc_param_view.rc_table.insertRow(current_row_count)
                    param_item = QTableWidgetItem(line.split(":")[0])
                    val_item = QTableWidgetItem(line.split(":")[1].strip("\n"))
                    self._rc_param_view.rc_table.setItem(current_row_count, 0, param_item)
                    self._rc_param_view.rc_table.setItem(current_row_count, 1, val_item)

    def populate_shortcuts(self):
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.Open).toString(), "Opens the add file menu.")
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.New).toString(), "Add a subplot to the figure.")
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.Find).toString(), "Opens the figure options menu.")
        self._add_shortcut("Ctrl+T", "Applies the tight layout format to the figure.")
        self._add_shortcut(QKeySequence(QKeySequence.StandardKey.Save).toString(), "Opens the save figure menu.")
        self._add_shortcut(
            QKeySequence(QKeySequence.StandardKey.SelectStartOfDocument).toString(),
            "Resets the figure to the default home view.",
        )
        self._add_shortcut("Ctrl+Up", "Moves a subplot up.")
        self._add_shortcut("Ctrl+Down", "Moves a subplot down.")
        self._add_shortcut(QKeySequence("Alt+h").toString(), "Opens the help/settings menu")

    def _add_shortcut(self, key, val):
        self._model.add_shortcut({key: val})
        current_row_count = self._keyshort_view.shortcut_table.rowCount()
        self._keyshort_view.shortcut_table.insertRow(current_row_count)
        key_item = QTableWidgetItem(key)
        val_item = QTableWidgetItem(val)
        self._keyshort_view.shortcut_table.setItem(current_row_count, 0, key_item)
        self._keyshort_view.shortcut_table.setItem(current_row_count, 1, val_item)


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


class KeyboardShortuctsView(View):
    def __init__(self, parent: QDialog = None, controller: SettingsController = None):
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


class MplRcParametersView(View):
    def __init__(self, parent: QDialog = None, controller: SettingsController = None):
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
