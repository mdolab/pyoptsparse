# External modules
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLineEdit,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

# First party modules
from pyoptsparse.postprocessing.baseclasses import Controller, Model, View
from pyoptsparse.postprocessing.general_widgets import Button, FileTreeWidgetItem


class OptTreeWidgetItem(QTreeWidgetItem):
    def __init__(self, *args, **kwargs):
        """
        Custom tree widget item for metadata options.
        """
        super().__init__(*args, **kwargs)


class OptTableWidgetItem(QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        """
        Custom table widget item for metatdata options.
        """
        super().__init__(*args, **kwargs)


class MetadataModel(Model):
    def __init__(self):
        """
        The metadata window model.
        """
        super(MetadataModel, self).__init__()


class MetadataController(Controller):
    def __init__(self, model: Model, parent_model: Model):
        """
        The controller for the metadata options view.

        Parameters
        ----------
        model : Model
            The metadata options model.
        files : List
            The list of files.
        """
        self._model = model
        self._parent_model = parent_model
        self._current_file = None

    def add_file(self):
        """
        Opens a file dialog to get user input for adding a new history
        file.
        """
        # --- Open file dialog and get selected user files ---
        file_names, _ = QFileDialog.getOpenFileNames(
            self._view, "Open History File", "", "History File (*.hst)", options=QFileDialog.Option.DontUseNativeDialog
        )

        # --- Load files into the model ---
        self._parent_model.load_files(file_names)

        # --- Populate the files ---
        self.populate_files()

    def populate_files(self):
        """
        Populates the file tree with loaded files from the tab view model.
        """
        for file in self._parent_model.files:
            file_item = FileTreeWidgetItem(self._view.file_tree)
            file_item.setFile(file)
            file_item.setText(0, file.short_name)
            self._view.file_tree.addTopLevelItem(file_item)

        if len(self._parent_model.files) > 0:
            self._current_file = self._parent_model.files[0]
            self.populate_opts()

    def file_selected(self, item: FileTreeWidgetItem, column: int):
        """
        Populates the options display widgets when a new file is
        selected.

        Parameters
        ----------
        item : FileTreeWidgetItem
            The selected file tree widget item
        column : int
            The file tree column index
        """
        self._current_file = item.file
        self.populate_opts()

    def search(self, s: str):
        """
        Searches the options for the string input.

        Parameters
        ----------
        s : str
            String search input.
        """
        table = self._view.opt_prob_table
        row_count = table.rowCount()
        sel_items = table.findItems(s, Qt.MatchFlag.MatchContains)

        rows_to_show = set()
        for item in sel_items:
            rows_to_show.add(item.row())

        for row in rows_to_show:
            table.setRowHidden(row, False)

        for row in range(row_count):
            if row not in rows_to_show:
                table.setRowHidden(row, True)

    def populate_opts(self):
        """
        Populates the options widgets with the history file's
        metadata.
        """
        self.clear_opts()

        metadata = self._current_file.metadata
        for key, val in metadata.items():
            if key != "optOptions":
                item = OptTreeWidgetItem(self._view.opt_tree)
                item.setText(0, key)
                item.setText(1, f"{val}")
                self._view.opt_tree.addTopLevelItem(item)

        self._view.opt_prob_table.setRowCount(len(metadata["optOptions"].keys()))
        self._view.opt_prob_table.setColumnCount(2)
        for i, (key, val) in enumerate(metadata["optOptions"].items()):
            option = OptTableWidgetItem(key)
            value = OptTableWidgetItem(f"{val}")

            self._view.opt_prob_table.setItem(i, 0, option)
            self._view.opt_prob_table.setItem(i, 1, value)
        self._view.opt_prob_table.resizeColumnsToContents()
        self._view.opt_prob_table.setHorizontalHeaderLabels(["Option", "Value"])
        self._view.opt_prob_table.verticalHeader().setVisible(False)

    def clear_opts(self):
        # Clear everything
        self._view.opt_tree.clear()
        self._view.opt_prob_table.clear()


class MetadataView(QDialog, View):
    def __init__(self, parent: QWidget, controller: MetadataController, name: str):
        """
        The view for the displaying metadata options.

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QWidget
            The parent tab view.
        controller : Controller
            The metadata view controller linked to this view.
        name : str
            The name of the window.
        """
        super(MetadataView, self).__init__(parent)
        self._center()
        self.setWindowTitle(name)
        self._controller = controller
        self._controller.view = self
        self.resize(1000, 800)
        self._initView()

    def _initView(self):
        """
        Initializes the view.
        """
        # --- Create layout ---
        layout = QHBoxLayout()

        left_layout = QVBoxLayout()
        layout.addLayout(left_layout)

        right_layout = QVBoxLayout()
        layout.addLayout(right_layout, 3)

        # ==============================================================
        # File List - Left Layout
        # ==============================================================
        # --- Add file button ---
        self.add_file_btn = Button("Add file(s)", self)
        self.add_file_btn.clicked.connect(self._controller.add_file)
        left_layout.addWidget(self.add_file_btn)

        # --- File list ---
        self.file_tree = QTreeWidget(self)
        self.file_tree.setColumnCount(1)
        self.file_tree.setHeaderLabels(["File Name"])
        self.file_tree.itemClicked.connect(self._controller.file_selected)
        left_layout.addWidget(self.file_tree)

        # ==============================================================
        # Options Table - Right Layout
        # ==============================================================
        self.opt_tree = QTreeWidget(self)
        self.opt_tree.setColumnCount(2)
        self.opt_tree.setHeaderLabels(["Name", "Value"])
        right_layout.addWidget(self.opt_tree)

        self.query = QLineEdit()
        self.query.setPlaceholderText("Search...")
        self.query.textChanged.connect(self._controller.search)
        right_layout.addWidget(self.query)

        self.opt_prob_table = QTableWidget(self)
        self.opt_prob_table.setShowGrid(False)
        right_layout.addWidget(self.opt_prob_table)

        self._controller.populate_files()

        self.setLayout(layout)

        self.show()

    def _center(self):
        """
        Centers the view on the screen.
        """
        qr = self.frameGeometry()
        self.move(qr.center())
