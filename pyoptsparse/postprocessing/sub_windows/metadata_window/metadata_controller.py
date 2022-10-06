# External modules
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QFileDialog

# First party modules
from pyoptsparse.postprocessing.baseclasses.controller import Controller
from pyoptsparse.postprocessing.baseclasses.model import Model
from pyoptsparse.postprocessing.utils.widgets import FileTreeWidgetItem, OptTableWidgetItem, OptTreeWidgetItem


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
        # --- Set file dialog options ---
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        # --- Open file dialog and get selected user files ---
        file_names, _ = QFileDialog.getOpenFileNames(
            self._view, "Open History File", "", "History File (*.hst)", options=options
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
        sel_items = table.findItems(s, Qt.MatchContains)

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
