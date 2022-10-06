# External modules
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QFileDialog

# First party modules
from pyoptsparse.postprocessing.sub_windows.configure_plot_window.configure_model import ConfigureModel
from pyoptsparse.postprocessing.sub_windows.variables.x_controller import XController
from pyoptsparse.postprocessing.sub_windows.variables.y_controller import YController
from pyoptsparse.postprocessing.baseclasses.controller import Controller
from pyoptsparse.postprocessing.baseclasses.model import Model
from pyoptsparse.postprocessing.utils.colors import GREEN
from pyoptsparse.postprocessing.utils.widgets import FileTreeWidgetItem, VarTableWidgetItem


class ConfigureController(Controller):
    def __init__(self, parent_model: Model, plot_model: Model):
        """
        Controller for the plot configuration view.  This controller
        facilitates adding x and y variables, file management, and
        plotting.

        Parameters
        ----------
        parent_model : Model
            The tab model that holds the files.
        plot_model : Model
            The plot model linked to the configuration view.
        """
        super(ConfigureController, self).__init__()
        self._parent_model = parent_model
        self._plot_model = plot_model
        self._model = ConfigureModel()
        self._current_file = None

        self._xtable_controller = None
        self._ytable_controller = None

    def setup_var_tables(self):
        """
        Sets up the x and y variable tables.  Creates the controllers
        for each variable table and passes the table views, plot model,
        and tab model.  Finally, sets the controllers for the variable
        table views.
        """
        self._xtable_controller = XController(self._plot_model, self._parent_model)
        self._ytable_controller = YController(self._plot_model, self._parent_model)

        self._view.x_table.setController(self._xtable_controller)
        self._view.y_table.setController(self._ytable_controller)

        self._xtable_controller.set_view(self._view.x_table)
        self._ytable_controller.set_view(self._view.y_table)

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
        Adds all the files to the file tree widget display.  If files
        are loaded, the data from the first file in the tab model is
        displayed by default.
        """
        self._view.file_tree.clear()
        for file in self._parent_model.files:
            file_item = FileTreeWidgetItem(self._view.file_tree)
            file_item.setFile(file)
            file_item.setText(0, file.short_name)
            self._view.file_tree.addTopLevelItem(file_item)

        if len(self._parent_model.files) > 0:
            self._current_file = self._parent_model.files[0]
            self.populate_vars()

    def file_selected(self, item: FileTreeWidgetItem, column: int):
        """
        Grabs the selected file, clears the plot associated with the
        configuration view, and then updates the variables in the
        configuration view based on the selected file's data.

        Parameters
        ----------
        item : FileTreeWidgetItem
            The file tree widget item that is selected.
        column : int
            The tree widget column number.
        """
        self._current_file = item.file
        self._plot_model.clear_vars()
        self._plot_model.clear_axis()
        self._parent_model.canvas.draw()
        self.populate_vars()

    def populate_vars(self):
        """
        Adds all the variables to the  x and y variable tables based
        on the current file selection.  If there is no x-variable in
        the plot model, the default is set to 'iter' and added to the
        plot model.
        """
        self.populate_x_var_combobox()
        self._ytable_controller.clear_vars()
        self._xtable_controller.clear_vars()
        self._ytable_controller.populate_vars(self._current_file)

        if self._plot_model.x_var is None:
            x_var = self._current_file.x_vars["iter"]
            self._view.x_cbox.setCurrentText(x_var.full_name)
            self._plot_model.add_var(x_var, "x")

            var_item = VarTableWidgetItem(x_var.name)
            var_item.var = x_var

            self._xtable_controller.add_row(var_item)

    def populate_x_var_combobox(self):
        """
        Adds the x-variable names to the selection box.
        """
        self._view.x_cbox.clear()

        for var in self._current_file.x_vars.values():
            self._view.x_cbox.addItem(var.full_name)

    def add_x_var(self):
        """
        Adds an x-variable to the plot model and the x-variable table.
        """
        self._xtable_controller.clear_vars()
        var_name = self._view.x_cbox.currentText()

        x_var = self._current_file.x_vars[var_name]
        self._plot_model.add_var(x_var, "x")
        var_item = VarTableWidgetItem(x_var.full_name)
        var_item.var = x_var
        self._xtable_controller.add_row(var_item)

    def y_var_search(self, s: str):
        """
        Searches the y-variable table for the string input.

        Parameters
        ----------
        s : str
            User string search input.
        """
        table = self._view.y_table
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

    def add_selected_vars(self):
        items = self._view.y_table.selectedItems()
        for item in items:
            if item.column() == 0:
                label = self._view.y_table.cellWidget(item.row(), 2).text()
                if str(label) != "":
                    item.getVar().set_label(str(label))
                self._plot_model.add_var(item.getVar(), "y")

                item.setRowColor(GREEN)
                self._ytable_controller.set_row_color(item.row())

        self._plot_model.plot()
        self._parent_model.canvas.draw()

    def rem_selected_vars(self):
        items = self._view.y_table.selectedItems()
        for item in items:
            if item.column() == 0:
                self._plot_model.remove_var(item.getVar(), "y")

                self._view.y_table.cellWidget(item.row(), 3).setChecked(False)
                self._view.y_table.cellWidget(item.row(), 4).setChecked(False)

                item.setRowColor(item.getDefaultRowColor())
                self._ytable_controller.set_row_color(item.row())

        # Update the plot
        self._plot_model.plot()
        self._parent_model.canvas.draw()
