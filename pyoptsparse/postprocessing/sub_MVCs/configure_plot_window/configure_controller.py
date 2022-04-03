# Standard Python modules

# External modules
from PyQt5.QtWidgets import QFileDialog, QDialog
from PyQt5.QtCore import Qt

# Local modules
from pyoptsparse.postprocessing.utils.widgets import FileTreeWidgetItem, VarTableWidgetItem, FileTableWidgetItem
from pyoptsparse.postprocessing.utils.data_structures import Variable
from pyoptsparse.postprocessing.utils.base_classes import Controller, Model
from pyoptsparse.postprocessing.sub_MVCs.configure_plot_window.configure_model import ConfigureModel
from pyoptsparse.postprocessing.sub_MVCs.variables.y_controller import YController
from pyoptsparse.postprocessing.sub_MVCs.variables.x_controller import XController


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
        self._view = None
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

    def set_view(self, view: QDialog):
        """
        Sets the view for the controller.

        Parameters
        ----------
        view : PyQt5.QtWidgets.QDialog
            The configure plot view to link to this controller.
        """
        self._view = view

    def add_file(self):
        """
        Opens a file dialog to get user input for adding a new history
        file.
        """
        # --- Set file dialog options ---
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        # --- Open file dialog and get selected user files ---
        file_names, _ = QFileDialog.getOpenFileNames(self._view, "Open History File", "", "", options=options)

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
            file_item.setText(0, file.name_short)
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
            new_var = Variable("iter")
            new_var.file = self._current_file
            self._plot_model.add_var(new_var, "x")

        self._xtable_controller.populate_vars()

    def populate_x_var_combobox(self):
        """
        Adds the x-variable names to the selection box.
        """
        self._view.x_cbox.clear()

        for var_name in self._current_file.get_all_x_var_names():
            self._view.x_cbox.addItem(var_name)

    def add_x_var(self):
        """
        Adds an x-variable to the plot model and the x-variable table.
        """
        self._xtable_controller.clear_vars()
        var_name = self._view.x_cbox.currentText()

        # Create a new variable and add to the plot model
        var = Variable(var_name)
        var.file = self._current_file
        self._plot_model.add_var(var, "x")

        # Create a new variable widget item for the tree view
        file_item = FileTableWidgetItem(var.file.name_short)
        var_item = VarTableWidgetItem(var.name)
        var_item.var = var

        self._xtable_controller.add_row(file_item, var_item)

    def y_var_search(self, s: str):
        """
        Searches the y-variable table for the string input.

        Parameters
        ----------
        s : str
            User string search input.
        """
        items = self._view.y_table.findItems(s, Qt.MatchContains)
        if items:
            item = items[0]
            self._view.y_table.setCurrentItem(item)
