# External modules
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)

# First party modules
from pyoptsparse.postprocessing.baseclasses import Controller, Model, View
from pyoptsparse.postprocessing.colors import BLUE, GREEN
from pyoptsparse.postprocessing.data_structures import File, Variable
from pyoptsparse.postprocessing.general_widgets import Button, ExtendedComboBox, FileTreeWidgetItem, Switch


class VarTableWidgetItem(QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        """
        Custom variable table widget item.
        """
        super(VarTableWidgetItem, self).__init__(*args, **kwargs)
        self._var = None
        self._row_color = None
        self._default_row_color = self.background()

    def __lt__(self, other):
        if self._var.name == other._var.name:
            return self._var.idx < other._var.idx
        else:
            return self._var.name < other._var.name

    @property
    def var(self) -> Variable:
        return self._var

    @var.setter
    def var(self, var: Variable) -> None:
        self._var = var

    @property
    def row_color(self) -> QColor:
        return self._row_color

    @row_color.setter
    def row_color(self, row_color: QColor) -> None:
        self._row_color = row_color

    @property
    def default_row_color(self) -> QColor:
        return self._default_row_color


class YTableWidget(QTableWidget, View):
    def __init__(self, parent: QDialog = None):
        """
        Defines the y-variable table view.

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QDialog, optional
            The parent view, by default None
        """
        super(YTableWidget, self).__init__(parent)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.setColumnCount(6)
        self.setHorizontalHeaderLabels(["Name", "Index", "Scaled", "Bounds", "Add", "Remove"])
        self._controller = None

    def setController(self, controller):
        self._controller = controller


class YController(Controller):
    def __init__(self, plot_model: Model, parent_model: Model):
        """
        The controller for the y-variable table view.

        Parameters
        ----------
        plot_model : Model
            The plot model to add the y-variables.
        parent_model : Model
            The parent model of the plot.
        """
        super(YController, self).__init__()
        self._view = None
        self._parent_model = parent_model
        self._plot_model = plot_model

    def populate_vars(self, current_file: File):
        """
        Populates the y-variable table view with all the y-variables
        in the current file.

        Parameters
        ----------
        current_file : File
            The current file selected by the user.
        """
        for var in current_file.y_vars.values():
            self.add_row(var)

        self._view.sortItems(0, Qt.SortOrder.AscendingOrder)

        # Find all variables that are in the plot, highlight them in
        # green, and link the table variable to the plot variable
        for row in range(self._view.rowCount()):
            var_item = self._view.item(row, 0)
            self.set_alternating_color(var_item, row)
            if any(var_item.var == y_var for y_var in self._plot_model.y_vars):
                var_item.row_color = GREEN
                self.set_row_color(row)

    def set_alternating_color(self, item, row):
        if row > 0:
            item_prev = self._view.item(row - 1, 0)
            color_prev = item_prev.row_color

            if color_prev == BLUE:
                if item.var.name == item_prev.var.name:
                    item.row_color = BLUE
                    self.set_row_color(row)

            else:
                if item.var.name != item_prev.var.name:
                    item.row_color = BLUE
                    self.set_row_color(row)

    def add_row(self, var: Variable):
        """
        Adds a row to the y-variable table view.

        Parameters
        ----------
        file_item : FileTableWidgetItem
            The file widget item being added.
        var_item : VarTableWidgetItem
            The variable widget item being added.
        """
        row = self._view.rowCount()
        self._view.setRowCount(row + 1)

        var_item = VarTableWidgetItem(var.name)
        var_item.var = var
        self._view.setItem(row, 0, var_item)

        idx_item = QTableWidgetItem(f"{var_item.var.idx}")
        idx_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        self._view.setItem(row, 1, idx_item)

        scaled_opt_item = QTableWidgetItem()
        scaled_opt_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        scaled_opt_chbx = QCheckBox(self._view)
        scaled_opt_chbx.setChecked(False)
        scaled_opt_chbx.stateChanged.connect(self.scale_opt_set)
        self._view.setItem(row, 2, scaled_opt_item)
        self._view.setCellWidget(row, 2, scaled_opt_chbx)

        bounds_opt_item = QTableWidgetItem()
        bounds_opt_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        bounds_opt_chbx = QCheckBox(self._view)
        bounds_opt_chbx.setChecked(False)
        bounds_opt_chbx.stateChanged.connect(self.bounds_opt_set)
        self._view.setItem(row, 3, bounds_opt_item)
        self._view.setCellWidget(row, 3, bounds_opt_chbx)

        add_btn = Button("Add", self._view)
        add_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        add_item = QTableWidgetItem()
        add_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        add_btn.clicked.connect(self.add_var_to_plot)
        self._view.setItem(row, 4, add_item)
        self._view.setCellWidget(row, 4, add_btn)

        rem_btn = Button("Remove", self._view)
        rem_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        rem_item = QTableWidgetItem()
        rem_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
        rem_btn.clicked.connect(self.remove_var_from_plot)
        self._view.setItem(row, 5, rem_item)
        self._view.setCellWidget(row, 5, rem_btn)

        self._view.setHorizontalHeaderLabels(["Name", "Index", "Scaled", "Bounds", "Add", "Remove"])
        self._view.resizeColumnsToContents()

    def clear_vars(self):
        """
        Clears the y-variable table view and resets the row count to
        zero.
        """
        self._view.clear()
        self._view.setRowCount(0)

    def scale_opt_set(self):
        """
        Sets the scale option for the selected variable and re-plots.
        """
        checkbox = self._view.sender()
        index = self._view.indexAt(checkbox.pos())
        selected_item = self._view.item(index.row(), 0)
        scaled_opt = checkbox.checkState()

        selected_item.var.options["scaled"] = scaled_opt

        self._plot_model.plot()
        self._parent_model.draw_canvas()

    def bounds_opt_set(self):
        """
        Sets the bounds options for the selected variable and re-plots
        """
        checkbox = self._view.sender()
        index = self._view.indexAt(checkbox.pos())
        selected_item = self._view.item(index.row(), 0)
        bounds_opt = checkbox.checkState()

        selected_item.var.options["bounds"] = bounds_opt

        self._plot_model.plot()
        self._parent_model.draw_canvas()

    def add_var_to_plot(self):
        """
        Adds a y-variable to the plot model and re-plots.
        """
        button = self._view.sender()
        index = self._view.indexAt(button.pos())
        selected_item = self._view.item(index.row(), 0)
        var = selected_item.var

        if not self._plot_model.x_var.options["major"]:
            var.options["major"] = False

        self._plot_model.add_var(var, "y")

        self._plot_model.plot()
        self._parent_model.draw_canvas()

        selected_item.row_color = GREEN
        self.set_row_color(index.row())

    def remove_var_from_plot(self):
        """
        Removes a y-variable from the plot model and re-plots
        """
        button = self._view.sender()
        index = self._view.indexAt(button.pos())
        selected_item = self._view.item(index.row(), 0)
        self._plot_model.remove_var(selected_item.var, "y")

        self._view.cellWidget(index.row(), 3).setChecked(False)
        self._view.cellWidget(index.row(), 4).setChecked(False)

        # Update the plot
        self._plot_model.plot()
        self._parent_model.draw_canvas()

        selected_item.row_color = selected_item.default_row_color
        self.set_row_color(index.row())

    def set_row_color(self, row: int):
        """
        Sets a given row to a specific color.

        Parameters
        ----------
        row : int
            The row being colored.
        color : QtGui.QColor
            The color for the row.
        """
        color = self._view.item(row, 0).row_color
        for j in range(self._view.columnCount()):
            self._view.item(row, j).setBackground(color)


class XTableWidget(QTableWidget, View):
    def __init__(self, parent: QDialog = None):
        """
        Defines a x-variable table view.

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QDialog, optional
            The parent view, by default None
        """
        super(XTableWidget, self).__init__(parent)
        self.setColumnCount(2)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self._controller = None

    def setController(self, controller: Controller):
        """
        Sets the controller for this view.

        Parameters
        ----------
        controller : Controller
            The controller for this view.
        """
        self._controller = controller


class XController(Controller):
    def __init__(self, plot_model: Model, parent_model: Model):
        """
        The controller for the x-variables view.

        Parameters
        ----------
        plot_model : Model
            The plot model where variables will be added.
        parent_model : Model
            The tab view model that contains the plot.
        """
        super(XController, self).__init__()
        self._view = None
        self._parent_model = parent_model
        self._plot_model = plot_model

    def add_row(self, var_item: VarTableWidgetItem):
        """
        Adds a row to the table view formatted specifically for
        x-variables.

        Parameters
        ----------
        file_item : FileTableWidgetItem
            The file table widget item being added.
        var_item : VarTableWidgetItem
            The variable table widget item being added.
        """
        row = self._view.rowCount()
        self._view.setRowCount(row + 1)
        self._view.setItem(row, 0, var_item)

        iter_switch = Switch(self._view)
        iter_switch.clicked.connect(self.iter_switch_togg)
        iter_switch.setToolTip("Turn on for minor iterations, off for major iterations")
        self._view.setCellWidget(row, 1, iter_switch)

        # Turn off the switch if the x-variable doesn't allow for
        # minor iterations.
        if var_item.var.data_minor is None and var_item.var.name != "iter":
            iter_switch.setEnabled(False)

        self._view.setHorizontalHeaderLabels(["Name", "Major <-> Minor"])
        self._view.resizeColumnsToContents()
        self._view.resizeRowsToContents()

    def iter_switch_togg(self):
        """
        Controls the functionality when the major/minor iteration switch
        is toggled on/off.

        In the off position, only the major iterations are included.
        If the switch is on, we attempt to plot minor iterations unless
        they do not exist for one or more of the x or y-variables.
        """
        switch = self._view.sender()
        x_var = self._plot_model.x_var

        if switch.isChecked():
            flag = False
            for y_var in self._plot_model.y_vars.values():
                if y_var.data_minor is not None:
                    y_var.options["major"] = False
                else:
                    flag = True

            if not flag:
                x_var.options["major"] = False
            else:
                msg_title = "Minor Iterations Warning"
                msg_text = (
                    "One of the y-variables does not support minor iterations.\n\nSwitching back to major iterations."
                )
                QMessageBox.warning(self._view, msg_title, msg_text)
                switch.setChecked(False)

        else:
            switch.setChecked(False)
            x_var.options["major"] = True
            for y_var in self._plot_model.y_vars.values():
                y_var.options["major"] = True

        self._plot_model.plot()
        self._parent_model.canvas.draw()

    def clear_vars(self):
        """
        Clears all the variables in the table view and resets the row
        count to zero.
        """
        self._view.clear()
        self._view.setRowCount(0)


class ConfigurePlotView(QDialog, QStackedWidget, View):
    def __init__(self, parent: QWidget, controller: Controller, name: str):
        """
        The view for the plot configuration window.

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QWidget
            The parent tab view.
        controller : Controller
            The configure plot controller linked to this view.
        name : str
            The name of the window, should be the same as the plot.
        """
        super(ConfigurePlotView, self).__init__(parent)
        self._center()  # centers the application in the middle of the screen
        self.setWindowTitle(name)  # sets the GUI title
        self._controller = controller
        self._controller.view = self
        self.resize(1200, 800)
        self._initView()

        # --- Anything that needs to be done upon re-opening the window ---
        self._controller.setup_var_tables()
        self._controller.populate_files()

    def _initView(self):
        """
        Initializes the view layout.
        """
        # --- Create top layout
        layout = QHBoxLayout()

        # --- Create sub layout for files ---
        left_layout = QVBoxLayout()
        layout.addLayout(left_layout, 1)

        # --- Create sub layout for variables ---
        right_layout = QVBoxLayout()
        layout.addLayout(right_layout, 2)

        # ==============================================================
        # File Management - Top Left Layout
        # ==============================================================
        # --- Add file(s) button ---
        self.add_file_btn = Button("Add file(s)", self)
        self.add_file_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.add_file_btn.clicked.connect(self._controller.add_file)
        left_layout.addWidget(self.add_file_btn)

        # --- File list ---
        self.file_tree = QTreeWidget(self)
        self.file_tree.setColumnCount(1)
        self.file_tree.setHeaderLabels(["File Name"])
        self.file_tree.itemClicked.connect(self._controller.file_selected)
        left_layout.addWidget(self.file_tree)

        # ==============================================================
        # Y-Variables - Top Right Layout
        # ==============================================================
        # --- Add y-vars combobox ---
        self.y_query = QLineEdit(self)
        self.y_query.setPlaceholderText("Search...")
        self.y_query.textChanged.connect(self._controller.y_var_search)
        right_layout.addWidget(self.y_query)

        # --- Create button layout for y-variables ---
        right_button_layout = QHBoxLayout()
        right_layout.addLayout(right_button_layout)

        # --- Add selected variables button ---
        self.add_sel_btn = Button("Add Selected Variables")
        self.add_sel_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.add_sel_btn.clicked.connect(self._controller.add_selected_vars)
        right_button_layout.addWidget(self.add_sel_btn)

        # --- Remove selected variables button ---
        self.rem_sel_btn = Button("Remove Selected Variables")
        self.rem_sel_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.rem_sel_btn.clicked.connect(self._controller.rem_selected_vars)
        right_button_layout.addWidget(self.rem_sel_btn)

        # --- Add y-vars variable list ---
        self.y_table = YTableWidget(self)
        right_layout.addWidget(self.y_table, 5)

        # ==============================================================
        # X-Variables - Middle Right Layout
        # ==============================================================
        self.x_cbox = ExtendedComboBox(self)
        self.x_cbox.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.x_cbox.setToolTip("Type to search for variables")
        self.x_cbox.activated.connect(self._controller.add_x_var)
        right_layout.addWidget(self.x_cbox)

        self.x_table = XTableWidget(self)
        right_layout.addWidget(self.x_table)

        # --- Set the main layout ---
        self.setLayout(layout)

    def _center(self):
        """
        Centers the window on the screen.
        """
        qr = self.frameGeometry()
        self.move(qr.center())


class ConfigureModel(Model):
    def __init__(self):
        """
        Model for the configure view and controller.
        """
        super(ConfigureModel, self).__init__()


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

        self._xtable_controller.view = self._view.x_table
        self._ytable_controller.view = self._view.y_table

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

        self._plot_model.plot()
        self._parent_model.canvas.draw()

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
        sel_items = table.findItems(s, Qt.MatchFlag.MatchContains)

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
                    item.var.set_label(str(label))
                self._plot_model.add_var(item.var, "y")

                item.row_color = GREEN
                self._ytable_controller.set_row_color(item.row())

        self._plot_model.plot()
        self._parent_model.canvas.draw()

    def rem_selected_vars(self):
        items = self._view.y_table.selectedItems()
        for item in items:
            if item.column() == 0:
                self._plot_model.remove_var(item.var, "y")

                self._view.y_table.cellWidget(item.row(), 3).setChecked(False)
                self._view.y_table.cellWidget(item.row(), 4).setChecked(False)

                item.row_color = item.default_row_color
                self._ytable_controller.set_row_color(item.row())

        # Update the plot
        self._plot_model.plot()
        self._parent_model.canvas.draw()
