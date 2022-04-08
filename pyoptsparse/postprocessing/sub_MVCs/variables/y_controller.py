# External modules
from PyQt5 import QtCore
from PyQt5.QtWidgets import QCheckBox, QLineEdit, QTableWidget, QTableWidgetItem

# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Controller, Model
from pyoptsparse.postprocessing.utils.button import Button
from pyoptsparse.postprocessing.utils.colors import BLUE, GREEN
from pyoptsparse.postprocessing.utils.data_structures import File, Variable
from pyoptsparse.postprocessing.utils.widgets import VarTableWidgetItem


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

    def set_view(self, view: QTableWidget):
        """
        Sets the view for the controller.

        Parameters
        ----------
        view : PyQt5.QtWidgets.QTableWidget
            The view for this controller.
        """
        self._view = view

    def populate_vars(self, current_file: File):
        """
        Populates the y-variable tabe view with all the y-variables
        in the current file.

        Parameters
        ----------
        current_file : File
            The current file selected by the user.
        """
        for var in current_file.y_vars:
            self.add_row(var)

        self._view.sortItems(0, QtCore.Qt.AscendingOrder)

        # Find all variables that are in the plot, highlight them in
        # green, and link the table variable to the plot variable
        for row in range(self._view.rowCount()):
            var_item = self._view.item(row, 0)
            self.set_alternating_color(var_item, row)
            if any(var_item.getVar() == y_var for y_var in self._plot_model.y_vars):
                var_item.setRowColor(GREEN)
                self.set_row_color(row)

    def set_alternating_color(self, item, row):
        if row > 0:
            item_prev = self._view.item(row - 1, 0)
            color_prev = item_prev.getRowColor()

            if color_prev == BLUE:
                if item.getVar().name == item_prev.getVar().name:
                    item.setRowColor(BLUE)
                    self.set_row_color(row)

            else:
                if item.getVar().name != item_prev.getVar().name:
                    item.setRowColor(BLUE)
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
        var_item.setVar(var)
        self._view.setItem(row, 0, var_item)

        idx_item = QTableWidgetItem(f"{var_item.getVar().idx}")
        idx_item.setFlags(QtCore.Qt.ItemIsEnabled)
        self._view.setItem(row, 1, idx_item)

        label_item = QTableWidgetItem()
        label_item.setFlags(QtCore.Qt.ItemIsEnabled)
        label = QLineEdit(self._view)
        self._view.setItem(row, 2, label_item)
        self._view.setCellWidget(row, 2, label)

        scaled_opt_item = QTableWidgetItem()
        scaled_opt_item.setFlags(QtCore.Qt.ItemIsEnabled)
        scaled_opt_chbx = QCheckBox(self._view)
        scaled_opt_chbx.setChecked(False)
        scaled_opt_chbx.stateChanged.connect(self.scale_opt_set)
        self._view.setItem(row, 3, scaled_opt_item)
        self._view.setCellWidget(row, 3, scaled_opt_chbx)

        bounds_opt_item = QTableWidgetItem()
        bounds_opt_item.setFlags(QtCore.Qt.ItemIsEnabled)
        bounds_opt_chbx = QCheckBox(self._view)
        bounds_opt_chbx.setChecked(False)
        bounds_opt_chbx.stateChanged.connect(self.bounds_opt_set)
        self._view.setItem(row, 4, bounds_opt_item)
        self._view.setCellWidget(row, 4, bounds_opt_chbx)

        add_btn = Button("Add", self._view)
        add_btn.setFocusPolicy(QtCore.Qt.NoFocus)
        add_item = QTableWidgetItem()
        add_item.setFlags(QtCore.Qt.ItemIsEnabled)
        add_btn.clicked.connect(self.add_var_to_plot)
        self._view.setItem(row, 5, add_item)
        self._view.setCellWidget(row, 5, add_btn)

        rem_btn = Button("Remove", self._view)
        rem_btn.setFocusPolicy(QtCore.Qt.NoFocus)
        rem_item = QTableWidgetItem()
        rem_item.setFlags(QtCore.Qt.ItemIsEnabled)
        rem_btn.clicked.connect(self.remove_var_from_plot)
        self._view.setItem(row, 6, rem_item)
        self._view.setCellWidget(row, 6, rem_btn)

        self._view.setHorizontalHeaderLabels(["Name", "Index", "Label", "Scaled", "Bounds", "Add", "Remove"])
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

        selected_item.getVar().options["scaled"] = scaled_opt

        self._plot_model.plot()
        self._parent_model.canvas.draw()

    def bounds_opt_set(self):
        """
        Sets the bounds options for the selected variable and re-plots
        """
        checkbox = self._view.sender()
        index = self._view.indexAt(checkbox.pos())
        selected_item = self._view.item(index.row(), 0)
        bounds_opt = checkbox.checkState()

        selected_item.getVar().options["bounds"] = bounds_opt

        self._plot_model.plot()
        self._parent_model.canvas.draw()

    def add_var_to_plot(self):
        """
        Adds a y-variable to the plot model and re-plots.
        """
        button = self._view.sender()
        index = self._view.indexAt(button.pos())
        selected_item = self._view.item(index.row(), 0)
        var = selected_item.getVar()
        label = self._view.cellWidget(index.row(), 2).text()
        if str(label) != "":
            var.set_label(str(label))

        self._plot_model.add_var(var, "y")

        self._plot_model.plot()
        self._parent_model.canvas.draw()

        selected_item.setRowColor(GREEN)
        self.set_row_color(index.row())

    def remove_var_from_plot(self):
        """
        Removes a y-variable from the plot model and re-plots
        """
        button = self._view.sender()
        index = self._view.indexAt(button.pos())
        selected_item = self._view.item(index.row(), 0)
        self._plot_model.remove_var(selected_item.getVar(), "y")

        self._view.cellWidget(index.row(), 3).setChecked(False)
        self._view.cellWidget(index.row(), 4).setChecked(False)

        # Update the plot
        self._plot_model.plot()
        self._parent_model.canvas.draw()

        selected_item.setRowColor(selected_item.getDefaultRowColor())
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
        color = self._view.item(row, 0).getRowColor()
        for j in range(self._view.columnCount()):
            self._view.item(row, j).setBackground(color)
