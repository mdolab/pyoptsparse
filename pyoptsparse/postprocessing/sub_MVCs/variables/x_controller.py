# Standard Python modules

# External modules
from PyQt5.QtWidgets import QMessageBox, QTableWidget

# Local modules
from pyoptsparse.postprocessing.utils.widgets import VarTableWidgetItem, FileTableWidgetItem, IterSwitchWidget
from pyoptsparse.postprocessing.utils.base_classes import Controller, Model


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

    def set_view(self, view: QTableWidget):
        """
        Sets the view for the controller.

        Parameters
        ----------
        view : PyQt5.QtWidgets.QTableWidget
            The view for this controller.
        """
        self._view = view

    def populate_vars(self):
        """
        Adds the x-variables to the the table view.
        """
        var = self._plot_model.x_var
        # Create a new variable widget item for the tree view
        var_item = VarTableWidgetItem(var.name)
        file_item = FileTableWidgetItem(var.file.name_short)
        var_item.var = var

        self.add_row(file_item, var_item)

    def add_row(self, file_item: FileTableWidgetItem, var_item: VarTableWidgetItem):
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
        self._view.setItem(row, 0, file_item)
        self._view.setItem(row, 1, var_item)

        iter_switch = IterSwitchWidget(row, self._view)
        iter_switch.clicked.connect(self.iter_switch_togg)
        iter_switch.setToolTip("Turn on for minor iterations, off for major iterations")
        self._view.setCellWidget(row, 2, iter_switch)

        self._view.setHorizontalHeaderLabels(["File", "Name", "Minor/Major"])
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
        sender = self._view.sender()
        x_var = self._plot_model.x_var
        if sender.isChecked():
            # Adjust the iteration option for the x-variable
            x_var.options["major"] = False
            for y_var in self._plot_model.y_vars:
                y_var.options["major"] = False
        else:
            x_var.options["major"] = True
            for y_var in self._plot_model.y_vars:
                y_var.options["major"] = True

        flag = self._plot_model.plot()
        if not flag:
            msg = QMessageBox(self._view)
            msg.setWindowTitle("Variable Warning")
            msg.setText("The x variable or y variables don't have minor iterations.\n\nReverting to major iterations.")
            msg.setIcon(QMessageBox.Warning)
            msg.exec_()

            x_var.options["major"] = True
            for y_var in self._plot_model.y_vars:
                y_var.options["major"] = True

            sender.setChecked(False)

            self._plot_model.plot()

        self._parent_model.canvas.draw()

    def clear_vars(self):
        """
        Clears all the variables in the table view and resets the row
        count to zero.
        """
        self._view.clear()
        self._view.setRowCount(0)
