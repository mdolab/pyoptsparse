# External modules
from PyQt5.QtWidgets import QMessageBox, QTableWidget

# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Controller, Model
from pyoptsparse.postprocessing.utils.switch import Switch
from pyoptsparse.postprocessing.utils.widgets import VarTableWidgetItem


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
