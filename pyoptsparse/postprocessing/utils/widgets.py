# External modules
from PyQt5.QtGui import QDropEvent, QKeySequence
from PyQt5.QtWidgets import (
    QHBoxLayout,
    QLineEdit,
    QListWidget,
    QPushButton,
    QShortcut,
    QTableWidgetItem,
    QTreeWidgetItem,
    QWidget,
)

# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Controller
from pyoptsparse.postprocessing.utils.data_structures import File, Variable


class PlotList(QListWidget):
    def __init__(self, parent: QWidget = None, controller: Controller = None) -> None:
        super(PlotList, self).__init__(parent)

        self.controller = controller

        self.plot_up_action = QShortcut(QKeySequence("Ctrl+Up"), self)
        self.plot_down_action = QShortcut(QKeySequence("Ctrl+Down"), self)

        self.plot_up_action.activated.connect(self.movePlotUp)
        self.plot_down_action.activated.connect(self.movePlotDown)

    def dropEvent(self, event: QDropEvent) -> None:
        super(PlotList, self).dropEvent(event)
        self.controller.reorder_plots()

    def movePlotUp(self) -> None:
        self.controller.move_plot_up()

    def movePlotDown(self) -> None:
        self.controller.move_plot_down()


class PlotListWidget(QWidget):
    def __init__(self, parent: QWidget = None, controller: Controller = None, idx: int = 0):
        """
        Custom list widget that adds functionality for configuring and
        removing plots.  The widget needs access to the tab window
        controller so the embedded buttons can call functions for
        configuring and removing plots.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget, optional
            The parent view, by default None
        controller : Controller, optional
            Tab controller linked to the tab view., by default None
        idx : int, optional
            The index of the plot in the tab model, by default 0
        """
        super(PlotListWidget, self).__init__(parent)

        self.controller = controller
        self.idx = idx

        # Set the plot title
        self.title = QLineEdit(f"Plot {idx}")

        # Add the configure plot button
        self.configure_button = QPushButton("Configure/Add Variables")
        self.configure_button.clicked.connect(self.configure)

        # Add the remove plot button
        self.remove_button = QPushButton("Remove Plot")
        self.remove_button.clicked.connect(self.remove)

        # Configure the layout
        layout = QHBoxLayout()
        layout.addWidget(self.title, 1)
        layout.addWidget(self.configure_button, 1)
        layout.addWidget(self.remove_button, 1)

        # Set the layout
        self.setLayout(layout)

    def remove(self):
        """
        Calls the remove_plot function in the tab controller.
        """
        self.controller.remove_plot(self.idx)

    def configure(self):
        """
        Calls teh configure_view function in the tab controller.
        """
        self.controller.configure_view(self.idx)


class FileTreeWidgetItem(QTreeWidgetItem):
    def __init__(self, *args, **kwargs):
        """
        Custom tree widget item that can store file objects.
        """
        self.file = None
        super().__init__(*args, **kwargs)

    def setFile(self, file: File):
        self.file = file


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


class FileTableWidgetItem(QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        """
        Custom file table widget item.
        """
        super(FileTableWidgetItem, self).__init__(*args, **kwargs)


class VarTableWidgetItem(QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        """
        Custom variable table widget item.
        """
        super(VarTableWidgetItem, self).__init__(*args, **kwargs)
        self._var = None
        self._default_row_color = self.background()
        self._row_color = None

    def __lt__(self, other):
        if self._var.name == other._var.name:
            return self._var.idx < other._var.idx
        else:
            return self._var.name < other._var.name

    def setVar(self, var: Variable):
        """
        Sets the variable associated with the widget item.

        Parameters
        ----------
        var : Variable
            The variable to be added to the widget item.
        """
        self._var = var

    def getVar(self):
        return self._var

    def setRowColor(self, color):
        self._row_color = color

    def getRowColor(self):
        return self._row_color

    def getDefaultRowColor(self):
        return self._default_row_color
