# External modules
from PyQt5.QtWidgets import QDialog, QTableWidget

# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Controller


class XTableWidget(QTableWidget):
    def __init__(self, parent: QDialog = None):
        """
        Defines a x-variable table view.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QDialog, optional
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
