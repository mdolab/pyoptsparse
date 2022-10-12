# External modules
from PyQt6.QtWidgets import QDialog, QTableWidget

# First party modules
from pyoptsparse.postprocessing.baseclasses.view import View


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
