# External modules
from PyQt5.QtWidgets import QDialog, QTableWidget


class YTableWidget(QTableWidget):
    def __init__(self, parent: QDialog = None):
        """
        Defines the y-variable table view.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QDialog, optional
            The parent view, by default None
        """
        super(YTableWidget, self).__init__(parent)
        self.setShowGrid(False)
        self.verticalHeader().setVisible(False)
        self.setColumnCount(6)
        self.setHorizontalHeaderLabels(["File", "Name", "Scaled", "Bounds", "Add", "Remove"])
        self._controller = None

    def setController(self, controller):
        self._controller = controller
