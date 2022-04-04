# External modules
from PyQt5.QtWidgets import QDialog, QHBoxLayout, QLineEdit, QTableWidget, QTreeWidget, QVBoxLayout, QWidget

# First party modules
from pyoptsparse.postprocessing.utils.base_classes import Controller
from pyoptsparse.postprocessing.utils.button import Button


class MetadataView(QDialog):
    def __init__(self, parent: QWidget, controller: Controller, name: str):
        """
        The view for the displaying metadata options.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget
            The parent tab view.
        controller : Controller
            The metadata view controller linked to this view.
        name : str
            The name of the window.
        """
        super(MetadataView, self).__init__(parent)
        self._center()
        self.setWindowTitle(name)
        self._controller = controller
        self._controller.set_view(self)
        self.resize(1000, 800)
        self._initView()

    def _initView(self):
        """
        Initializes the view.
        """
        # --- Create layout ---
        layout = QHBoxLayout()

        # ==============================================================
        # File List - Left Layout
        # ==============================================================
        # --- Add file button ---
        self.add_file_btn = Button("Add file(s)", self)
        self.add_file_btn.clicked.connect(self._controller.add_file)
        layout.addWidget(self.add_file_btn)

        # --- File list ---
        self.file_tree = QTreeWidget(self)
        self.file_tree.setColumnCount(1)
        self.file_tree.setHeaderLabels(["File Name"])
        self.file_tree.itemClicked.connect(self._controller.file_selected)
        layout.addWidget(self.file_tree)

        # ==============================================================
        # Options Table - Right Layout
        # ==============================================================
        right_layout = QVBoxLayout()
        layout.addLayout(right_layout, 3)

        self.opt_tree = QTreeWidget(self)
        self.opt_tree.setColumnCount(2)
        self.opt_tree.setHeaderLabels(["Name", "Value"])
        right_layout.addWidget(self.opt_tree)

        self.query = QLineEdit()
        self.query.setPlaceholderText("Search...")
        self.query.textChanged.connect(self._controller.search)
        right_layout.addWidget(self.query)

        self.opt_prob_table = QTableWidget(self)
        right_layout.addWidget(self.opt_prob_table)

        self._controller.populate_files()

        self.setLayout(layout)

        self.show()

    def _center(self):
        """
        Centers the view on the screen.
        """
        qr = self.frameGeometry()
        self.move(qr.center())
