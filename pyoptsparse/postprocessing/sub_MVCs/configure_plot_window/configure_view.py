# Standard Python modules

# External modules
from PyQt5.QtWidgets import QWidget, QDialog, QHBoxLayout, QVBoxLayout, QTreeWidget, QLineEdit, QComboBox

# Local modules
from pyoptsparse.postprocessing.utils.combo_box import ExtendedComboBox
from pyoptsparse.postprocessing.utils.button import Button
from pyoptsparse.postprocessing.utils.base_classes import Controller
from pyoptsparse.postprocessing.sub_MVCs.variables.y_view import YTableWidget
from pyoptsparse.postprocessing.sub_MVCs.variables.x_view import XTableWidget


class ConfigurePlotView(QDialog):
    def __init__(self, parent: QWidget, controller: Controller, name: str):
        """
        The view for the plot configuration window.

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget
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
        self._controller.set_view(self)
        self.resize(1000, 800)
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

        # --- Add y-vars variable list ---
        self.y_table = YTableWidget(self)
        right_layout.addWidget(self.y_table)

        # ==============================================================
        # X-Variables - Middle Right Layout
        # ==============================================================
        self.x_cbox = ExtendedComboBox(self)
        self.x_cbox.setInsertPolicy(QComboBox.NoInsert)
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