# External modules
from PyQt6.QtCore import QPropertyAnimation, QRectF, QSize, QSortFilterProxyModel, Qt, pyqtProperty
from PyQt6.QtGui import QColor, QDropEvent, QKeySequence, QPainter, QShortcut
from PyQt6.QtWidgets import (
    QAbstractButton,
    QComboBox,
    QCompleter,
    QHBoxLayout,
    QLineEdit,
    QListWidget,
    QPushButton,
    QSizePolicy,
    QTableWidgetItem,
    QTreeWidgetItem,
    QWidget,
)

# First party modules
from pyoptsparse.postprocessing.baseclasses import Controller, View
from pyoptsparse.postprocessing.data_structures import File

# --- Color constants ---
GREEN = QColor(0, 255, 0, 20)
BLUE = QColor(0, 150, 255, 20)
WHITE = QColor(255, 255, 255)


class Button(QPushButton):
    def __init__(self, *args, **kwargs):
        """
        Inherits the PyQt6 push button class and implements a custom
        button format.
        """
        super().__init__(*args, **kwargs)
        self.resize(150, 30)


class ExtendedComboBox(QComboBox):
    def __init__(self, parent: QWidget = None):
        """
        Combobox view with custom autocomplete functionality for storing and
        searching for variables

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QWidget, optional
            The parent view, by default None
        """
        super(ExtendedComboBox, self).__init__(parent)

        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setEditable(True)

        # add a filter model to filter matching items
        self.pFilterModel = QSortFilterProxyModel(self)
        self.pFilterModel.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.pFilterModel.setSourceModel(self.model())

        # add a completer, which uses the filter model
        self.completer = QCompleter(self.pFilterModel, self)
        # always show all (filtered) completions
        self.completer.setCompletionMode(QCompleter.CompletionMode.UnfilteredPopupCompletion)
        self.setCompleter(self.completer)

        # connect signals
        self.lineEdit().textEdited.connect(self.pFilterModel.setFilterFixedString)
        self.completer.activated.connect(self.on_completer_activated)

    # on selection of an item from the completer, select the
    # corresponding item from combobox
    def on_completer_activated(self, text):
        if text:
            index = self.findText(text)
            self.setCurrentIndex(index)

    # on model change, update the models of the filter and completer
    # as well
    def setModel(self, model):
        super(ExtendedComboBox, self).setModel(model)
        self.pFilterModel.setSourceModel(model)
        self.completer.setModel(self.pFilterModel)

    # on model column change, update the model column of the filter
    # and completer as well
    def setModelColumn(self, column):
        self.completer.setCompletionColumn(column)
        self.pFilterModel.setFilterKeyColumn(column)
        super(ExtendedComboBox, self).setModelColumn(column)


class Switch(QAbstractButton):
    def __init__(self, parent: QWidget = None, track_radius: int = 10, thumb_radius: int = 8):
        """
        On/off slider switch Widget.
        Source code found at:
        https://stackoverflow.com/questions/14780517/toggle-switch-in-qt/51023362

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QWidget
            The parent view, default = None
        track_radius : int
            The radius of the slider track, default = 10
        thumb_radius : int
            The radius of the thumb slider, default = 8
        """
        super().__init__(parent=parent)
        self.setCheckable(True)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

        self._track_radius = track_radius
        self._thumb_radius = thumb_radius

        self._margin = max(0, self._thumb_radius - self._track_radius)
        self._base_offset = max(self._thumb_radius, self._track_radius)
        self._end_offset = {
            True: lambda: self.width() - self._base_offset,
            False: lambda: self._base_offset,
        }
        self._offset = self._base_offset

        palette = self.palette()
        if self._thumb_radius > self._track_radius:
            self._track_color = {
                True: palette.highlight(),
                False: palette.dark(),
            }
            self._thumb_color = {
                True: palette.highlight(),
                False: palette.light(),
            }
            self._text_color = {
                True: palette.highlightedText().color(),
                False: palette.dark().color(),
            }
            self._thumb_text = {
                True: "",
                False: "",
            }
            self._track_opacity = 0.5
        else:
            self._thumb_color = {
                True: palette.highlightedText(),
                False: palette.light(),
            }
            self._track_color = {
                True: palette.highlight(),
                False: palette.dark(),
            }
            self._text_color = {
                True: palette.highlight().color(),
                False: palette.dark().color(),
            }
            self._thumb_text = {
                True: "\u2713",
                False: "\u2716",
            }
            self._track_opacity = 1

    @pyqtProperty(int)
    def offset(self):
        """
        Returns the offset.
        """
        return self._offset

    @offset.setter
    def offset(self, value):
        """
        Sets the offset.

        Parameters
        ----------
        value
            Offset value.
        """
        self._offset = value
        self.update()

    def sizeHint(self):
        """
        Returns the size.
        """
        return QSize(
            4 * self._track_radius + 2 * self._margin,
            2 * self._track_radius + 2 * self._margin,
        )

    def setChecked(self, checked: bool):
        """
        Sets the widget to the checked state.

        Parameters
        ----------
        checked : bool
            True for checked, False for unchecked.
        """
        super().setChecked(checked)
        self.offset = self._end_offset[checked]()

    def resizeEvent(self, event):
        """
        Handles resizing the widget.

        Parameters
        ----------
        event
            The event that triggered the resize call.
        """
        super().resizeEvent(event)
        self.offset = self._end_offset[self.isChecked()]()

    def paintEvent(self, event):
        """
        Handles painting the widget (setting the style).

        Parameters
        ----------
        event
            The event that triggered the paint call.
        """
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.setPen(Qt.PenStyle.NoPen)
        track_opacity = self._track_opacity
        thumb_opacity = 1.0
        text_opacity = 1.0
        if self.isEnabled():
            track_brush = self._track_color[self.isChecked()]
            thumb_brush = self._thumb_color[self.isChecked()]
            text_color = self._text_color[self.isChecked()]
        else:
            track_opacity *= 0.8
            track_brush = self.palette().shadow()
            thumb_brush = self.palette().mid()
            text_color = self.palette().shadow().color()

        p.setBrush(track_brush)
        p.setOpacity(track_opacity)
        p.drawRoundedRect(
            self._margin,
            self._margin,
            self.width() - 2 * self._margin,
            self.height() - 2 * self._margin,
            self._track_radius,
            self._track_radius,
        )
        p.setBrush(thumb_brush)
        p.setOpacity(thumb_opacity)
        p.drawEllipse(
            self.offset - self._thumb_radius,
            self._base_offset - self._thumb_radius,
            2 * self._thumb_radius,
            2 * self._thumb_radius,
        )
        p.setPen(text_color)
        p.setOpacity(text_opacity)
        font = p.font()
        font.setPixelSize(1.5 * self._thumb_radius)
        p.setFont(font)
        p.drawText(
            QRectF(
                self.offset - self._thumb_radius,
                self._base_offset - self._thumb_radius,
                2 * self._thumb_radius,
                2 * self._thumb_radius,
            ),
            Qt.AlignmentFlag.AlignCenter,
            self._thumb_text[self.isChecked()],
        )

    def mouseReleaseEvent(self, event):
        """
        Specifies what should happen when the mouse is released while
        using the slider.

        Parameters
        ----------
        event
            The mouse release event that triggered the call.
        """
        super().mouseReleaseEvent(event)
        if event.button() == Qt.MouseButton.LeftButton:
            anim = QPropertyAnimation(self, b"offset", self)
            anim.setDuration(120)
            anim.setStartValue(self.offset)
            anim.setEndValue(self._end_offset[self.isChecked()]())
            anim.start()

    def enterEvent(self, event):
        """
        Specifies what should happen when the enter key is pressed while
        using the slider.

        Parameters
        ----------
        event
            The enter key event that triggered the call.
        """
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        super().enterEvent(event)


class PlotList(QListWidget, View):
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


class PlotListWidget(View):
    def __init__(self, parent: QWidget = None, controller: Controller = None, idx: int = 0):
        """
        Custom list widget that adds functionality for configuring and
        removing plots.  The widget needs access to the tab window
        controller so the embedded buttons can call functions for
        configuring and removing plots.

        Parameters
        ----------
        parent : PyQt6.QtWidgets.QWidget, optional
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


class FileTableWidgetItem(QTableWidgetItem):
    def __init__(self, *args, **kwargs):
        """
        Custom file table widget item.
        """
        super(FileTableWidgetItem, self).__init__(*args, **kwargs)
