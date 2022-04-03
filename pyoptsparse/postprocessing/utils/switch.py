# --- Python 3.8 ---

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from PyQt5.QtCore import QPropertyAnimation, QRectF, QSize, Qt, pyqtProperty
from PyQt5.QtGui import QPainter
from PyQt5.QtWidgets import QAbstractButton, QSizePolicy, QWidget

# ==============================================================================
# Extension modules
# ==============================================================================


class Switch(QAbstractButton):
    def __init__(self, parent: QWidget = None, track_radius: int = 10, thumb_radius: int = 8):
        """
        On/off slider switch Widget.
        Source code found at:
        https://stackoverflow.com/questions/14780517/toggle-switch-in-qt/51023362

        Parameters
        ----------
        parent : PyQt5.QtWidgets.QWidget
            The parent view, default = None
        track_radius : int
            The radius of the slider track, default = 10
        thumb_radius : int
            The radius of the thumb slider, default = 8
        """
        super().__init__(parent=parent)
        self.setCheckable(True)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

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
        p.setRenderHint(QPainter.Antialiasing, True)
        p.setPen(Qt.NoPen)
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
            Qt.AlignCenter,
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
        if event.button() == Qt.LeftButton:
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
        self.setCursor(Qt.PointingHandCursor)
        super().enterEvent(event)
