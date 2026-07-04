"""Shared in-widget Help box (so beginners can learn each widget in place)."""
from AnyQt.QtCore import QUrl
from AnyQt.QtGui import QDesktopServices

from Orange.widgets import gui

TUTORIAL_URL = "https://tai-shengyeh.github.io/spectraview/assay.html"


def add_help(widget, text: str, anchor: str = "") -> None:
    """Add a compact 'How to use / 說明' box at the top of a widget's controls."""
    box = gui.widgetBox(widget.controlArea, "ℹ 說明 How to use")
    label = gui.label(box, widget, text)
    try:
        label.setWordWrap(True)
    except Exception:
        pass
    url = TUTORIAL_URL + (f"#{anchor}" if anchor else "")
    gui.button(box, widget, "📖 開啟線上教學 Open tutorial",
               callback=lambda: QDesktopServices.openUrl(QUrl(url)))
