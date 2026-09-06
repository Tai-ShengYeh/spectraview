"""Shared in-widget Help box (so beginners can learn each widget in place)."""
from AnyQt.QtCore import QUrl
from AnyQt.QtGui import QDesktopServices
from Orange.widgets import gui

TUTORIAL_URL = "https://tai-shengyeh.github.io/spectraview/orange_en.html"


def add_help(widget, text: str, anchor: str = "") -> None:
    """Add a compact 'How to use' box at the top of a widget's controls.

    ``text`` is short English guidance; the button opens the online tutorial
    (optionally at ``#anchor``) in the user's browser.
    """
    box = gui.widgetBox(widget.controlArea, "ℹ How to use")
    label = gui.label(box, widget, text)
    try:
        label.setWordWrap(True)
    except Exception:
        pass
    url = TUTORIAL_URL + (f"#{anchor}" if anchor else "")
    gui.button(box, widget, "📖 Open tutorial",
               callback=lambda: QDesktopServices.openUrl(QUrl(url)))
