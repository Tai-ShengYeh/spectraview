"""Application entry point: configure Qt + pyqtgraph and show the main window."""
from __future__ import annotations

import os
import sys

# Make sure pyqtgraph binds to PySide6 (not a stray PyQt install).
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PySide6")

import pyqtgraph as pg  # noqa: E402
from PySide6 import QtWidgets  # noqa: E402

from . import __app_name__  # noqa: E402
from .ui import MainWindow  # noqa: E402


def configure_pyqtgraph() -> None:
    pg.setConfigOptions(antialias=True, background="w", foreground="k",
                        imageAxisOrder="row-major")


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv if argv is None else argv)
    configure_pyqtgraph()

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(argv)
    app.setApplicationName(__app_name__)

    win = MainWindow()
    win.show()

    # Load any spectra passed on the command line.
    files = [a for a in argv[1:] if os.path.isfile(a)]
    if files:
        win._load_paths(files)

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
