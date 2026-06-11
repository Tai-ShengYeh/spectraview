#!/usr/bin/env python
"""Launcher for SpectraView.

Usage:
    python run.py                # open the app
    python run.py file1.dx ...   # open the app and load these spectra
"""
import sys

from specview.app import main

if __name__ == "__main__":
    sys.exit(main())
