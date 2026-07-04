"""Orange widget category: Spectra (光譜)."""

NAME = "Spectra"
DESCRIPTION = "Spectroscopy: import from IRUG/SOPRANO, similarity, library, mixtures."
ICON = "icons/category.svg"
BACKGROUND = "#e8f4f8"
PRIORITY = 500

# F1 / "?" help: the html-index provider fetches this page, reads the links
# under the element with id="widgets" (link text must equal the widget name)
# and opens the matching page (see orange.canvas.help entry point).
WIDGET_HELP_PATH = (
    ("https://tai-shengyeh.github.io/spectraview/widgets/index.html", None),
)
