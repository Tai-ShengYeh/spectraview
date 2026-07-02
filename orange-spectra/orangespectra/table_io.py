"""Orange Table ↔ spectrum-dict converters.

Follows the Orange-Spectroscopy add-on convention: each attribute (column) is a
wavelength/wavenumber value (the column name is the number), each row is one
spectrum, and a "name" meta column labels the rows. Tables produced here can be
viewed directly with Orange-Spectroscopy's own Spectra widget.
"""
from __future__ import annotations

import numpy as np
from Orange.data import ContinuousVariable, Domain, StringVariable, Table

from .core import make_spectrum, merge_spectra


def table_from_spectra(spectra: list, x_label: str = "") -> Table:
    """Build a Table (rows = spectra) on a shared x-grid."""
    gx, ys = merge_spectra(spectra)
    attrs = [ContinuousVariable.make(f"{v:.6g}") for v in gx]
    metas = [StringVariable.make("name"), StringVariable.make("source")]
    domain = Domain(attrs, metas=metas)
    X = np.vstack([np.asarray(y, float) for y in ys])
    M = np.array([[s.get("name", f"spectrum {i}"), s.get("source", "")]
                  for i, s in enumerate(spectra)], dtype=object)
    table = Table.from_numpy(domain, X, metas=M)
    table.name = x_label or (spectra[0].get("x_label", "") if spectra else "")
    table.attributes["x_label"] = x_label or (
        spectra[0].get("x_label", "x") if spectra else "x")
    return table


def spectra_from_table(table: Table) -> list:
    """Read rows of a spectral Table back into spectrum dicts.

    Columns whose names parse as numbers form the x-axis; a "name" meta (or any
    first string meta) labels each row. Non-numeric columns are ignored.
    """
    if table is None or len(table) == 0:
        return []
    xs, cols = [], []
    for i, var in enumerate(table.domain.attributes):
        try:
            xs.append(float(str(var.name)))
            cols.append(i)
        except ValueError:
            continue
    if len(xs) < 2:
        raise ValueError(
            "No spectral columns found: column names must be wavelength/"
            "wavenumber values (e.g. output of Import Spectrum URL, or a "
            "matrix CSV with the axis in the header row)."
        )
    x = np.asarray(xs, float)
    X = table.X[:, cols].astype(float)
    name_idx = None
    for j, mv in enumerate(table.domain.metas):
        if isinstance(mv, StringVariable):
            name_idx = j
            if mv.name.lower() == "name":
                break
    x_label = table.attributes.get("x_label", "x") if hasattr(table, "attributes") else "x"
    out = []
    for r in range(len(table)):
        name = (str(table.metas[r, name_idx]) if name_idx is not None
                else f"spectrum {r + 1}")
        out.append(make_spectrum(x, X[r], name=name or f"spectrum {r + 1}",
                                 x_label=x_label))
    return out
