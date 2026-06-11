"""XRF characteristic X-ray emission lines and peak→element identification.

Energies are in keV (standard reference values, X-ray Data Booklet). This is a
curated set of the lines most used in XRF — K lines for light/mid-Z elements and
L lines for heavy elements (whose K lines fall outside the usual detector range).
"""
from __future__ import annotations

# Each element: symbol, name, Z, and the main emission-line energies (keV).
ELEMENTS = [
    ("Na", "Sodium", 11, {"Ka1": 1.041}),
    ("Mg", "Magnesium", 12, {"Ka1": 1.254}),
    ("Al", "Aluminium", 13, {"Ka1": 1.487, "Kb1": 1.557}),
    ("Si", "Silicon", 14, {"Ka1": 1.740, "Kb1": 1.836}),
    ("P", "Phosphorus", 15, {"Ka1": 2.014, "Kb1": 2.139}),
    ("S", "Sulfur", 16, {"Ka1": 2.308, "Kb1": 2.464}),
    ("Cl", "Chlorine", 17, {"Ka1": 2.622, "Kb1": 2.816}),
    ("K", "Potassium", 19, {"Ka1": 3.314, "Kb1": 3.590}),
    ("Ca", "Calcium", 20, {"Ka1": 3.692, "Kb1": 4.013}),
    ("Sc", "Scandium", 21, {"Ka1": 4.091, "Kb1": 4.461}),
    ("Ti", "Titanium", 22, {"Ka1": 4.511, "Kb1": 4.932}),
    ("V", "Vanadium", 23, {"Ka1": 4.952, "Kb1": 5.427}),
    ("Cr", "Chromium", 24, {"Ka1": 5.415, "Kb1": 5.947}),
    ("Mn", "Manganese", 25, {"Ka1": 5.899, "Kb1": 6.490}),
    ("Fe", "Iron", 26, {"Ka1": 6.404, "Kb1": 7.058}),
    ("Co", "Cobalt", 27, {"Ka1": 6.930, "Kb1": 7.649}),
    ("Ni", "Nickel", 28, {"Ka1": 7.478, "Kb1": 8.265}),
    ("Cu", "Copper", 29, {"Ka1": 8.048, "Kb1": 8.905}),
    ("Zn", "Zinc", 30, {"Ka1": 8.639, "Kb1": 9.572}),
    ("Ga", "Gallium", 31, {"Ka1": 9.252, "Kb1": 10.264}),
    ("Ge", "Germanium", 32, {"Ka1": 9.886, "Kb1": 10.982}),
    ("As", "Arsenic", 33, {"Ka1": 10.544, "Kb1": 11.726}),
    ("Se", "Selenium", 34, {"Ka1": 11.222, "Kb1": 12.496}),
    ("Br", "Bromine", 35, {"Ka1": 11.924, "Kb1": 13.291}),
    ("Rb", "Rubidium", 37, {"Ka1": 13.395, "Kb1": 14.961}),
    ("Sr", "Strontium", 38, {"Ka1": 14.165, "Kb1": 15.836}),
    ("Y", "Yttrium", 39, {"Ka1": 14.958, "Kb1": 16.738}),
    ("Zr", "Zirconium", 40, {"Ka1": 15.775, "Kb1": 17.668}),
    ("Nb", "Niobium", 41, {"Ka1": 16.615, "Kb1": 18.623}),
    ("Mo", "Molybdenum", 42, {"Ka1": 17.479, "Kb1": 19.608}),
    ("Ag", "Silver", 47, {"Ka1": 22.163, "Kb1": 24.942, "La1": 2.984}),
    ("Cd", "Cadmium", 48, {"Ka1": 23.174, "Kb1": 26.096, "La1": 3.134}),
    ("Sn", "Tin", 50, {"Ka1": 25.271, "Kb1": 28.486, "La1": 3.444}),
    ("Sb", "Antimony", 51, {"Ka1": 26.359, "Kb1": 29.726, "La1": 3.605}),
    ("I", "Iodine", 53, {"Ka1": 28.612, "Kb1": 32.295, "La1": 3.938}),
    ("Cs", "Caesium", 55, {"Ka1": 30.973, "Kb1": 34.987, "La1": 4.286}),
    ("Ba", "Barium", 56, {"Ka1": 32.194, "Kb1": 36.378, "La1": 4.466, "Lb1": 4.828}),
    ("La", "Lanthanum", 57, {"La1": 4.651, "Lb1": 5.042}),
    ("Ce", "Cerium", 58, {"La1": 4.840, "Lb1": 5.262}),
    ("Nd", "Neodymium", 60, {"La1": 5.230, "Lb1": 5.722}),
    ("Sm", "Samarium", 62, {"La1": 5.636, "Lb1": 6.205}),
    ("Gd", "Gadolinium", 64, {"La1": 6.057, "Lb1": 6.714}),
    ("Hf", "Hafnium", 72, {"La1": 7.899, "Lb1": 9.023}),
    ("Ta", "Tantalum", 73, {"La1": 8.146, "Lb1": 9.343}),
    ("W", "Tungsten", 74, {"La1": 8.398, "Lb1": 9.672}),
    ("Pt", "Platinum", 78, {"La1": 9.442, "Lb1": 11.071}),
    ("Au", "Gold", 79, {"La1": 9.713, "Lb1": 11.443}),
    ("Hg", "Mercury", 80, {"La1": 9.989, "Lb1": 11.823}),
    ("Tl", "Thallium", 81, {"La1": 10.268, "Lb1": 12.213}),
    ("Pb", "Lead", 82, {"La1": 10.551, "Lb1": 12.614}),
    ("Bi", "Bismuth", 83, {"La1": 10.839, "Lb1": 13.023}),
    ("Th", "Thorium", 90, {"La1": 12.968, "Lb1": 15.624}),
    ("U", "Uranium", 92, {"La1": 13.615, "Lb1": 16.428}),
]

# Pretty names for the line keys.
LINE_NAMES = {"Ka1": "Kα1", "Kb1": "Kβ1", "La1": "Lα1", "Lb1": "Lβ1"}

# Flat search table: (symbol, name, line_key, energy_keV).
LINES = [(sym, name, line, energy)
         for sym, name, _z, d in ELEMENTS
         for line, energy in d.items()]

ENERGY_RANGE = (min(e for *_, e in LINES), max(e for *_, e in LINES))


def identify_energy(energy_kev: float, tol: float = 0.10,
                    line_filter=None):
    """Return element/line matches within ``tol`` keV of ``energy_kev``.

    Each match: dict(symbol, name, line, line_label, energy, delta). Sorted by
    |delta| (closest first).
    """
    matches = []
    for sym, name, line, e in LINES:
        if line_filter and line not in line_filter:
            continue
        delta = energy_kev - e
        if abs(delta) <= tol:
            matches.append({"symbol": sym, "name": name, "line": line,
                            "line_label": LINE_NAMES.get(line, line),
                            "energy": e, "delta": delta})
    matches.sort(key=lambda m: abs(m["delta"]))
    return matches


def identify_peaks(energies, tol: float = 0.10, line_filter=None):
    """Identify a list of peak energies (keV). Returns one result per energy:
    dict(energy, matches=[...], best=<closest match or None>).
    """
    out = []
    for e in energies:
        ms = identify_energy(float(e), tol, line_filter)
        out.append({"energy": float(e), "matches": ms,
                    "best": ms[0] if ms else None})
    return out
