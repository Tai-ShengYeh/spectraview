"""Synthetic demo spectra so the app is useful immediately, with no data files."""
from __future__ import annotations

import numpy as np

from .spectrum import Spectrum


def _gauss(x, center, amp, width):
    return amp * np.exp(-0.5 * ((x - center) / width) ** 2)


def _bands(x, peaks):
    y = np.zeros_like(x)
    for c, a, w in peaks:
        y += _gauss(x, c, a, w)
    return y


def demo_ftir() -> Spectrum:
    """A synthetic mid-IR absorbance spectrum with a sloping baseline + noise."""
    x = np.linspace(4000, 400, 1800)
    bands = _bands(x, [(3400, 0.45, 120), (2920, 0.6, 30), (2850, 0.4, 25),
                       (1745, 0.8, 18), (1640, 0.3, 30), (1460, 0.35, 20),
                       (1160, 0.5, 22), (720, 0.25, 15)])
    baseline = 0.05 + 0.00003 * (4000 - x)
    rng = np.random.default_rng(1)
    y = bands + baseline + rng.normal(0, 0.004, x.size)
    return Spectrum(x=x, y=y, name="Demo FTIR (oil)", x_unit="cm-1",
                    y_unit="absorbance")


def demo_raman() -> Spectrum:
    """A synthetic Raman spectrum (sharp bands on a fluorescence background)."""
    x = np.linspace(200, 3200, 1600)
    bands = _bands(x, [(1003, 9000, 6), (1300, 4000, 25), (1450, 5000, 18),
                       (1660, 7000, 14), (2850, 6000, 22), (2930, 8000, 20)])
    fluor = 3000 * np.exp(-((x - 2200) ** 2) / (2 * 1400 ** 2)) + 500
    rng = np.random.default_rng(2)
    y = bands + fluor + rng.normal(0, 80, x.size)
    return Spectrum(x=x, y=y, name="Demo Raman", x_unit="raman_cm-1",
                    y_unit="counts", meta={"laser_nm": 785.0})


def demo_uvvis() -> Spectrum:
    """A synthetic UV/Vis absorbance spectrum."""
    x = np.linspace(200, 800, 1200)
    bands = _bands(x, [(270, 0.9, 18), (320, 0.5, 30), (520, 0.7, 40)])
    rng = np.random.default_rng(3)
    y = bands + 0.02 + rng.normal(0, 0.003, x.size)
    return Spectrum(x=x, y=y, name="Demo UV/Vis", x_unit="nm", y_unit="absorbance")


def demo_nir() -> Spectrum:
    """A synthetic NIR reflectance-style spectrum with scatter offset."""
    x = np.linspace(1100, 2500, 1400)
    bands = _bands(x, [(1450, 0.25, 40), (1730, 0.3, 35), (1940, 0.4, 30),
                       (2100, 0.2, 45), (2310, 0.35, 25)])
    rng = np.random.default_rng(4)
    y = 0.4 + 0.0002 * (x - 1100) + bands + rng.normal(0, 0.004, x.size)
    return Spectrum(x=x, y=y, name="Demo NIR", x_unit="nm", y_unit="absorbance")


def demo_xrf() -> Spectrum:
    """A synthetic XRF spectrum (keV) with common food-mineral lines."""
    x = np.linspace(1.0, 20.0, 2200)

    def line(c, a, w=0.06):
        return a * np.exp(-(x - c) ** 2 / (2 * w ** 2))

    peaks = (line(3.314, 2000) + line(3.692, 5000)     # K, Ca
             + line(6.404, 8000) + line(7.058, 1500)    # Fe Kα, Kβ
             + line(8.048, 3000) + line(8.639, 2500))    # Cu, Zn
    rng = np.random.default_rng(9)
    y = peaks + 300 * np.exp(-x / 4.0) + 50 + rng.normal(0, 28, x.size)
    return Spectrum(x=x, y=y, name="Demo XRF (minerals)", x_unit="keV", y_unit="counts")


def demo_eem():
    """A synthetic fluorescence EEM: two fluorophores + a Rayleigh scatter ridge."""
    from .eem import EEM
    ex = np.linspace(240.0, 450.0, 44)
    em = np.linspace(280.0, 560.0, 130)
    EX, EM = np.meshgrid(ex, em, indexing="ij")   # (n_ex, n_em)

    def blob(ex0, em0, a, sx, sy):
        return a * np.exp(-((EX - ex0) ** 2 / (2 * sx ** 2) + (EM - em0) ** 2 / (2 * sy ** 2)))

    Z = blob(275, 330, 1000, 14, 24) + blob(350, 440, 720, 20, 34)   # two fluorophores
    Z += 1600 * np.exp(-((EM - EX) ** 2) / (2 * 6.0 ** 2))           # Rayleigh ridge (em≈ex)
    Z += 350 * np.exp(-((EM - 2 * EX) ** 2) / (2 * 7.0 ** 2))        # 2nd-order Rayleigh
    Z += 300.0
    return EEM(ex=ex, em=em, Z=Z, name="Demo EEM")


def demo_eem_stack(n_samples: int = 7, seed: int = 0):
    """A stack of EEMs for PARAFAC: 3 fluorophores at varying concentrations."""
    from .eem import EEM
    ex = np.linspace(240.0, 450.0, 40)
    em = np.linspace(280.0, 560.0, 90)
    EX, EM = np.meshgrid(ex, em, indexing="ij")

    def blob(ex0, em0, sx, sy):
        return np.exp(-((EX - ex0) ** 2 / (2 * sx ** 2) + (EM - em0) ** 2 / (2 * sy ** 2)))

    profiles = [blob(275, 330, 14, 24), blob(350, 440, 20, 34), blob(300, 400, 16, 28)]
    rng = np.random.default_rng(seed)
    conc = rng.uniform(0.2, 1.0, size=(n_samples, 3))
    eems = []
    for i in range(n_samples):
        Z = sum(conc[i, f] * profiles[f] for f in range(3))
        Z = Z + rng.normal(0, 0.004, Z.shape)            # mild noise
        eems.append(EEM(ex=ex, em=em, Z=np.clip(Z, 0, None), name=f"sample {i + 1}"))
    return eems


def demo_cos_series(n: int = 14):
    """A perturbation series for 2D-COS.

    Three bands change at DIFFERENT points along the perturbation (sigmoidal),
    so both the synchronous and the asynchronous maps show clear cross-peaks:
    the 1200 band falls early, 1350 mid, 1500 rises late.
    """
    x = np.linspace(1000.0, 1800.0, 400)

    def band(c, a, w):
        return a * np.exp(-((x - c) ** 2) / (2 * w ** 2))

    def sig(t, c, k=14.0):
        return 1.0 / (1.0 + np.exp(-k * (t - c)))

    specs = []
    for j in range(n):
        t = j / (n - 1)
        y = ((1 - sig(t, 0.35)) * band(1200, 1.0, 25)
             + sig(t, 0.50) * band(1350, 0.6, 20)
             + sig(t, 0.68) * band(1500, 0.9, 30))
        specs.append(Spectrum(x=x, y=y, name=f"t={t:.2f}", x_unit="cm-1",
                              y_unit="absorbance", meta={"perturbation": t}))
    return specs


def load_demo_set():
    """Return a small variety pack of demo spectra."""
    return [demo_ftir(), demo_raman(), demo_uvvis(), demo_nir(), demo_xrf()]
