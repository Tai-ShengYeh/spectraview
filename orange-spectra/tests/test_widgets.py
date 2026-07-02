"""GUI smoke tests for the four Orange widgets (offscreen).

Run:  QT_QPA_PLATFORM=offscreen python -m unittest tests.test_widgets -v
"""
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Orange.widgets.tests.base import WidgetTest  # noqa: E402

from orangespectra import core  # noqa: E402
from orangespectra.table_io import table_from_spectra  # noqa: E402
from orangespectra.widgets import owimporturl, owlibrary  # noqa: E402
from orangespectra.widgets.owimporturl import OWImportSpectrumURL  # noqa: E402
from orangespectra.widgets.owlibrary import OWSpectralLibrary  # noqa: E402
from orangespectra.widgets.owmixture import OWMixtureAnalysis  # noqa: E402
from orangespectra.widgets.owsimilarity import OWSpectraSimilarity  # noqa: E402

X = np.linspace(400, 1800, 200)


def _g(c, w):
    return np.exp(-((X - c) ** 2) / (2.0 * w * w))


def _table(*name_y):
    return table_from_spectra(
        [core.make_spectrum(X, y, name=n, x_label="wavenumber (cm-1)")
         for n, y in name_y])


_IRUG_PAGE = ("<html><script>var jqPlotData={};jqPlotData.series=[{data:["
              + ",".join(f'"{x:.1f}:{y:.4f}"' for x, y in
                         zip(np.linspace(1900, 100, 50), np.cos(np.linspace(0, 6, 50))))
              + "]}];</script></html>").encode()


class TestImportURL(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWImportSpectrumURL)
        self._orig = owimporturl.load_spectrum_url
        owimporturl.load_spectrum_url = lambda u: core.load_spectrum_url(
            u, fetch=lambda _url: (_IRUG_PAGE, "text/html"))

    def tearDown(self):
        owimporturl.load_spectrum_url = self._orig
        super().tearDown()

    def test_fetch_outputs_table(self):
        self.widget.url_edit.setText("4119")
        self.widget.fetch()
        out = self.get_output(self.widget.Outputs.data)
        self.assertIsNotNone(out)
        self.assertEqual(len(out), 1)
        self.assertEqual(len(out.domain.attributes), 50)
        self.assertEqual(str(out.metas[0][0]), "IRUG 4119")
        # second fetch accumulates
        self.widget.fetch()
        out2 = self.get_output(self.widget.Outputs.data)
        self.assertEqual(len(out2), 2)
        self.widget.clear_all()
        self.assertIsNone(self.get_output(self.widget.Outputs.data))

    def test_fetch_error_reported(self):
        owimporturl.load_spectrum_url = lambda u: (_ for _ in ()).throw(
            ValueError("boom"))
        self.widget.url_edit.setText("https://x.org/none")
        self.widget.fetch()
        self.assertTrue(self.widget.Error.fetch_failed.is_shown())


class TestSimilarity(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWSpectraSimilarity)

    def test_data_vs_references(self):
        data = _table(("q", _g(600, 30)))
        refs = _table(("match", _g(600, 30) * 2), ("other", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.data, data)
        self.send_signal(self.widget.Inputs.references, refs)
        out = self.get_output(self.widget.Outputs.scores)
        self.assertEqual(len(out), 2)
        self.assertEqual(str(out.metas[0][1]), "match")   # best first
        self.assertAlmostEqual(float(out.X[0][0]), 1.0, places=6)

    def test_within_data_pairs(self):
        data = _table(("a", _g(600, 30)), ("b", _g(600, 30)), ("c", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.data, data)
        out = self.get_output(self.widget.Outputs.scores)
        self.assertEqual(len(out), 3)                     # C(3,2) pairs

    def test_non_spectral_table_errors(self):
        from Orange.data import Table
        self.send_signal(self.widget.Inputs.data, Table("iris"))
        self.assertTrue(self.widget.Error.bad_table.is_shown())


class TestLibrary(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWSpectralLibrary)

    def test_build_search_save_load(self):
        refs = _table(("R1", _g(600, 30)), ("R2", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.spectra, refs)
        self.widget.add_pending()
        lib_out = self.get_output(self.widget.Outputs.library)
        self.assertEqual(len(lib_out), 2)

        query = _table(("unknown", _g(600, 30) + 0.01))
        self.send_signal(self.widget.Inputs.query, query)
        hits = self.get_output(self.widget.Outputs.hits)
        self.assertEqual(len(hits), 2)
        self.assertEqual(str(hits.metas[0][1]), "R1")     # best match first
        best = self.get_output(self.widget.Outputs.best)
        self.assertEqual(str(best.metas[0][0]), "R1")

        # save via core (file dialog is interactive) and reload through widget API
        path = os.path.join(tempfile.mkdtemp(), "t.speclib")
        core.save_library(self.widget._library, path)
        self.widget.clear()
        self.assertIsNone(self.get_output(self.widget.Outputs.library))
        loaded = core.load_library(path)
        self.assertEqual([e["name"] for e in loaded], ["R1", "R2"])


class TestMixture(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWMixtureAnalysis)

    def test_nnls_composition(self):
        mix = _table(("mix", 3.0 * _g(600, 30) + 1.0 * _g(1200, 40) + 0.05))
        refs = _table(("R1", _g(600, 30)), ("R2", _g(1200, 40)),
                      ("R3", _g(900, 25)))
        self.send_signal(self.widget.Inputs.mixture, mix)
        self.send_signal(self.widget.Inputs.references, refs)
        comp = self.get_output(self.widget.Outputs.composition)
        self.assertEqual(len(comp), 3)
        # sorted by fraction: R1 (3x) first, then R2 (1x)
        self.assertEqual(str(comp.metas[0][0]), "R1")
        self.assertEqual(str(comp.metas[1][0]), "R2")
        self.assertAlmostEqual(float(comp.X[0][0]), 3.0, places=1)
        self.assertGreater(comp.attributes["r_squared"], 0.999)
        fit = self.get_output(self.widget.Outputs.fit)
        self.assertEqual(len(fit), 3)                     # mixture / fit / residual

    def test_no_overlap_error(self):
        mix = _table(("mix", _g(600, 30)))
        x2 = np.linspace(5000, 6000, 50)
        refs = table_from_spectra(
            [core.make_spectrum(x2, np.ones(50), name="far")])
        self.send_signal(self.widget.Inputs.mixture, mix)
        self.send_signal(self.widget.Inputs.references, refs)
        self.assertTrue(self.widget.Error.analysis_failed.is_shown())
        self.assertIsNone(self.get_output(self.widget.Outputs.composition))


if __name__ == "__main__":
    unittest.main()
