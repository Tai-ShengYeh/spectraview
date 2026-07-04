"""Offscreen widget smoke tests for orange-assay.

Run:  QT_QPA_PLATFORM=offscreen python tests/test_widgets.py
"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from Orange.data import ContinuousVariable, Domain, Table  # noqa: E402
from Orange.widgets.tests.base import WidgetTest  # noqa: E402

from orangeassay import core  # noqa: E402
from orangeassay.widgets.owcoffeering import OWCoffeeRing  # noqa: E402
from orangeassay.widgets.owdoseresponse import OWDoseResponse  # noqa: E402
from orangeassay.widgets.owmicroplate import OWMicroplate  # noqa: E402


def synthetic_plate(n_rows=3, n_cols=8, cell=36, gap=24):
    pitch = cell + gap
    plate = np.full((n_rows * pitch, n_cols * pitch), 0.08)
    for r in range(n_rows):
        for c in range(n_cols):
            b = 0.5 + 0.45 * (0.12 if c == 0 else
                              0.15 + 0.8 / (1 + np.exp(-(c - 3))))
            y0, x0 = r * pitch + gap // 2, c * pitch + gap // 2
            plate[y0:y0 + cell, x0:x0 + cell] = b
    return core.load_grayscale(plate)


class TestOWCoffeeRing(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWCoffeeRing)
        self.widget._gray = synthetic_plate()

    def test_normalized_output(self):
        self.widget.n_rows, self.widget.n_cols = 3, 8
        self.widget.threshold_idx, self.widget.manual_threshold = 0, 0.3
        self.widget.normalize = True
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.cells)
        self.assertEqual(len(out), 7)                 # 7 treatment columns
        names = [v.name for v in out.domain.attributes]
        self.assertIn("concentration", names)
        self.assertIn("ratio", names)
        self.assertTrue(np.all(np.isfinite(out.get_column("ratio"))))

    def test_raw_output(self):
        self.widget.n_rows, self.widget.n_cols = 3, 8
        self.widget.normalize = False
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.cells)
        self.assertEqual(len(out), 24)

    def test_otsu_threshold_mode(self):
        self.widget.threshold_idx = 1                 # Otsu
        self.widget._recompute()
        self.assertIsNotNone(self.get_output(self.widget.Outputs.cells))

    def test_no_image(self):
        w = self.create_widget(OWCoffeeRing)
        self.assertTrue(w.Information.no_image.is_shown())


class TestOWMicroplate(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWMicroplate)
        self.widget._gray = synthetic_plate(4, 6)

    def test_well_table(self):
        self.widget.plate_idx = list(core.PLATE_FORMATS).index("24 (4×6)")
        self.widget.metric_idx = 0                     # mean
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.wells)
        self.assertEqual(len(out), 24)
        self.assertEqual(str(out.metas[0][0]), "A1")
        self.assertEqual(str(out.metas[-1][0]), "D6")

    def test_count_metric(self):
        self.widget.plate_idx = list(core.PLATE_FORMATS).index("24 (4×6)")
        self.widget.metric_idx = 1                     # thresholded count
        self.widget.threshold_idx = 0                  # Otsu
        self.widget._recompute()
        self.assertEqual(len(self.get_output(self.widget.Outputs.wells)), 24)

    def test_96_format(self):
        self.widget.plate_idx = list(core.PLATE_FORMATS).index("96 (8×12)")
        self.widget._recompute()
        self.assertEqual(len(self.get_output(self.widget.Outputs.wells)), 96)


class TestOWDoseResponse(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWDoseResponse)

    def _table(self):
        conc = np.array([0.01, 0.09, 0.23, 0.46, 0.91, 3.63, 7.27])
        lx = np.log10(conc)
        y = 2.0 / (1 + np.exp(-2.5 * (lx - np.log10(0.5))))
        dom = Domain([ContinuousVariable("concentration"),
                      ContinuousVariable("signal")])
        return Table.from_numpy(dom, np.column_stack([conc, y]))

    def test_fit_outputs(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        self.widget.conc_var = self.widget._data.domain["concentration"]
        self.widget.signal_var = self.widget._data.domain["signal"]
        self.widget.log_x = True
        self.widget._recompute()
        fit = self.get_output(self.widget.Outputs.fit)
        self.assertIsNotNone(fit)
        cols = {v.name for v in fit.domain.attributes}
        self.assertTrue({"EC50", "LOD", "R2"} <= cols)
        self.assertGreater(float(fit.get_column("R2")[0]), 0.99)
        curve = self.get_output(self.widget.Outputs.curve)
        self.assertEqual(len(curve), 200)

    def test_4pl(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        self.widget.conc_var = self.widget._data.domain["concentration"]
        self.widget.signal_var = self.widget._data.domain["signal"]
        self.widget.model_idx = 1                      # 4PL
        self.widget._recompute()
        self.assertIsNotNone(self.get_output(self.widget.Outputs.fit))

    def test_clear(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        self.send_signal(self.widget.Inputs.data, None)
        self.assertIsNone(self.get_output(self.widget.Outputs.fit))


if __name__ == "__main__":
    unittest.main()
