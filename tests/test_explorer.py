# -*- coding: utf-8 -*-
"""Smoke test: every mode of explore/donut_explorer.py runs with the defaults, headless."""
import contextlib
import io
import os
import sys
import tempfile
import unittest

import matplotlib

matplotlib.use("Agg")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "explore"))

import donut_explorer as ex  # noqa: E402


class TestExplorer(unittest.TestCase):
    def test_all_modes_run(self):
        with tempfile.TemporaryDirectory() as d:
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                ex.main(["--all", "--yes", "--no-show", "--save", d])
            pngs = [f for f in os.listdir(d) if f.endswith(".png")]
            self.assertEqual(len(pngs), len(ex.MODES))
        text = out.getvalue()
        self.assertIn("valor puntual (S27) = 1.802 nm", text)
        self.assertIn("STD/CRB = 1.00", text)

    def test_bad_mode(self):
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            ex.main(["--mode", "99", "--yes", "--no-show"])


if __name__ == "__main__":
    unittest.main()
