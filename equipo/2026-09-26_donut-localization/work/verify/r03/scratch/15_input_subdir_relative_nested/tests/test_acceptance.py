"""Acceptance test -- the definition of "done" for this project.

Written BEFORE the team started, and hash-pinned: the team may make it pass, never edit it.
It deliberately re-computes the key physics with its OWN minimal implementation (independent of
src/donutloc), so the reported numbers cannot pass by construction.

Run:  python -m unittest tests.test_acceptance -v
"""

import json
import os
import subprocess
import sys
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NUMBERS = os.path.join(ROOT, "data", "paper_numbers.json")
FIGURES = os.path.join(ROOT, "structure", "figures.json")


def _numbers():
    with open(NUMBERS, encoding="utf-8") as fh:
        data = json.load(fh)
    return {k: (v["value"] if isinstance(v, dict) else v) for k, v in data.items()}


# --- independent reference implementation (LG donut + TCP, no background) ---------------------

def _lg(x, y, fwhm):
    r2 = x * x + y * y
    a = 4.0 * np.log(2.0) / fwhm ** 2
    return np.e * a * r2 * np.exp(-a * r2)          # = 4 e ln2 r^2/fwhm^2 exp(-4 ln2 r^2/fwhm^2)


def _tcp(L):
    ang = np.pi / 2 + 2 * np.pi * np.arange(3) / 3
    ring = np.stack([0.5 * L * np.cos(ang), 0.5 * L * np.sin(ang)], axis=1)
    return np.vstack([ring, [[0.0, 0.0]]])


def _probs(x, y, L, fwhm):
    c = _tcp(L)
    lam = np.array([_lg(x - cx, y - cy, fwhm) for cx, cy in c])
    return lam / lam.sum()


def _crb_center(L, N, fwhm, r=1e-3):
    """At the exact center the central exposure has p = 0 (perfect zero, no background), so
    the CRB there is defined as the limit r -> 0; averaged over directions of approach (for the
    symmetric TCP it is isotropic to < 1e-4 relative)."""
    vals = [_crb(r * np.cos(a), r * np.sin(a), L, N, fwhm, h=r * 1e-2)
            for a in np.linspace(0, 2 * np.pi, 12, endpoint=False)]
    return float(np.mean(vals))


def _crb(x, y, L, N, fwhm, h=1e-4):
    p = _probs(x, y, L, fwhm)
    dpx = (_probs(x + h, y, L, fwhm) - _probs(x - h, y, L, fwhm)) / (2 * h)
    dpy = (_probs(x, y + h, L, fwhm) - _probs(x, y - h, L, fwhm)) / (2 * h)
    F = N * np.array([[np.sum(dpx * dpx / p), np.sum(dpx * dpy / p)],
                      [np.sum(dpy * dpx / p), np.sum(dpy * dpy / p)]])
    cov = np.linalg.inv(F)
    return float(np.sqrt(0.5 * (cov[0, 0] + cov[1, 1])))


class Acceptance(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(NUMBERS):
            raise unittest.SkipTest("data/paper_numbers.json does not exist yet")
        cls.n = _numbers()

    def test_fwhm_convention_is_sane(self):
        self.assertTrue(200 <= self.n["fwhm_nm"] <= 500, self.n["fwhm_nm"])

    def test_crb_at_tcp_center_matches_independent_computation(self):
        ref = _crb_center(50.0, 100, self.n["fwhm_nm"])
        got = self.n["crb_center_lg_L50_N100_nm"]
        self.assertAlmostEqual(got / ref, 1.0, delta=0.01, msg=f"reported {got}, independent {ref}")

    def test_crb_scales_linearly_with_L_and_as_inverse_sqrt_N(self):
        self.assertAlmostEqual(self.n["crb_exponent_L"], 1.0, delta=0.05)
        self.assertAlmostEqual(self.n["crb_exponent_N"], -0.5, delta=0.01)

    def test_mle_is_efficient_at_the_center(self):
        eff = self.n["mle_efficiency_center"]          # std(MLE) / CRB, Monte Carlo
        self.assertTrue(0.9 <= eff <= 1.2, eff)

    def test_vectorial_donut_zero_depth_depends_on_handedness(self):
        self.assertLess(self.n["vectorial_zero_depth_correct"], 1e-3)
        self.assertGreater(self.n["vectorial_zero_depth_wrong_handedness"], 1e-2)

    def test_iterative_minflux_beats_camera_at_equal_photons(self):
        self.assertLess(self.n["iterative_sigma_nm"], self.n["camera_sigma_nm"])

    def test_every_declared_figure_exists_with_its_script(self):
        with open(FIGURES, encoding="utf-8") as fh:
            figs = json.load(fh)
        self.assertGreaterEqual(len(figs), 6, "expected at least six figures (R1-R5 + schematic)")
        for f in figs:
            self.assertTrue(os.path.exists(os.path.join(ROOT, f["pdf"])), f["pdf"])
            self.assertTrue(os.path.exists(os.path.join(ROOT, f["script"])), f["script"])

    def test_manuscript_provenance_passes(self):
        proc = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "check_provenance.py")],
                              cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertIn("all checks pass", proc.stdout)


if __name__ == "__main__":
    unittest.main()
