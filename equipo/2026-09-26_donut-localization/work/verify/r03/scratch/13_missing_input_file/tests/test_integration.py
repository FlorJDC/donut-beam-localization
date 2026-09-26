# -*- coding: utf-8 -*-
"""Integration test: the real glue path patterns -> beams -> photons.make_model -> fisher /
closed_forms / estimators / montecarlo, against a local COPY of the formulas of
``tests/test_acceptance.py::_crb_center`` (copied, not imported)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import beams, closed_forms as cf, estimators as est, fisher, montecarlo as mc
from donutloc import patterns, photons


# --- local copy of the acceptance-test reference (LG donut + TCP, no background) -------------
def _lg(x, y, fwhm):
    r2 = x * x + y * y
    a = 4.0 * np.log(2.0) / fwhm ** 2
    return np.e * a * r2 * np.exp(-a * r2)


def _tcp(L):
    ang = np.pi / 2 + 2 * np.pi * np.arange(3) / 3
    ring = np.stack([0.5 * L * np.cos(ang), 0.5 * L * np.sin(ang)], axis=1)
    return np.vstack([ring, [[0.0, 0.0]]])


def _probs(x, y, L, fwhm):
    c = _tcp(L)
    lam = np.array([_lg(x - cx, y - cy, fwhm) for cx, cy in c])
    return lam / lam.sum()


def _crb(x, y, L, N, fwhm, h=1e-4):
    p = _probs(x, y, L, fwhm)
    dpx = (_probs(x + h, y, L, fwhm) - _probs(x - h, y, L, fwhm)) / (2 * h)
    dpy = (_probs(x, y + h, L, fwhm) - _probs(x, y - h, L, fwhm)) / (2 * h)
    F = N * np.array([[np.sum(dpx * dpx / p), np.sum(dpx * dpy / p)],
                      [np.sum(dpy * dpx / p), np.sum(dpy * dpy / p)]])
    cov = np.linalg.inv(F)
    return float(np.sqrt(0.5 * (cov[0, 0] + cov[1, 1])))


def _crb_center(L, N, fwhm, r=1e-3):
    vals = [_crb(r * np.cos(a), r * np.sin(a), L, N, fwhm, h=r * 1e-2)
            for a in np.linspace(0, 2 * np.pi, 12, endpoint=False)]
    return float(np.mean(vals))


def _model(L, fwhm=300.0, sbr=None):
    return photons.make_model(patterns.tcp_centers(L), beams.make_beam("donut", fwhm=fwhm), sbr=sbr)


class TestCenterLimitThreeRoutes(unittest.TestCase):

    def test_package_closed_form_and_acceptance_copy_agree(self):
        for L, fwhm in ((50.0, 300.0), (100.0, 300.0), (50.0, 360.0)):
            p = _model(L, fwhm)
            num = fisher.crb_limit(p, 100)
            ref = _crb_center(L, 100, fwhm)
            cl = cf.crb_tcp_center_limit(L, 100, fwhm)
            self.assertAlmostEqual(num / ref, 1.0, delta=1e-6, msg=(L, fwhm))
            self.assertAlmostEqual(num / cl, 1.0, delta=1e-6, msg=(L, fwhm))
            # crb_map (default zero_policy="limit") reads the same number at the origin
            m = fisher.crb_map(p, [-1.0, 0.0, 1.0], [0.0], 100)
            self.assertAlmostEqual(m[0, 1] / ref, 1.0, delta=1e-6)
        self.assertAlmostEqual(fisher.crb_limit(_model(50.0), 100), 1.605096, delta=1e-5)
        self.assertAlmostEqual(float(fisher.crb(_model(50.0), [0.0, 0.0], 100)), 1.802472, delta=1e-5)


class TestEstimatorsOnMakeModel(unittest.TestCase):

    def test_mle_efficient_at_center_sbr10(self):
        p = _model(50.0, sbr=10.0)
        res = mc.run_mc(lambda C: est.mle(C, p, search_radius=50.0), p, [0.0, 0.0], 100, 1000, seed=42)
        crb = fisher.crb_limit(p, 100)
        self.assertAlmostEqual(crb / cf.crb_tcp_center_point(50.0, 100, 300.0, sbr=10.0), 1.0, delta=1e-6)
        eff = res["sigma"] / crb
        self.assertTrue(0.9 <= eff <= 1.1, eff)
        self.assertTrue(np.all(np.abs(res["bias"]) < 4 * res["std"] / np.sqrt(1000)))

    def test_lms_general_equals_lms_tcp(self):
        for sbr in (None, 10.0):
            p = _model(50.0, sbr=sbr)
            counts = 1e6 * p(np.array([0.5, -0.3]))              # expected counts
            a = est.lms(counts, p)
            b = est.lms_tcp(counts, 50.0, 300.0, sbr=sbr)
            np.testing.assert_allclose(a, b, atol=1e-8)
            rng = np.random.default_rng(42)
            C = rng.multinomial(200, p(np.array([2.0, 1.0])), size=50)
            np.testing.assert_allclose(est.lms(C, p), est.lms_tcp(C, 50.0, 300.0, sbr=sbr), atol=1e-7)


if __name__ == "__main__":
    unittest.main()
