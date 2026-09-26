# -*- coding: utf-8 -*-
"""Tests for donutloc.fisher (numerical Fisher information / CRB)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import fisher


def tcp_model(L, fwhm=300.0, sbr=None, rotation=np.pi / 2):
    """Minimal inline LG-donut + TCP model (centre last), optional SBR (Balzarotti2017 Eq. S30)."""
    a = 4.0 * np.log(2.0) / fwhm ** 2
    ang = rotation + 2 * np.pi * np.arange(3) / 3
    c = np.vstack([0.5 * L * np.stack([np.cos(ang), np.sin(ang)], 1), [[0.0, 0.0]]])

    def p_fn(r):
        d2 = np.sum((np.asarray(r, float)[..., None, :] - c) ** 2, -1)
        lam = np.e * a * d2 * np.exp(-a * d2)
        p = lam / lam.sum(-1, keepdims=True)
        if sbr is not None and np.isfinite(sbr):
            p = sbr / (sbr + 1.0) * p + 1.0 / (sbr + 1.0) / 4.0
        return p
    return p_fn


def s27(L, N, fwhm):
    return L / (2 * np.sqrt(2 * N)) / (1 - L ** 2 * np.log(2) / fwhm ** 2)


# --- re-implementation of tests/test_acceptance.py::_crb_center (not imported) ---------------
def _ref_crb(x, y, p_fn, N, h):
    def P(x, y):
        return p_fn(np.array([x, y]))
    p = P(x, y)
    dpx = (P(x + h, y) - P(x - h, y)) / (2 * h)
    dpy = (P(x, y + h) - P(x, y - h)) / (2 * h)
    F = N * np.array([[np.sum(dpx * dpx / p), np.sum(dpx * dpy / p)],
                      [np.sum(dpy * dpx / p), np.sum(dpy * dpy / p)]])
    cov = np.linalg.inv(F)
    return float(np.sqrt(0.5 * (cov[0, 0] + cov[1, 1])))


def _ref_crb_center(p_fn, N, r=1e-3):
    return float(np.mean([_ref_crb(r * np.cos(a), r * np.sin(a), p_fn, N, h=r * 1e-2)
                          for a in np.linspace(0, 2 * np.pi, 12, endpoint=False)]))


class TestPointValue(unittest.TestCase):

    def test_center_point_equals_S27(self):
        for L in (50.0, 100.0):
            got = float(fisher.crb(tcp_model(L), [0.0, 0.0], 100))
            self.assertAlmostEqual(got / s27(L, 100, 300.0), 1.0, delta=1e-4)

    def test_exact_zero_excluded_even_with_pmin_zero(self):
        got = float(fisher.crb(tcp_model(50.0), [0.0, 0.0], 100, p_min=0.0))
        self.assertAlmostEqual(got / s27(50.0, 100, 300.0), 1.0, delta=1e-4)

    def test_center_with_sbr10_equals_S31(self):
        L, N, fwhm, sbr = 50.0, 100, 300.0, 10.0
        ref = s27(L, N, fwhm) * np.sqrt((1 + 1 / sbr) * (1 + 3 / (4 * sbr)))
        got = float(fisher.crb(tcp_model(L, fwhm, sbr), [0.0, 0.0], N))
        self.assertAlmostEqual(got / ref, 1.0, delta=1e-4)
        # with background the CRB is continuous: limit == point value
        self.assertAlmostEqual(fisher.crb_limit(tcp_model(L, fwhm, sbr), N) / ref, 1.0, delta=1e-4)


class TestLimit(unittest.TestCase):

    def test_limit_matches_acceptance_definition(self):
        for L, fwhm in ((50.0, 300.0), (100.0, 300.0), (50.0, 360.0)):
            p = tcp_model(L, fwhm)
            got = fisher.crb_limit(p, 100)
            ref = _ref_crb_center(p, 100)
            self.assertAlmostEqual(got / ref, 1.0, delta=1e-3)

    def test_limit_over_point_quadratic_is_2_over_sqrt5(self):
        p = tcp_model(50.0, fwhm=1e6)   # quadratic regime
        ratio = fisher.crb_limit(p, 100) / float(fisher.crb(p, [0.0, 0.0], 100))
        self.assertAlmostEqual(ratio, 2 / np.sqrt(5), delta=1e-3)

    def test_limit_is_not_point_value(self):
        p = tcp_model(50.0)
        self.assertLess(fisher.crb_limit(p, 100), 0.95 * float(fisher.crb(p, [0.0, 0.0], 100)))


class TestScaling(unittest.TestCase):

    def test_inverse_sqrt_N_exact(self):
        p = tcp_model(50.0)
        r = np.array([[0.0, 0.0], [7.0, -3.0]])
        np.testing.assert_allclose(fisher.crb(p, r, 400), 0.5 * fisher.crb(p, r, 100), rtol=1e-12)
        self.assertAlmostEqual(fisher.crb_limit(p, 400) / fisher.crb_limit(p, 100), 0.5, delta=1e-10)

    def test_linear_in_L_for_small_L(self):
        a = float(fisher.crb(tcp_model(5.0), [0.0, 0.0], 100))
        b = float(fisher.crb(tcp_model(10.0), [0.0, 0.0], 100))
        self.assertAlmostEqual(b / a, 2.0, delta=0.01)
        a = fisher.crb_limit(tcp_model(5.0), 100)
        b = fisher.crb_limit(tcp_model(10.0), 100)
        self.assertAlmostEqual(b / a, 2.0, delta=0.01)


class TestShapesAndAxes(unittest.TestCase):

    def test_fisher_matrix_shape_symmetry_and_offcenter_value(self):
        p = tcp_model(50.0)
        r = np.random.default_rng(42).uniform(-20, 20, size=(3, 5, 2))
        F = fisher.fisher_matrix(p, r, 100)
        self.assertEqual(F.shape, (3, 5, 2, 2))
        np.testing.assert_allclose(F[..., 0, 1], F[..., 1, 0])
        c = fisher.crb(p, r, 100, h=1e-4)
        self.assertEqual(c.shape, (3, 5))
        ref = _ref_crb(r[1, 2, 0], r[1, 2, 1], p, 100, h=1e-4)
        self.assertAlmostEqual(c[1, 2] / ref, 1.0, delta=1e-8)

    def test_N_broadcasts(self):
        p = tcp_model(50.0)
        r = np.array([[1.0, 2.0], [1.0, 2.0]])
        c = fisher.crb(p, r, np.array([100.0, 400.0]))
        self.assertAlmostEqual(c[1] / c[0], 0.5, delta=1e-12)

    def test_crb_axes(self):
        p = tcp_model(50.0)
        sx, sy, iso = fisher.crb_axes(p, [0.0, 0.0], 100)
        self.assertAlmostEqual(float(sx) / float(sy), 1.0, delta=1e-6)
        self.assertAlmostEqual(float(iso), 1.0, delta=1e-6)
        sx, sy, iso = fisher.crb_axes(p, [15.0, 5.0], 100)
        self.assertLess(float(iso), 0.99)
        c = float(fisher.crb(p, [15.0, 5.0], 100))
        self.assertAlmostEqual(np.sqrt(0.5 * (sx ** 2 + sy ** 2)) / c, 1.0, delta=1e-12)

    def test_crb_map_shape(self):
        xs = np.linspace(-10, 10, 7)
        ys = np.linspace(-5, 5, 4)
        m = fisher.crb_map(tcp_model(50.0), xs, ys, 100)
        self.assertEqual(m.shape, (4, 7))
        self.assertAlmostEqual(m[1, 3] / float(fisher.crb(tcp_model(50.0), [xs[3], ys[1]], 100)), 1.0,
                               delta=1e-12)

    def test_singular_fisher_gives_inf(self):
        flat = lambda r: np.full(np.asarray(r).shape[:-1] + (4,), 0.25)
        self.assertTrue(np.isinf(fisher.crb(flat, [0.0, 0.0], 100)))

    def test_bad_shape_raises(self):
        with self.assertRaises(ValueError):
            fisher.fisher_matrix(tcp_model(50.0), [0.0, 0.0, 0.0], 100)


if __name__ == "__main__":
    unittest.main()
