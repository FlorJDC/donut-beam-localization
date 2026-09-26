# -*- coding: utf-8 -*-
"""Tests for donutloc.closed_forms against an INLINE numerical Fisher matrix (Eq. S11).

The numerical reference is self-contained on purpose (it does not import donutloc.fisher,
beams, patterns or photons): TCP with L = diameter, centre last, LG donut of Eq. S17 with
optional residual zero eps (constant / gaussian), multiphoton power and Eq. S30 background.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import closed_forms as cf


# ---------------- inline reference model + numerical Fisher --------------------------------

def _centers(L):
    ang = np.pi / 2 + 2 * np.pi * np.arange(3) / 3
    ring = np.stack([0.5 * L * np.cos(ang), 0.5 * L * np.sin(ang)], axis=1)
    return np.vstack([ring, [[0.0, 0.0]]])


def _probs(r, L, fwhm, eps=0.0, zero_model="constant", sbr=np.inf, power=1, sbr_ref="beam"):
    d2 = ((np.asarray(r, float)[None, :] - _centers(L)) ** 2).sum(axis=1)
    if np.isinf(fwhm):
        lg = d2.copy()
        g = np.ones(4)
    else:
        a = 4 * np.log(2) / fwhm ** 2
        g = np.exp(-a * d2)
        lg = np.e * a * d2 * g
    lam = lg + (eps if zero_model == "constant" else eps * g)
    lam = lam ** power
    if not np.isinf(sbr):
        sig = lam.sum() if sbr_ref == "beam" else lg.sum()
        lam = lam + sig / (4.0 * sbr)          # Eq. S28-S30 with equal background per exposure
    return lam / lam.sum()


def _crb_num(r, L, N, fwhm, h=1e-3, p_min=1e-12, **kw):
    r = np.asarray(r, float)
    p = _probs(r, L, fwhm, **kw)
    grad = np.array([(_probs(r + h * e, L, fwhm, **kw) - _probs(r - h * e, L, fwhm, **kw)) / (2 * h)
                     for e in np.eye(2)])
    m = p > p_min                                   # point convention: drop p_i = 0 terms
    F = N * np.einsum("ik,jk->ij", grad[:, m] / p[m], grad[:, m])
    cov = np.linalg.inv(F)
    return float(np.sqrt(0.5 * np.trace(cov))), cov


def _crb_lim_num(L, N, fwhm, r0=1e-3, n_dir=12, **kw):
    vals = [_crb_num(r0 * np.array([np.cos(t), np.sin(t)]), L, N, fwhm, h=r0 * 1e-2, **kw)[0]
            for t in np.linspace(0, 2 * np.pi, n_dir, endpoint=False)]
    return float(np.mean(vals)), vals


LS = (50.0, 100.0)
FWHMS = (300.0, 360.0, np.inf)


class TestPointValue(unittest.TestCase):
    def test_point_matches_numerical_fisher(self):
        for L in LS:
            for f in FWHMS:
                num = _crb_num([0, 0], L, 100, f)[0]
                cl = cf.crb_tcp_center_point(L, 100, f)
                self.assertAlmostEqual(num / cl, 1.0, delta=1e-6, msg=(L, f, num, cl))

    def test_point_is_eq_s27(self):
        L, N, f = 50.0, 100, 300.0
        s27 = L / (2 * np.sqrt(2 * N)) / (1 - L ** 2 * np.log(2) / f ** 2)
        self.assertAlmostEqual(cf.crb_tcp_center_point(L, N, f), s27, places=12)

    def test_point_is_isotropic(self):
        cov = _crb_num([0, 0], 100.0, 100, 300.0)[1]
        self.assertAlmostEqual(cov[0, 0] / cov[1, 1], 1.0, delta=1e-8)
        self.assertLess(abs(cov[0, 1]) / cov[0, 0], 1e-8)


class TestLimit(unittest.TestCase):
    def test_limit_matches_numerical_average(self):
        for L in LS:
            for f in FWHMS:
                num, _ = _crb_lim_num(L, 100, f)
                cl = cf.crb_tcp_center_limit(L, 100, f)
                self.assertAlmostEqual(num / cl, 1.0, delta=1e-3, msg=(L, f, num, cl))

    def test_limit_is_direction_independent(self):
        _, vals = _crb_lim_num(100.0, 100, 300.0)
        self.assertLess((max(vals) - min(vals)) / np.mean(vals), 1e-4)

    def test_quadratic_ratio_is_two_over_sqrt5(self):
        self.assertAlmostEqual(cf.limit_to_point_ratio(50.0, np.inf), 2 / np.sqrt(5), places=14)
        self.assertAlmostEqual(cf.crb_tcp_center_limit(100.0, 100) ** 2, 0.1 * 100.0 ** 2 / 100, places=10)

    def test_finite_fwhm_ratios(self):
        # ratios seen numerically at fwhm = 300: 0.8944 (L=5), 0.8905 (L=50), 0.8780 (L=100)
        for L, want in ((5.0, 0.8944), (50.0, 0.8905), (100.0, 0.8780)):
            self.assertAlmostEqual(cf.limit_to_point_ratio(L, 300.0), want, delta=6e-5)
            num = _crb_lim_num(L, 100, 300.0)[0] / _crb_num([0, 0], L, 100, 300.0)[0]
            self.assertAlmostEqual(num, cf.limit_to_point_ratio(L, 300.0), delta=1e-4)

    def test_limit_axes(self):
        spar, sperp = cf.crb_tcp_center_limit_axes(100.0, 100)
        cov = _crb_num([1e-3, 0], 100.0, 100, np.inf, h=1e-5)[1]
        self.assertAlmostEqual(np.sqrt(cov[0, 0]) / spar, 1.0, delta=1e-4)
        self.assertAlmostEqual(np.sqrt(cov[1, 1]) / sperp, 1.0, delta=1e-4)
        self.assertAlmostEqual(spar, 2.739, delta=5e-4)            # note A, check 12
        self.assertAlmostEqual(np.sqrt(0.5 * (spar ** 2 + sperp ** 2)),
                               cf.crb_tcp_center_limit(100.0, 100), places=12)


class TestBackground(unittest.TestCase):
    def test_eq_s31_point_and_continuity(self):
        for L in LS:
            for f in (300.0, 360.0):
                for sbr in (5.0, 10.0):
                    cl = cf.crb_tcp_center_point(L, 100, f, sbr=sbr)
                    num = _crb_num([0, 0], L, 100, f, sbr=sbr)[0]
                    self.assertAlmostEqual(num / cl, 1.0, delta=1e-6)
                    lim = _crb_lim_num(L, 100, f, sbr=sbr)[0]
                    self.assertAlmostEqual(lim / cl, 1.0, delta=1e-4)   # no discontinuity

    def test_non_commuting_limits(self):
        L, N, f = 50.0, 100, 300.0
        s27 = cf.crb_tcp_center_point(L, N, f)
        lim0 = cf.crb_tcp_center_limit(L, N, f)
        # order 1: r -> 0 at fixed SBR, then SBR -> inf  => POINT value (Eq. S27)
        big = cf.crb_tcp_center_point(L, N, f, sbr=1e4)
        self.assertLess(big / s27 - 1.0, 2e-4)
        self.assertGreater(big / lim0, 1.1)
        rc = cf.crossover_radius(L, f, 0.0, 1e4)
        self.assertGreater(rc, 100 * 1e-3)                 # r0 = 1e-3 nm << r_c
        num = _crb_lim_num(L, N, f, sbr=1e4)[0]            # numerically r -> 0 at SBR = 1e4
        self.assertAlmostEqual(num / big, 1.0, delta=1e-4)
        # order 2: SBR -> inf at fixed small r (r0 >> r_c), then r -> 0  => LIMIT value
        num2 = _crb_lim_num(L, N, f, sbr=1e14)[0]
        self.assertAlmostEqual(num2 / lim0, 1.0, delta=1e-3)


class TestEps(unittest.TestCase):
    def test_eps_matches_numerical(self):
        for zm in ("constant", "gaussian"):
            for eps in (0.01, 0.1):
                for L in LS:
                    for f in (300.0, 360.0):
                        for sbr in (np.inf, 5.0):
                            for ref in ("beam", "lg"):
                                cl = cf.crb_tcp_center_eps(L, 100, f, eps, sbr, zm, ref)
                                num = _crb_num([0, 0], L, 100, f, eps=eps, zero_model=zm,
                                               sbr=sbr, sbr_ref=ref)[0]
                                self.assertAlmostEqual(num / cl, 1.0, delta=1e-6,
                                                       msg=(zm, eps, L, f, sbr, ref))

    def test_eps_is_continuous_at_origin(self):
        cl = cf.crb_tcp_center_eps(100.0, 100, 300.0, 0.01, zero_model="gaussian")
        lim = _crb_lim_num(100.0, 100, 300.0, eps=0.01, zero_model="gaussian")[0]
        self.assertAlmostEqual(lim / cl, 1.0, delta=1e-4)

    def test_constant_pedestal_is_effective_sbr(self):
        for ref in ("beam", "lg"):
            for sbr in (np.inf, 5.0, 20.0):
                for eps in (0.002, 0.05):
                    se = cf.sbr_eff_constant_pedestal(100.0, 300.0, eps, sbr, ref)
                    a = cf.crb_tcp_center_point(100.0, 100, 300.0, sbr=se)
                    b = cf.crb_tcp_center_eps(100.0, 100, 300.0, eps, sbr, "constant", ref)
                    self.assertAlmostEqual(a / b, 1.0, delta=1e-12)

    def test_total_alias(self):
        a = cf.crb_tcp_center_eps(100.0, 100, 300.0, 0.03, 5.0, "gaussian", "total")
        b = cf.crb_tcp_center_eps(100.0, 100, 300.0, 0.03, 5.0, "gaussian", "beam")
        self.assertEqual(a, b)
        with self.assertRaises(ValueError):
            cf.crb_tcp_center_eps(100.0, 100, 300.0, 0.03, 5.0, "gaussian", "bad")

    def test_eps_zero_gives_point_value(self):
        for zm in ("constant", "gaussian"):
            self.assertAlmostEqual(cf.crb_tcp_center_eps(50.0, 100, 300.0, 0.0, zero_model=zm),
                                   cf.crb_tcp_center_point(50.0, 100, 300.0), places=12)

    def test_degradation_monotonic(self):
        eps = np.array([0.0, 0.002, 0.01, 0.03, 0.05, 0.1, 0.15])
        for zm in ("constant", "gaussian"):
            v = cf.crb_tcp_center_eps(100.0, 100, 300.0, eps, zero_model=zm)
            self.assertTrue(np.all(np.diff(v) > 0))


class TestNoteAChecks(unittest.TestCase):
    def test_numbers(self):
        self.assertAlmostEqual(cf.crb_tcp_center_point(100, 100), 3.536, delta=5e-4)
        self.assertAlmostEqual(cf.crb_tcp_center_point(100, 100, 360), 3.735, delta=5e-4)
        self.assertAlmostEqual(cf.crb_tcp_center_point(50, 100, 360), 1.792, delta=5e-4)
        self.assertAlmostEqual(cf.crb_tcp_center_point(100, 500, 360, sbr=10), 1.817, delta=5e-4)
        for L, w360, w300 in ((50, 0.941, 0.947), (100, 1.962, 2.012), (150, 3.167, 3.370)):
            self.assertAlmostEqual(cf.crb_tcp_center_point(L, 500, 360, sbr=5), w360, delta=5e-4)
            self.assertAlmostEqual(cf.crb_tcp_center_point(L, 500, 300, sbr=5), w300, delta=5e-4)
        self.assertAlmostEqual(cf.crb_1d_center(300, 100, 300, "gaussian"), 10.82, delta=5e-3)
        self.assertAlmostEqual(cf.crb_1d_center(50, 100, kind="quadratic"), 1.25, places=12)
        self.assertAlmostEqual(cf.crb_1d_center(50, 100, np.inf, "donut"), 1.25, places=12)
        self.assertAlmostEqual(cf.sbr_center_vs_L(50, 100, 10, 360), 2.60, delta=5e-3)
        self.assertAlmostEqual(cf.sbr_center_vs_L(25, 100, 10, 360), 0.66, delta=5e-3)
        for c, w in ((1, 3.536), (2, 1.768), (3, 1.179)):
            self.assertAlmostEqual(cf.crb_tcp_center_point(100, 100, power=c), w, delta=5e-4)

    def test_multiphoton_no_discontinuity(self):
        for c in (2, 3):
            lim = _crb_lim_num(100.0, 100, 300.0, power=c)[0]
            self.assertAlmostEqual(lim / cf.crb_tcp_center_limit(100.0, 100, 300.0, power=c), 1.0,
                                   delta=1e-3)
            pt = _crb_num([0, 0], 100.0, 100, 300.0, power=c)[0]
            self.assertAlmostEqual(pt / cf.crb_tcp_center_point(100.0, 100, 300.0, power=c), 1.0,
                                   delta=1e-6)

    def test_1d_donut_matches_numerical(self):
        L, N, f = 100.0, 100, 300.0
        a = 4 * np.log(2) / f ** 2

        def p0(x):
            i0 = np.e * a * (x + L / 2) ** 2 * np.exp(-a * (x + L / 2) ** 2)
            i1 = np.e * a * (x - L / 2) ** 2 * np.exp(-a * (x - L / 2) ** 2)
            return i0 / (i0 + i1)
        h = 1e-3
        d = (p0(h) - p0(-h)) / (2 * h)
        num = np.sqrt(p0(0) * (1 - p0(0))) / abs(d) / np.sqrt(N)       # Eq. S21
        self.assertAlmostEqual(num / cf.crb_1d_center(L, N, f), 1.0, delta=1e-6)


class TestApi(unittest.TestCase):
    def test_vectorised(self):
        L = np.array([50.0, 100.0, 150.0])
        out = cf.crb_tcp_center_limit(L, 100, 300.0)
        self.assertEqual(out.shape, (3,))
        for i, l in enumerate(L):
            self.assertAlmostEqual(out[i], cf.crb_tcp_center_limit(l, 100, 300.0), places=12)
        self.assertIsInstance(cf.crb_tcp_center_point(50, 100), float)

    def test_errors(self):
        with self.assertRaises(ValueError):
            cf.crb_tcp_center_point(400.0, 100, 300.0)        # x >= 1
        with self.assertRaises(ValueError):
            cf.crb_tcp_center_eps(50.0, 100, np.inf, 0.01)
        with self.assertRaises(ValueError):
            cf.crb_tcp_center_eps(50.0, 100, 300.0, 0.01, zero_model="bad")
        with self.assertRaises(ValueError):
            cf.crb_tcp_center_eps(50.0, 100, 300.0, -0.01)
        with self.assertRaises(ValueError):
            cf.crb_1d_center(50.0, 100, kind="bad")
        with self.assertRaises(ValueError):
            cf.crb_1d_center(50.0, 100, kind="gaussian")
        with self.assertRaises(ValueError):
            cf.crb_tcp_center_point(50.0, 100, power=0.5)


if __name__ == "__main__":
    unittest.main()
