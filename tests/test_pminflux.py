# -*- coding: utf-8 -*-
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import beams, estimators, fisher, patterns, photons, pminflux as pm


class TestCrosstalkMatrix(unittest.TestCase):

    def test_closed_form_vs_monte_carlo(self):
        for order in [(0, 1, 2, 3), (3, 2, 0, 1)]:
            for tau in (2.0, 5.0, 20.0):
                M = pm.crosstalk_matrix(tau, order=order)
                Mh = pm.crosstalk_mc(tau, order=order, n=200000, rng=42)
                se = np.sqrt(M * (1 - M) / 200000)
                self.assertTrue(np.all(np.abs(M - Mh) <= 5 * se + 1e-12),
                                (order, tau, np.abs(M - Mh).max()))

    def test_stochastic(self):
        for tau in (0.0, 0.5, 3.0, 1e3):
            M = pm.crosstalk_matrix(tau, order=(1, 3, 0, 2))
            np.testing.assert_allclose(M.sum(axis=0), 1.0, atol=1e-13)
            np.testing.assert_allclose(M.sum(axis=1), 1.0, atol=1e-13)
            self.assertTrue(np.all(M >= 0))

    def test_limits(self):
        np.testing.assert_array_equal(pm.crosstalk_matrix(0.0), np.eye(4))
        np.testing.assert_allclose(pm.crosstalk_matrix(1e-3), np.eye(4), atol=1e-300)
        # tau -> infinity: uniform mixing
        np.testing.assert_allclose(pm.crosstalk_matrix(1e7), 0.25, atol=1e-6)

    def test_explicit_entries_and_order(self):
        tau, T = 4.0, 12.5
        q = np.exp(-T / tau)
        c = (1 - q) * q ** np.arange(4) / (1 - q ** 4)
        M = pm.crosstalk_matrix(tau, T=T, order=(0, 1, 2, 3))
        # photons of exposure 2 (slot 2) leak first into the window of exposure 3 (slot 3)
        self.assertAlmostEqual(M[3, 2], c[1], places=14)
        self.assertAlmostEqual(M[0, 2], c[2], places=14)
        self.assertAlmostEqual(M[1, 2], c[3], places=14)
        # centre fired first: its photons leak into exposure 0
        M2 = pm.crosstalk_matrix(tau, T=T, order=(3, 0, 1, 2))
        self.assertAlmostEqual(M2[0, 3], c[1], places=14)
        # a cyclic shift of the order gives the same matrix
        np.testing.assert_allclose(M, M2)
        M3 = pm.crosstalk_matrix(tau, T=T, order=(0, 2, 1, 3))
        self.assertFalse(np.allclose(M, M3))

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            pm.crosstalk_matrix(-1.0)
        with self.assertRaises(ValueError):
            pm.crosstalk_matrix(1.0, order=(0, 1, 1, 3))
        with self.assertRaises(ValueError):
            pm.crosstalk_matrix(1.0, T=0.0)


class TestCrosstalkModel(unittest.TestCase):

    def setUp(self):
        self.c = patterns.tcp_centers(100.0)
        self.beam = beams.make_beam("donut", fwhm=300.0)

    def test_model_is_M_p(self):
        r = np.array([[0.0, 0.0], [20.0, -5.0], [40.0, 10.0]])
        f = pm.make_model_crosstalk(self.c, self.beam, 3.0)
        p = photons.probabilities(r, self.c, self.beam)
        np.testing.assert_allclose(f(r), p @ f.M.T, atol=1e-15)
        np.testing.assert_allclose(f(r).sum(-1), 1.0, atol=1e-13)
        self.assertEqual(f.K, 4)

    def test_uniform_background_invariant(self):
        f = pm.make_model_crosstalk(self.c, self.beam, 5.0, sbr=0.0)
        np.testing.assert_allclose(f(np.array([10.0, 3.0])), 0.25, atol=1e-14)

    def test_crb_centre_no_zero(self):
        # with cross-talk the centre exposure is never zero: point value == limit
        f = pm.make_model_crosstalk(self.c, self.beam, 3.0)
        a = fisher.crb(f, np.zeros(2), 100.0)
        b = fisher.crb(f, np.zeros(2), 100.0, zero_policy="limit")
        self.assertAlmostEqual(float(a), float(b), places=6)
        ideal = fisher.crb(photons.make_model(self.c, self.beam), np.zeros(2), 100.0,
                           zero_policy="limit")
        self.assertGreater(float(a), float(ideal))

    def test_mirror_symmetry_between_orders(self):
        # x -> -x swaps exposures 1 and 2 of the TCP; the pulse orders (0,1,2,3) and (0,2,1,3)
        # are mirror images, so p_A(x, y) equals p_B(-x, y) with exposures 1 and 2 swapped.
        fA = pm.make_model_crosstalk(self.c, self.beam, 4.0, order=(0, 1, 2, 3))
        fB = pm.make_model_crosstalk(self.c, self.beam, 4.0, order=(0, 2, 1, 3))
        r = np.array([[17.0, 4.0], [-30.0, 12.0]])
        rm = r * np.array([-1.0, 1.0])
        np.testing.assert_allclose(fA(r), fB(rm)[:, [0, 2, 1, 3]], atol=1e-14)

    def test_mle_with_model_recovers_noise_free(self):
        f = pm.make_model_crosstalk(self.c, self.beam, 5.0)
        r0 = np.array([12.0, -7.0])
        est = estimators.mle(1e6 * f(r0), f, 60.0, tol=1e-5)
        np.testing.assert_allclose(est, r0, atol=1e-3)


class TestFlicker(unittest.TestCase):

    def test_telegraph_duty_and_bounds(self):
        a, b = pm.telegraph_on_intervals(4000, 400.0, 100.0, 100.0, rng=42)
        self.assertTrue(np.all(b >= a) and np.all(a >= 0) and np.all(b <= 400.0))
        frac = (b - a).sum(axis=1) / 400.0
        self.assertAlmostEqual(frac.mean(), 0.5, delta=0.03)
        a, b = pm.telegraph_on_intervals(10, 400.0, 1e9, 1.0, rng=1, start="on")
        np.testing.assert_allclose((b - a).sum(axis=1), 400.0)

    def test_sequential_blocks(self):
        # one on-interval [30, 130) with t_total = 400, 4 blocks of 100
        a = np.array([[30.0]])
        b = np.array([[130.0]])
        np.testing.assert_allclose(pm.sequential_on_times(a, b, 400.0, 1), [[70.0, 30.0, 0, 0]])
        # 2 repetitions: blocks of 50: [0,50)->0 ... [200,250)->0
        np.testing.assert_allclose(pm.sequential_on_times(a, b, 400.0, 2),
                                   [[20.0, 50.0, 30.0, 0.0]])
        # order (3, 0, 1, 2): first block belongs to exposure 3
        np.testing.assert_allclose(pm.sequential_on_times(a, b, 400.0, 1, order=(3, 0, 1, 2)),
                                   [[30.0, 0.0, 0.0, 70.0]])

    def test_interleaved_pulse_count_brute_force(self):
        rng = np.random.default_rng(42)
        P, T, tt = 0.05, 0.0125, 20.0
        a, b = pm.telegraph_on_intervals(30, tt, 2.0, 3.0, rng=rng)
        got = pm.interleaved_on_pulses(a, b, tt, period=P, T=T, order=(2, 0, 3, 1))
        n_per = int(round(tt / P))
        for m in range(30):
            for s, j in enumerate((2, 0, 3, 1)):
                t = s * T + P * np.arange(n_per)
                on = np.any((t[:, None] >= a[m]) & (t[:, None] < b[m]), axis=1)
                self.assertEqual(got[m, j], int(on.sum()))

    def test_interleaved_weights_nearly_equal(self):
        w = pm.exposure_weights("interleaved", 2000, rng=42)
        tot = w.sum(axis=1)
        g = tot > 0.2      # at least ~10 us on
        rel = w[g] / w[g].mean(axis=1, keepdims=True)
        self.assertLess(np.abs(rel - 1).max(), 0.01)
        ws = pm.exposure_weights("sequential", 2000, reps=1, rng=42)
        self.assertGreater(np.std(ws[g] / ws[g].mean(axis=1, keepdims=True)), 0.3)
        np.testing.assert_allclose(w.mean(axis=0), 1.0, atol=0.1)

    def test_simulate_counts(self):
        p = np.array([0.2, 0.3, 0.5, 0.0])
        w = np.ones((5000, 4))
        c = pm.simulate_flicker_counts(p, w, n_mean=100.0, rng=42)
        self.assertAlmostEqual(c.sum(axis=1).mean(), 100.0, delta=1.0)
        self.assertTrue(np.all(c[:, 3] == 0))
        c = pm.simulate_flicker_counts(p, np.vstack([w[:3], np.zeros((1, 4))]), fixed_N=100, rng=1)
        np.testing.assert_array_equal(c.sum(axis=1), [100, 100, 100, 0])
        with self.assertRaises(ValueError):
            pm.exposure_weights("bogus", 3)


if __name__ == "__main__":
    unittest.main()
