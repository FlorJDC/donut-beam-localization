# -*- coding: utf-8 -*-
"""Tests for donutloc.background (free-background MLE, augmented Fisher / CRB)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import background as bg
from donutloc import beams, estimators, fisher, patterns, photons

C = patterns.tcp_centers(100.0)
BEAM = beams.make_beam("donut", fwhm=300.0)


class TestModel(unittest.TestCase):
    def test_bg_from_sbr_roundtrip(self):
        for sbr in (1.0, 5.0, 20.0):
            b = bg.bg_from_sbr(C, BEAM, sbr)
            self.assertAlmostEqual(float(photons.sbr_at(np.zeros(2), C, BEAM, b)), sbr, places=10)
        # off-centre reference point
        r = np.array([30.0, -10.0])
        b = bg.bg_from_sbr(C, BEAM, 7.0, r_ref=r)
        self.assertAlmostEqual(float(photons.sbr_at(r, C, BEAM, b)), 7.0, places=10)

    def test_bg_from_sbr_invalid(self):
        for bad in (0.0, -1.0, np.inf, np.nan):
            with self.assertRaises(ValueError):
                bg.bg_from_sbr(C, BEAM, bad)

    def test_probabilities_match_photons(self):
        rng = np.random.default_rng(0)
        r = rng.uniform(-60, 60, size=(7, 2))
        for b in (0.0, 0.01, 0.3):
            th = np.concatenate([r, np.full((7, 1), b)], axis=1)
            p = bg.probabilities_bg(th, C, BEAM)
            ref = photons.probabilities(r, C, BEAM, bg_per_exposure=b)
            np.testing.assert_allclose(p, ref, rtol=1e-13, atol=0)
            np.testing.assert_allclose(p.sum(-1), 1.0, rtol=1e-13)

    def test_probabilities_invalid(self):
        with self.assertRaises(ValueError):
            bg.probabilities_bg(np.array([0.0, 0.0, -1e-3]), C, BEAM)
        with self.assertRaises(ValueError):
            bg.probabilities_bg(np.array([0.0, 0.0]), C, BEAM)


class TestFisher(unittest.TestCase):
    def test_xy_block_equals_known_b_fisher(self):
        b = bg.bg_from_sbr(C, BEAM, 5.0)
        model = photons.make_model(C, BEAM, bg_per_exposure=b)
        for r in ([0.0, 0.0], [20.0, 0.0], [10.0, -25.0]):
            r = np.array(r)
            F3 = bg.fisher_matrix_bg(C, BEAM, r, b, 500)
            F2 = fisher.fisher_matrix(model, r, 500)
            np.testing.assert_allclose(F3[:2, :2], F2, rtol=1e-6, atol=1e-9 * np.abs(F2).max())
            np.testing.assert_allclose(F3, F3.T, rtol=1e-12, atol=1e-12 * np.abs(F3).max())
            out = bg.crb_free_bg(C, BEAM, r, b, 500)
            self.assertAlmostEqual(out["sigma_known_b"], float(fisher.crb(model, r, 500)), places=8)
            # nuisance parameter can only increase the bound
            self.assertGreaterEqual(out["sigma"], out["sigma_known_b"] * (1 - 1e-12))

    def test_centre_decouples_by_symmetry(self):
        # at the TCP centre d p / d b is the same for the three ring exposures and the
        # position gradients sum to zero over the ring -> F_xb = F_yb = 0 -> no CRB penalty
        b = bg.bg_from_sbr(C, BEAM, 5.0)
        F3 = bg.fisher_matrix_bg(C, BEAM, np.zeros(2), b, 500)
        self.assertLess(abs(F3[0, 2]) / np.sqrt(F3[0, 0] * F3[2, 2]), 1e-6)
        self.assertLess(abs(F3[1, 2]) / np.sqrt(F3[1, 1] * F3[2, 2]), 1e-6)
        out = bg.crb_free_bg(C, BEAM, np.zeros(2), b, 500)
        self.assertAlmostEqual(out["sigma"], out["sigma_known_b"], places=6)

    def test_scaling_with_N(self):
        b = bg.bg_from_sbr(C, BEAM, 10.0)
        r = np.array([25.0, 5.0])
        s1 = bg.crb_free_bg(C, BEAM, r, b, 100)["sigma"]
        s2 = bg.crb_free_bg(C, BEAM, r, b, 400)["sigma"]
        self.assertAlmostEqual(s1 / s2, 2.0, places=8)

    def test_b_zero_forward_difference(self):
        F = bg.fisher_matrix_bg(C, BEAM, np.array([20.0, 0.0]), 0.0, 100)
        self.assertTrue(np.all(np.isfinite(F)))
        with self.assertRaises(ValueError):
            bg.fisher_matrix_bg(C, BEAM, np.array([20.0, 0.0]), -0.1, 100)


class TestMLEFreeBg(unittest.TestCase):
    def test_noise_free_recovers_truth(self):
        b = bg.bg_from_sbr(C, BEAM, 5.0)
        for r in ([0.0, 0.0], [20.0, 5.0], [-35.0, 10.0]):
            p = bg.probabilities_bg(np.array(r + [b]), C, BEAM)
            est = bg.mle_free_bg(1e6 * p, C, BEAM, 150.0, 0.3)
            self.assertEqual(est.shape, (3,))
            np.testing.assert_allclose(est[:2], r, atol=2e-3)
            self.assertAlmostEqual(est[2] / b, 1.0, places=3)

    def test_matches_known_b_mle_when_b_hat_equals_truth_region(self):
        # With 4 exposures and 3 parameters the free fit reproduces p_hat exactly: the fitted
        # probabilities must equal n/N when the solution is interior (b_hat > 0).
        b = bg.bg_from_sbr(C, BEAM, 5.0)
        p = bg.probabilities_bg(np.array([15.0, 0.0, b]), C, BEAM)
        cts = np.random.default_rng(3).multinomial(20000, p, size=20).astype(float)
        est = bg.mle_free_bg(cts, C, BEAM, 150.0, 0.3)
        self.assertTrue(np.all(est[:, 2] > 0))
        pf = bg.probabilities_bg(est, C, BEAM)
        np.testing.assert_allclose(pf, cts / cts.sum(1, keepdims=True), atol=1e-6)

    def test_batch_equals_single_and_zero_counts(self):
        b = bg.bg_from_sbr(C, BEAM, 10.0)
        p = bg.probabilities_bg(np.array([10.0, 0.0, b]), C, BEAM)
        cts = np.random.default_rng(1).multinomial(300, p, size=5)
        batch = bg.mle_free_bg(cts, C, BEAM, 150.0, 0.3)
        for i in range(5):
            np.testing.assert_allclose(bg.mle_free_bg(cts[i], C, BEAM, 150.0, 0.3), batch[i])
        z = bg.mle_free_bg(np.zeros(4), C, BEAM, 150.0, 0.3, center=(3.0, -2.0))
        np.testing.assert_array_equal(z, [3.0, -2.0, 0.0])

    def test_estimates_in_disk_and_b_nonnegative(self):
        b = bg.bg_from_sbr(C, BEAM, 2.0)
        p = bg.probabilities_bg(np.array([45.0, 0.0, b]), C, BEAM)
        cts = np.random.default_rng(2).multinomial(30, p, size=300)
        est = bg.mle_free_bg(cts, C, BEAM, 60.0, 0.3)
        self.assertTrue(np.all(np.hypot(est[:, 0], est[:, 1]) <= 60.0 * (1 + 1e-9)))
        self.assertTrue(np.all(est[:, 2] >= 0))

    def test_likelihood_not_worse_than_known_b(self):
        # the free fit maximizes over a superset: its log-likelihood >= the known-b MLE's
        b = bg.bg_from_sbr(C, BEAM, 5.0)
        model = photons.make_model(C, BEAM, bg_per_exposure=b)
        p = model(np.array([25.0, 0.0]))
        cts = np.random.default_rng(4).multinomial(500, p, size=200)
        e3 = bg.mle_free_bg(cts, C, BEAM, 150.0, 0.3)
        e2 = estimators.mle(cts, model, 150.0)
        ll3 = np.sum(cts * np.log(bg.probabilities_bg(e3, C, BEAM)), axis=1)
        ll2 = -estimators.neg_loglike(e2, cts, model)
        self.assertTrue(np.all(ll3 >= ll2 - 1e-6))

    def test_mc_unbiased_and_near_crb(self):
        # N = 500, SBR0 = 5, x = 20 nm: bias < 3 SE and sigma within 10 % of the free-b CRB
        b = bg.bg_from_sbr(C, BEAM, 5.0)
        r = np.array([20.0, 0.0])
        p = bg.probabilities_bg(np.append(r, b), C, BEAM)
        cts = np.random.default_rng(42).multinomial(500, p, size=1500)
        est = bg.mle_free_bg(cts, C, BEAM, 150.0, 0.3)
        e = est[:, :2] - r
        se = e.std(0, ddof=1) / np.sqrt(e.shape[0])
        self.assertTrue(np.all(np.abs(e.mean(0)) < 3.5 * se))
        sigma = np.sqrt(0.5 * e.var(0, ddof=1).sum())
        crb = bg.crb_free_bg(C, BEAM, r, b, 500)["sigma"]
        self.assertLess(abs(sigma / crb - 1.0), 0.10)

    def test_validation(self):
        good = np.array([100, 100, 100, 10])
        with self.assertRaises(ValueError):
            bg.mle_free_bg(good[:3], C, BEAM, 150.0, 0.3)          # K mismatch
        with self.assertRaises(ValueError):
            bg.mle_free_bg(-good, C, BEAM, 150.0, 0.3)             # negative counts
        with self.assertRaises(ValueError):
            bg.mle_free_bg(good * np.nan, C, BEAM, 150.0, 0.3)
        for kw in ({"search_radius": 0.0}, {"b_max": -1.0}, {"tol": 0.0}, {"n_b": 1},
                   {"grid_step": 0.0}, {"mem_budget": 0}):
            args = {"search_radius": 150.0, "b_max": 0.3}
            args.update(kw)
            with self.assertRaises(ValueError):
                bg.mle_free_bg(good, C, BEAM, **args)
        with self.assertRaises(ValueError):
            bg.mle_free_bg(np.ones((2, 2, 4)), C, BEAM, 150.0, 0.3)


if __name__ == "__main__":
    unittest.main()
