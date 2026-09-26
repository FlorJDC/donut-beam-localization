import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import beams, patterns, photons


def _lg_ref(x, y, fwhm):
    r2 = x * x + y * y
    a = 4.0 * np.log(2.0) / fwhm ** 2
    return np.e * a * r2 * np.exp(-a * r2)


class TestPhotons(unittest.TestCase):

    def setUp(self):
        self.L, self.fwhm = 100.0, 300.0
        self.c = patterns.tcp_centers(self.L)
        self.beam = beams.make_beam("donut", fwhm=self.fwhm)

    def test_intensities_shape_and_values(self):
        r = np.array([[3.0, -7.0], [20.0, 10.0], [0.0, 0.0]])
        lam = photons.intensities(r, self.c, self.beam)
        self.assertEqual(lam.shape, (3, 4))
        for m in range(3):
            for k in range(4):
                ref = _lg_ref(r[m, 0] - self.c[k, 0], r[m, 1] - self.c[k, 1], self.fwhm)
                self.assertAlmostEqual(lam[m, k], ref, places=14)
        # arbitrary leading shape
        grid = np.zeros((5, 6, 2))
        self.assertEqual(photons.intensities(grid, self.c, self.beam).shape, (5, 6, 4))
        self.assertEqual(photons.intensities([1.0, 2.0], self.c, self.beam).shape, (4,))

    def test_probabilities_sum_to_one(self):
        rng = np.random.default_rng(42)
        r = rng.uniform(-60, 60, size=(50, 2))
        for kw in ({}, {"sbr": 5.0}, {"sbr": np.inf}, {"bg_per_exposure": 0.01}):
            p = photons.probabilities(r, self.c, self.beam, **kw)
            self.assertEqual(p.shape, (50, 4))
            np.testing.assert_allclose(p.sum(axis=-1), 1.0, rtol=1e-13)
            self.assertTrue(np.all(p >= 0))

    def test_center_prob_zero_without_background(self):
        p = photons.probabilities([0.0, 0.0], self.c, self.beam)
        self.assertEqual(p[3], 0.0)
        np.testing.assert_allclose(p[:3], 1.0 / 3.0, rtol=1e-14)
        np.testing.assert_array_equal(photons.probabilities([0.0, 0.0], self.c, self.beam,
                                                            sbr=np.inf), p)

    def test_sbr_eq_s30_at_center(self):
        for sbr in (1.0, 5.0, 10.0, 13.6):
            p = photons.probabilities([0.0, 0.0], self.c, self.beam, sbr=sbr)
            # Balzarotti2017 Eq. S30 with p^(0)(0) = (1/3, 1/3, 1/3, 0), K = 4
            np.testing.assert_allclose(p[:3], sbr / (sbr + 1) / 3 + 1 / (sbr + 1) / 4, rtol=1e-14)
            self.assertAlmostEqual(p[3], 1 / (sbr + 1) / 4, places=15)
        np.testing.assert_allclose(photons.probabilities([5.0, 1.0], self.c, self.beam, sbr=0),
                                   0.25, rtol=1e-15)

    def test_bg_per_exposure_eq_s28_s29(self):
        r = np.array([[4.0, -3.0], [0.0, 0.0]])
        b = 0.02
        lam = photons.intensities(r, self.c, self.beam)
        p = photons.probabilities(r, self.c, self.beam, bg_per_exposure=b)
        np.testing.assert_allclose(p, (lam + b) / (lam + b).sum(-1, keepdims=True), rtol=1e-14)
        sbr = photons.sbr_at(r, self.c, self.beam, b)
        np.testing.assert_allclose(sbr, lam.sum(-1) / (4 * b), rtol=1e-14)
        # fixed bg <=> Eq. S30 with the local SBR from Eq. S29
        for m in range(2):
            p30 = photons.probabilities(r[m], self.c, self.beam, sbr=sbr[m])
            np.testing.assert_allclose(p[m], p30, rtol=1e-12)
        self.assertTrue(np.isinf(photons.sbr_at([1.0, 1.0], self.c, self.beam, 0.0)))

    def test_constant_pedestal_equals_effective_background(self):
        # I_LG + eps (constant) is exactly a background eps per exposure.
        rng = np.random.default_rng(42)
        r = rng.uniform(-40, 40, size=(20, 2))
        r[0] = 0.0
        K = 4
        for eps in (0.01, 0.1):
            ped = beams.make_beam("donut", fwhm=self.fwhm, eps=eps, zero_model="constant")
            p_ped = photons.probabilities(r, self.c, ped)
            p_bg = photons.probabilities(r, self.c, self.beam, bg_per_exposure=eps)
            np.testing.assert_allclose(p_ped, p_bg, rtol=1e-13, atol=1e-16)
            S = photons.intensities(r, self.c, self.beam).sum(-1)   # sum of pure-LG intensities
            # (i) no other background: 1/SBR_eff = K eps / sum I_LG
            sbr_eff = S / (K * eps)
            for m in range(len(r)):
                np.testing.assert_allclose(
                    p_ped[m], photons.probabilities(r[m], self.c, self.beam, sbr=sbr_eff[m]),
                    rtol=1e-12)
            # (ii) with a fixed SBR defined on the pure-LG signal (bg total = S/SBR):
            #      1/SBR_eff = 1/SBR + K eps/S   (exact)
            for sbr in (5.0, 10.0):
                lam_ped = photons.intensities(r, self.c, ped)
                p_true = (lam_ped + (S / sbr / K)[:, None])
                p_true /= p_true.sum(-1, keepdims=True)
                sbr_eff = 1.0 / (1.0 / sbr + K * eps / S)
                for m in range(len(r)):
                    np.testing.assert_allclose(
                        p_true[m], photons.probabilities(r[m], self.c, self.beam, sbr=sbr_eff[m]),
                        rtol=1e-12)
                # (iii) photons.probabilities(sbr=...) on the pedestal beam defines the SBR on the
                #       full beam signal (S + K eps); then there is a cross term:
                #       1/SBR_eff = 1/SBR + (K eps/S)(1 + 1/SBR)
                p_impl = photons.probabilities(r, self.c, ped, sbr=sbr)
                sbr_eff = 1.0 / (1.0 / sbr + (K * eps / S) * (1.0 + 1.0 / sbr))
                for m in range(len(r)):
                    np.testing.assert_allclose(
                        p_impl[m], photons.probabilities(r[m], self.c, self.beam, sbr=sbr_eff[m]),
                        rtol=1e-12)

    def test_make_model_contract(self):
        p_fn = photons.make_model(self.c, self.beam, sbr=10.0)
        r = np.zeros((3, 2, 2))
        r[..., 0] = 7.0
        out = p_fn(r)
        self.assertEqual(out.shape, (3, 2, 4))
        np.testing.assert_allclose(out.sum(-1), 1.0, rtol=1e-14)
        np.testing.assert_array_equal(out[0, 0], photons.probabilities([7.0, 0.0], self.c,
                                                                       self.beam, sbr=10.0))
        self.assertEqual(p_fn.K, 4)
        self.assertEqual(p_fn([0.0, 0.0]).shape, (4,))
        # centers are copied: later mutation of the caller's array does not change the model
        c = self.c.copy()
        q_fn = photons.make_model(c, self.beam)
        before = q_fn([1.0, 2.0])
        c[:] = 0.0
        np.testing.assert_array_equal(q_fn([1.0, 2.0]), before)

    def test_errors(self):
        with self.assertRaises(ValueError):
            photons.probabilities([0, 0], self.c, self.beam, sbr=5, bg_per_exposure=0.1)
        with self.assertRaises(ValueError):
            photons.make_model(self.c, self.beam, sbr=5, bg_per_exposure=0.1)
        with self.assertRaises(ValueError):
            photons.probabilities([0, 0], self.c, self.beam, sbr=-1)
        with self.assertRaises(ValueError):
            photons.probabilities([0, 0], self.c, self.beam, bg_per_exposure=-1)
        with self.assertRaises(ValueError):
            photons.intensities([0, 0, 0], self.c, self.beam)
        with self.assertRaises(ValueError):
            photons.intensities([0, 0], np.ones(4), self.beam)

    def test_sample_counts_multinomial_mean(self):
        p = photons.probabilities([10.0, -5.0], self.c, self.beam, sbr=10.0)
        N, M = 100, 20000
        n = photons.sample_counts(p, N, size=M, rng=42)
        self.assertEqual(n.shape, (M, 4))
        self.assertTrue(np.issubdtype(n.dtype, np.integer))
        self.assertTrue(np.all(n.sum(axis=1) == N))
        se = np.sqrt(N * p * (1 - p) / M)
        self.assertTrue(np.all(np.abs(n.mean(0) - N * p) < 5 * se), (n.mean(0), N * p))
        # reproducible
        np.testing.assert_array_equal(n, photons.sample_counts(p, N, size=M,
                                                               rng=np.random.default_rng(42)))
        self.assertEqual(photons.sample_counts(p, N, rng=42).shape, (4,))

    def test_sample_counts_poisson_mean(self):
        p = np.array([0.1, 0.2, 0.3, 0.4])
        N, M = 150.0, 20000
        n = photons.sample_counts(p, N, size=M, rng=42, mode="poisson")
        self.assertEqual(n.shape, (M, 4))
        se = np.sqrt(N * p / M)
        self.assertTrue(np.all(np.abs(n.mean(0) - N * p) < 5 * se))
        self.assertTrue(np.all(np.abs(n.var(0) / (N * p) - 1) < 0.05))
        self.assertGreater(n.sum(1).std(), 1.0)       # total is not fixed

    def test_sample_counts_zero_prob_and_errors(self):
        p = photons.probabilities([0.0, 0.0], self.c, self.beam)
        n = photons.sample_counts(p, 50, size=100, rng=42)
        self.assertTrue(np.all(n[:, 3] == 0))
        with self.assertRaises(ValueError):
            photons.sample_counts([0.5, 0.6], 10)
        with self.assertRaises(ValueError):
            photons.sample_counts([1.5, -0.5], 10)
        with self.assertRaises(ValueError):
            photons.sample_counts([0.5, 0.5], 10.5)
        with self.assertRaises(ValueError):
            photons.sample_counts([0.5, 0.5], 10, mode="bogus")
        with self.assertRaises(ValueError):
            photons.sample_counts(np.full((2, 2), 0.25), 10)


class TestFarField(unittest.TestCase):

    def test_all_intensities_underflow(self):
        c = patterns.tcp_centers(50.0)
        g = beams.make_beam("gaussian", fwhm=300.0)
        r = np.array([[6000.0, 0.0], [10.0, 0.0]])
        self.assertTrue(np.all(photons.intensities(r[0], c, g) == 0.0))
        p = photons.probabilities(r, c, g)
        self.assertTrue(np.all(np.isnan(p[0])))
        np.testing.assert_allclose(p[1].sum(), 1.0)
        self.assertTrue(np.all(np.isnan(photons.probabilities(r[0], c, g, sbr=10.0))))
        np.testing.assert_allclose(photons.probabilities(r[0], c, g, bg_per_exposure=1e-3), 0.25)


if __name__ == "__main__":
    unittest.main()
