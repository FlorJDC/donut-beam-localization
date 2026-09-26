# -*- coding: utf-8 -*-
"""Tests for donutloc.montecarlo (and MLE efficiency)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import time
import unittest

import numpy as np

from donutloc import estimators as est
from donutloc import fisher
from donutloc import montecarlo as mc


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


def mle_for(p, R):
    return lambda C: est.mle(C, p, search_radius=R)


class TestRunMC(unittest.TestCase):

    def test_keys_shapes_reproducibility(self):
        p = tcp_model(50.0)
        f = lambda C: est.lms_tcp(C, 50.0, 300.0)
        a = mc.run_mc(f, p, [1.0, 2.0], 100, 300)
        b = mc.run_mc(f, p, [1.0, 2.0], 100, 300)
        for k in ("bias", "std", "sigma", "rmse", "sigma_err", "n_rep", "seed", "estimates"):
            self.assertIn(k, a)
        self.assertEqual(a["estimates"].shape, (300, 2))
        self.assertEqual(a["bias"].shape, (2,))
        self.assertEqual(a["seed"], 42)
        np.testing.assert_array_equal(a["estimates"], b["estimates"])
        c = mc.run_mc(f, p, [1.0, 2.0], 100, 300, seed=7)
        self.assertFalse(np.array_equal(a["estimates"], c["estimates"]))

    def test_statistics_identities(self):
        p = tcp_model(50.0)
        r = np.array([3.0, -2.0])
        res = mc.run_mc(lambda C: est.lms_tcp(C, 50.0, 300.0), p, r, 100, 1000)
        e, n = res["estimates"], res["n_rep"]
        self.assertAlmostEqual(res["sigma"], np.sqrt(0.5 * np.sum(res["std"] ** 2)), places=12)
        # rmse^2 = sigma^2 (n-1)/n + |bias|^2 / 2
        self.assertAlmostEqual(res["rmse"] ** 2,
                               res["sigma"] ** 2 * (n - 1) / n + 0.5 * np.sum(res["bias"] ** 2), places=10)
        self.assertAlmostEqual(res["sigma_err"], res["sigma"] / (2.0 * np.sqrt(n)), places=12)
        np.testing.assert_allclose(res["bias"], e.mean(0) - r)

    def test_sampling_modes(self):
        p = tcp_model(50.0)
        r = np.array([5.0, 5.0])
        res = mc.run_mc(lambda C: est.lms_tcp(C, 50.0, 300.0), p, r, 200, 4000)
        self.assertTrue(np.all(res["counts"].sum(1) == 200))
        np.testing.assert_allclose(res["counts"].mean(0), 200 * p(r), atol=0.5)
        res = mc.run_mc(lambda C: est.lms_tcp(C, 50.0, 300.0), p, r, 200, 4000, mode="poisson")
        self.assertGreater(res["counts"].sum(1).std(), 5.0)
        self.assertAlmostEqual(res["counts"].sum(1).mean() / 200.0, 1.0, delta=0.01)
        with self.assertRaises(ValueError):
            mc.run_mc(lambda C: est.lms_tcp(C, 50.0, 300.0), p, r, 200, 10, mode="bogus")

    def test_sigma_err_matches_seed_to_seed_scatter(self):
        # sigma_err = sigma / (2 sqrt(n)) (V13 correction of the R1 sigma/sqrt(2n))
        p = tcp_model(50.0)
        f = lambda C: est.lms_tcp(C, 50.0, 300.0)
        n = 400
        res = [mc.run_mc(f, p, [3.0, -2.0], 200, n, seed=s) for s in range(60)]
        sig = np.array([r["sigma"] for r in res])
        err = np.mean([r["sigma_err"] for r in res])
        self.assertAlmostEqual(sig.std(ddof=1) / err, 1.0, delta=0.3)
        # the old formula would be sqrt(2) too large
        self.assertLess(sig.std(ddof=1), 0.85 * np.mean(sig) / np.sqrt(2 * n))


class TestMLEEfficiency(unittest.TestCase):

    def test_mle_efficient_at_center_with_background(self):
        # regular model (p_centre > 0): std(MLE) / CRB in [0.9, 1.2]
        p = tcp_model(50.0, 300.0, sbr=10.0)
        res = mc.run_mc(mle_for(p, 50.0), p, [0.0, 0.0], 100, 500)
        eff = res["sigma"] / fisher.crb_limit(p, 100)
        self.assertTrue(0.9 <= eff <= 1.2, eff)

    def test_mle_efficient_off_center_without_background(self):
        p = tcp_model(50.0)
        r = [10.0, 0.0]
        res = mc.run_mc(mle_for(p, 50.0), p, r, 1000, 500)
        eff = res["sigma"] / float(fisher.crb(p, r, 1000))
        self.assertTrue(0.9 <= eff <= 1.15, eff)

    def test_mle_below_limit_crb_at_perfect_zero(self):
        # Documented finding (see reports/r01-worker-B.md): with a perfect zero and no background
        # the model is non-regular at the centre; the MLE is biased towards the centre (n_centre = 0
        # acts as a Gaussian pull) and its std falls BELOW the r->0 CRB, ~0.83 x crb_limit for
        # all N (measured 0.831-0.837 with 5000 reps for N = 100, 1e3, 1e4).
        p = tcp_model(50.0)
        res = mc.run_mc(mle_for(p, 50.0), p, [0.0, 0.0], 100, 500)
        eff = res["sigma"] / fisher.crb_limit(p, 100)
        self.assertTrue(0.75 <= eff <= 0.92, eff)
        # and the MLE is biased towards the centre slightly off-centre
        res = mc.run_mc(mle_for(p, 50.0), p, [2.0, 0.0], 100, 2000)
        self.assertLess(res["bias"][0], -4 * res["std"][0] / np.sqrt(2000))

    def test_speed_2000_reps(self):
        p = tcp_model(50.0)
        t = time.time()
        mc.run_mc(mle_for(p, 50.0), p, [0.0, 0.0], 100, 2000)
        self.assertLess(time.time() - t, 30.0)


if __name__ == "__main__":
    unittest.main()
