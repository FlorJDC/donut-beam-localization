# -*- coding: utf-8 -*-
"""Tests for donutloc.experiments (iterative MINFLUX, eps x L, misalignment drivers).

Small configurations only (n_rep <= 200); reference-size numbers go in the round report.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import experiments as ex


class TestIterative(unittest.TestCase):

    def test_beats_camera_no_background_and_sbr10(self):
        # N_total = 1000, L0 = 150 -> L_min = 25 nm, 4 equal photon splits
        for sbr in (None, 10):
            out = ex.iterative_minflux(1000, L0=150.0, L_min=25.0, n_iter=4, sbr=sbr, n_rep=200)
            self.assertAlmostEqual(out["camera_sigma"], 100.0 / np.sqrt(1000.0), places=12)
            self.assertLess(out["sigma"], out["camera_sigma"])
            self.assertLess(out["rmse"], out["camera_sigma"])
            np.testing.assert_allclose(out["L"], 150.0 * (25.0 / 150.0) ** (np.arange(4) / 3.0))
            self.assertEqual(int(out["N_k"].sum()), 1000)
            # final spread is close to the centre CRB of the last iteration (250 photons, L=25)
            self.assertLess(abs(out["sigma"] / out["crb_center_iter"][-1] - 1.0), 0.15)
            # zooming in helps: sigma decreases along the iterations
            self.assertTrue(np.all(np.diff(out["sigma_iter"]) < 0))

    def test_sigma_decreases_with_N(self):
        res = ex.iterative_vs_photons([250, 1000, 4000], n_rep=150)
        self.assertTrue(np.all(np.diff(res["sigma"]) < 0))
        self.assertTrue(-0.7 < res["slope"] < -0.3, res["slope"])
        np.testing.assert_allclose(res["camera_sigma"], 100.0 / np.sqrt([250, 1000, 4000]))
        self.assertNotIn("estimates", res["runs"][0])

    def test_reproducible_seed(self):
        a = ex.iterative_minflux(400, n_rep=50, seed=42)
        b = ex.iterative_minflux(400, n_rep=50, seed=42)
        c = ex.iterative_minflux(400, n_rep=50, seed=7)
        np.testing.assert_array_equal(a["estimates"], b["estimates"])
        self.assertFalse(np.array_equal(a["estimates"], c["estimates"]))
        self.assertEqual(a["params"]["seed"], 42)

    def test_no_recenter_is_worse(self):
        a = ex.iterative_minflux(1000, n_rep=100, recenter=True)
        b = ex.iterative_minflux(1000, n_rep=100, recenter=False)
        self.assertGreater(b["rmse"], 3.0 * a["rmse"])

    def test_adaptive_schedule(self):
        Nk = np.array([250, 250, 250, 250])
        L = ex.l_schedule(4, 150.0, 25.0, rule="adaptive", N_k=Nk)
        self.assertEqual(L[0], 150.0)
        self.assertTrue(np.all(np.diff(L) <= 0) and np.all(L >= 25.0))
        out = ex.iterative_minflux(1000, rule="adaptive", n_rep=100)
        self.assertLess(out["sigma"], out["camera_sigma"])

    def test_fixed_background_lowers_sbr(self):
        out = ex.iterative_minflux(1000, bg_per_exposure=0.03, n_rep=60)
        self.assertTrue(np.all(np.diff(out["sbr_center"]) < 0))   # Eq. S32: SBR drops with L

    def test_explicit_schedule_and_split(self):
        out = ex.iterative_minflux(1001, L_schedule=[120.0, 40.0], photon_split=[1, 3],
                                   n_rep=30)
        np.testing.assert_array_equal(out["L"], [120.0, 40.0])
        np.testing.assert_array_equal(out["N_k"], [250, 751])
        self.assertEqual(out["params"]["rule"], "explicit")

    def test_errors(self):
        with self.assertRaises(ValueError):
            ex.iterative_minflux(3, n_iter=4, n_rep=10)
        with self.assertRaises(ValueError):
            ex.iterative_minflux(1000, rule="bogus", n_rep=10)
        with self.assertRaises(ValueError):
            ex.iterative_minflux(1000, photon_split="log", n_rep=10)
        with self.assertRaises(ValueError):
            ex.iterative_minflux(1000, photon_split=[1, 2], n_rep=10)
        with self.assertRaises(ValueError):
            ex.iterative_minflux(1000, n_rep=1)
        with self.assertRaises(ValueError):
            ex.l_schedule(3, rule="adaptive")


class TestEpsL(unittest.TestCase):

    def test_eps0_linear_in_L(self):
        L = np.array([10.0, 20.0, 40.0])
        s = ex.eps_L_sweep([0.0], L, N=100)
        c = s["crb_center"][0]
        self.assertTrue(np.all(np.diff(c) > 0))
        self.assertTrue(all(m == "limit" for m in s["method"][0]))
        # quadratic-zero limit: sqrt(0.1) L / sqrt(N) (limit value, note A section 4.3)
        np.testing.assert_allclose(c, np.sqrt(0.1) * L / 10.0, rtol=0.01)
        # L = 50, fwhm = 300, N = 100 limit value used elsewhere in the project
        self.assertAlmostEqual(ex.crb_center(50.0), 1.6051, places=3)

    def test_eps_interior_minimum(self):
        L = np.array([5.0, 10.0, 25.0, 50.0, 100.0, 200.0])
        s = ex.eps_L_sweep([0.05], L, N=100)
        c = s["crb_center"][0]
        j = int(np.argmin(c))
        self.assertTrue(0 < j < len(L) - 1)
        self.assertEqual(s["method"][0, 0], "point")

    def test_optimal_L_finite_and_growing(self):
        Ls = []
        for e in (0.002, 0.01, 0.05, 0.15):
            o = ex.optimal_L(e, N=100)
            self.assertTrue(o["interior"])
            Ls.append(o["L_opt"])
        self.assertTrue(np.all(np.diff(Ls) > 0), Ls)
        o0 = ex.optimal_L(0.0, N=100)
        self.assertFalse(o0["interior"])
        self.assertEqual(o0["L_opt"], 1.0)
        o10 = ex.optimal_L(0.05, N=100, sbr=10)
        self.assertTrue(o10["interior"])

    def test_point_equals_limit_for_eps(self):
        # with eps > 0 the CRB is continuous at the centre
        from donutloc import beams, fisher, patterns, photons
        p = photons.make_model(patterns.tcp_centers(50.0), beams.make_beam(eps=0.05))
        self.assertAlmostEqual(ex.crb_center(50.0, eps=0.05), fisher.crb_limit(p, 100), places=6)

    def test_fov_mean(self):
        s = ex.eps_L_sweep([0.0, 0.05], [50.0], fov=True)
        self.assertTrue(np.all(s["crb_fov_mean"] > s["crb_center"]))


class TestMisalignment(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.m = ex.misalignment_study([0.0, 10.0], L=100.0, N=500, sbr=10, n_patterns=6,
                                      n_rep=100)

    def test_delta0_identical(self):
        m = self.m
        self.assertTrue(m["identical"][0])
        for k in ("bias_abs_mean", "sigma_mean", "rmse"):
            np.testing.assert_array_equal(m["honest"][k][0], m["naive"][k][0])

    def test_naive_biased_honest_not(self):
        m = self.m
        h, n = m["honest"], m["naive"]
        self.assertFalse(m["identical"][1])
        # naive bias >> honest bias and outside 3 SE; honest compatible with MC noise
        self.assertTrue(np.all(n["bias_abs_mean"][1] > h["bias_abs_mean"][1]
                               + 3 * n["bias_abs_se"][1]))
        self.assertTrue(np.all(h["bias_abs_mean"][1] < m["bias_noise_floor"][1]
                               + 3 * h["bias_abs_se"][1]))
        self.assertTrue(np.all(h["bias_chi2"] < 1.0 + 3.0 / np.sqrt(6) * 1.5))
        self.assertTrue(np.all(n["bias_chi2"][1] > 20))
        self.assertTrue(np.all(n["rmse"][1] > 1.5 * h["rmse"][1]))
        # honest MLE is efficient (SBR = 10)
        np.testing.assert_allclose(h["sigma_mean"], m["crb_honest"], rtol=0.12)

    def test_reproducible_and_lms(self):
        a = ex.misalignment_study([5.0], n_patterns=2, n_rep=40, positions=[(0.0, 0.0)])
        b = ex.misalignment_study([5.0], n_patterns=2, n_rep=40, positions=[(0.0, 0.0)])
        np.testing.assert_array_equal(a["naive"]["rmse"], b["naive"]["rmse"])
        c = ex.misalignment_study([0.0], n_patterns=2, n_rep=40, positions=[(0.0, 0.0)],
                                  estimator="lms")
        self.assertTrue(c["identical"][0])
        with self.assertRaises(ValueError):
            ex.misalignment_study([1.0], n_patterns=1, n_rep=10, estimator="bogus")
        with self.assertRaises(ValueError):
            ex.misalignment_study([-1.0], n_patterns=1, n_rep=10)


class TestErrorStats(unittest.TestCase):

    def test_values(self):
        est = np.array([[1.0, 0.0], [-1.0, 0.0], [1.0, 2.0], [-1.0, 2.0]])
        st = ex.error_stats(est, np.array([0.0, 1.0]))
        np.testing.assert_allclose(st["bias"], [0.0, 0.0])
        var = 4.0 / 3.0
        self.assertAlmostEqual(st["sigma"], np.sqrt(var))
        self.assertAlmostEqual(st["sigma_se_gauss"], np.sqrt(var) / 4.0)
        self.assertTrue(np.isfinite(st["sigma_se"]) and st["sigma_se"] > 0)   # bootstrap
        self.assertAlmostEqual(ex.error_stats(est, np.array([0.0, 1.0]), n_boot=0)["sigma_se"],
                               np.sqrt(var) / 4.0)
        self.assertAlmostEqual(st["rmse"], 1.0)
        with self.assertRaises(ValueError):
            ex.error_stats(est[:1], [0.0, 0.0])


class TestR3Fixes(unittest.TestCase):

    def test_error_stats_bootstrap_matches_gaussian_formula(self):
        rng = np.random.default_rng(11)
        e = rng.normal(0.0, 2.0, size=(3000, 2))
        st = ex.error_stats(e, [0.0, 0.0])
        self.assertAlmostEqual(st["sigma_se"] / st["sigma_se_gauss"], 1.0, delta=0.10)

    def test_iterative_returns_bootstrap_se(self):
        out = ex.iterative_minflux(1000, n_rep=300)
        self.assertIn("sigma_se_gauss", out)
        self.assertEqual(out["sigma_se_iter"].shape, (4,))
        self.assertAlmostEqual(out["sigma_se"], out["sigma_se_iter"][-1])
        self.assertTrue(0.5 < out["sigma_se"] / out["sigma_se_gauss"] < 3.0)

    def test_misalignment_sigma_se_between_patterns(self):
        """R2 refuted scenario: delta = 10 at the centre, 20 patterns x 200 reps.  The old formula
        (now sigma_se_within) ignored the pattern-to-pattern variance (2.6x / 3.6x too small)."""
        m = ex.misalignment_study([10.0], L=100.0, N=500, sbr=10, n_patterns=20, n_rep=200,
                                  positions=[(0.0, 0.0)])
        for est in ("honest", "naive"):
            se, within = m[est]["sigma_se"][0, 0], m[est]["sigma_se_within"][0, 0]
            self.assertAlmostEqual(within, m[est]["sigma_mean"][0, 0] / (2 * np.sqrt(200 * 20)))
            self.assertGreater(se, 1.5 * within, est)

    def test_misalignment_sigma_se_formula_exact(self):
        """sigma_se = std over patterns (ddof 1) of per-pattern sigmas / sqrt(P), recomputed
        independently with the documented seeding (default_rng(seed) for the patterns,
        SeedSequence(seed) for the MC seeds)."""
        from donutloc import beams, estimators, montecarlo, patterns, photons
        P, R, L, N, d = 3, 60, 100.0, 500, 5.0
        m = ex.misalignment_study([d], L=L, N=N, sbr=10, n_patterns=P, n_rep=R,
                                  positions=[(0.0, 0.0)], seed=42)
        beam = beams.make_beam("donut", fwhm=300.0)
        ideal = patterns.tcp_centers(L)
        p_naive = photons.make_model(ideal, beam, sbr=10)
        rng = np.random.default_rng(42)
        seeds = np.random.SeedSequence(42).generate_state(P).reshape(P, 1)
        sig = []
        for j in range(P):
            p_true = photons.make_model(patterns.perturb_centers(ideal, d, rng=rng), beam, sbr=10)
            o = montecarlo.run_mc(lambda C: estimators.mle(C, p_naive, search_radius=0.75 * L),
                                  p_true, [0.0, 0.0], N, R, seed=int(seeds[j, 0]), n_boot=0)
            sig.append(o["sigma"])
        sig = np.array(sig)
        self.assertAlmostEqual(m["naive"]["sigma_mean"][0, 0], sig.mean(), places=10)
        self.assertAlmostEqual(m["naive"]["sigma_se"][0, 0], sig.std(ddof=1) / np.sqrt(P),
                               places=10)
        one = ex.misalignment_study([d], n_patterns=1, n_rep=20, positions=[(0.0, 0.0)])
        self.assertTrue(np.isnan(one["naive"]["sigma_se"][0, 0]))

    def test_adaptive_schedule_coverage_is_about_97_percent(self):
        """R2 refuted docstring: with the defaults the next pattern radius L_1/2 contains ~97 %
        of the emitters (not 99 %), because sigma_k is the centre CRB and underestimates the
        real error of iteration 0 (emitters up to L0/4 off-centre)."""
        Nk = np.array([250, 250, 250, 250])
        L = ex.l_schedule(4, 150.0, rule="adaptive", N_k=Nk)
        self.assertAlmostEqual(L[1], 25.0)
        it0 = ex.iterative_minflux(250, L_schedule=[150.0], n_rep=3000)
        err = np.hypot(*(it0["estimates"] - it0["r_true"]).T)
        outside = float(np.mean(err > L[1] / 2.0))
        self.assertTrue(0.015 < outside < 0.045, outside)
        self.assertGreater(it0["sigma"], 1.1 * it0["crb_center_iter"][0])
        self.assertIn("97 %", ex.l_schedule.__doc__)
        self.assertNotIn("99 % probability", ex.l_schedule.__doc__)


class TestMisalignmentPopulation(unittest.TestCase):
    """Noise-free (population) naive-MLE bias under misalignment (inbox r4)."""

    def test_zero_displacement_gives_zero_bias(self):
        m = ex.misalignment_population_bias([0.0], n_patterns=20)
        self.assertLess(float(np.max(m["bias_abs_mean"])), 1e-3)
        self.assertTrue(np.all(np.isnan(m["ratio"])))
        self.assertTrue(np.all(np.isnan(m["slope"])))

    def test_converges_as_delta_to_zero(self):
        """|bias|/delta tends to a finite constant as delta -> 0 (the bias is linear in delta
        for small delta): the ratios at delta = 0.25 and 0.5 agree to < 2 % and lie close to
        the delta = 2 value; the bias itself vanishes."""
        m = ex.misalignment_population_bias([0.25, 0.5, 2.0], n_patterns=300)
        r = m["ratio"]
        for q in range(2):
            self.assertLess(abs(r[0, q] / r[1, q] - 1.0), 0.02, r[:, q])
            self.assertLess(abs(r[1, q] / r[2, q] - 1.0), 0.05, r[:, q])
        self.assertTrue(np.all(m["bias_abs_mean"][0] < 0.25))
        # the ratio at the centre is ~0.75 and larger at (L/4, 0)
        self.assertTrue(0.70 < r[0, 0] < 0.80, r[0, 0])
        self.assertGreater(r[0, 1], r[0, 0])

    def test_slope_definition_and_same_patterns_as_study(self):
        d = [2.0, 5.0, 10.0]
        m = ex.misalignment_population_bias(d, n_patterns=30, positions=[(0.0, 0.0)])
        b = m["bias_abs_mean"][:, 0]
        dd = np.asarray(d)
        self.assertAlmostEqual(m["slope"][0], float(np.sum(dd * b) / np.sum(dd ** 2)), places=12)
        np.testing.assert_allclose(m["ratio"][:, 0], b / dd)
        self.assertTrue(np.isfinite(m["slope_se"][0]) and m["slope_se"][0] > 0)
        # pattern k is the same draw as misalignment_study's pattern k (rng re-created per delta)
        from donutloc import beams, estimators, patterns, photons
        L = 100.0
        rng = np.random.default_rng(42)
        ideal = patterns.tcp_centers(L)
        beam = beams.make_beam("donut", fwhm=300.0)
        true_c = patterns.perturb_centers(ideal, 5.0, rng=rng)
        cnt = 500.0 * photons.make_model(true_c, beam, sbr=10)(np.zeros(2))
        est = estimators.mle(cnt, photons.make_model(ideal, beam, sbr=10), search_radius=0.75 * L)
        one = ex.misalignment_population_bias([5.0], n_patterns=2, positions=[(0.0, 0.0)])
        b0 = float(np.hypot(*est))
        self.assertGreater(one["bias_abs_se"][0, 0], 0.0)
        self.assertAlmostEqual(2 * one["bias_abs_mean"][0, 0] - b0,
                               float(np.hypot(*_second_pattern_bias())), places=8)

    def test_agrees_with_monte_carlo_centre_within_3_se(self):
        """Paired check against the 400 x 200 Monte Carlo stored in data/paper_numbers.json:
        with the SAME 400 patterns (same seed and draw order as misalignment_study), the
        noise-free naive |bias| at the centre agrees with the MC |bias| within 3 SE for
        delta = 2, 5, 10, and the MC-weighted slope of the population values agrees with the MC
        slope within 3 SE."""
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "data", "paper_numbers.json")
        if not os.path.exists(path):
            self.skipTest("data/paper_numbers.json missing")
        import json
        with open(path, encoding="utf-8") as fh:
            pn = json.load(fh)
        d = np.array([2.0, 5.0, 10.0])
        m = ex.misalignment_population_bias(d, n_patterns=400, positions=[(0.0, 0.0)])
        mc_b = np.array([pn["misalignment_naive_bias_abs_center_d%d_nm" % k]["value"]
                         for k in (2, 5, 10)])
        mc_se = np.array([pn["misalignment_naive_bias_abs_center_d%d_nm" % k]["se"]
                          for k in (2, 5, 10)])
        z = (m["bias_abs_mean"][:, 0] - mc_b) / mc_se
        self.assertTrue(np.all(np.abs(z) < 3.0), z)
        w = 1.0 / mc_se ** 2
        s_pop = np.sum(w * d * m["bias_abs_mean"][:, 0]) / np.sum(w * d ** 2)
        mc = pn["misalignment_naive_bias_over_delta_center"]
        self.assertLess(abs(s_pop - mc["value"]) / mc["se"], 3.0, (s_pop, mc["value"]))

    def test_validation(self):
        with self.assertRaises(ValueError):
            ex.misalignment_population_bias([-1.0], n_patterns=5)
        with self.assertRaises(ValueError):
            ex.misalignment_population_bias([1.0], n_patterns=1)


def _second_pattern_bias():
    """Naive-MLE error of the 2nd random pattern (delta = 5, centre), drawn independently."""
    from donutloc import beams, estimators, patterns, photons
    L = 100.0
    rng = np.random.default_rng(42)
    ideal = patterns.tcp_centers(L)
    beam = beams.make_beam("donut", fwhm=300.0)
    patterns.perturb_centers(ideal, 5.0, rng=rng)
    true_c = patterns.perturb_centers(ideal, 5.0, rng=rng)
    cnt = 500.0 * photons.make_model(true_c, beam, sbr=10)(np.zeros(2))
    return estimators.mle(cnt, photons.make_model(ideal, beam, sbr=10), search_radius=0.75 * L)


if __name__ == "__main__":
    unittest.main()
