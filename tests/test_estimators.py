# -*- coding: utf-8 -*-
"""Tests for donutloc.estimators (MLE, LMS, mLMS)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import estimators as est


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


class TestNegLogLike(unittest.TestCase):

    def test_value_and_shape(self):
        p = tcp_model(50.0)
        n = np.array([10, 20, 30, 0])
        r = np.array([[1.0, 2.0], [3.0, -4.0], [0.0, 0.0]])
        v = est.neg_loglike(r, n, p)
        self.assertEqual(v.shape, (3,))
        ref = -np.sum(n * np.log(p(r[1])))
        self.assertAlmostEqual(v[1], ref, places=10)
        # p_centre = 0 at the origin with n_centre > 0: finite thanks to the 1e-300 clip
        self.assertTrue(np.isfinite(est.neg_loglike([0.0, 0.0], [1, 1, 1, 1], p)))


class TestMLE(unittest.TestCase):

    def test_noise_free_recovery(self):
        L = 50.0
        p = tcp_model(L)
        truths = np.array([[0.0, 0.0], [3.3, -1.7], [-10.0, 6.0], [0.0, 15.0], [12.0, -12.0]])
        counts = np.round(1e7 * p(truths))            # expected counts, large integers
        got = est.mle(counts, p, search_radius=L)
        self.assertEqual(got.shape, (5, 2))
        np.testing.assert_allclose(got, truths, atol=5e-3)

    def test_noise_free_recovery_with_background_and_rotation(self):
        p = tcp_model(100.0, fwhm=360.0, sbr=5.0, rotation=0.3)
        truth = np.array([8.0, -21.0])
        got = est.mle(np.round(1e7 * p(truth)), p, search_radius=100.0)
        self.assertEqual(got.shape, (2,))
        np.testing.assert_allclose(got, truth, atol=5e-3)

    def test_single_equals_batch_row(self):
        p = tcp_model(50.0)
        C = np.random.default_rng(42).multinomial(100, p(np.array([2.0, 1.0])), size=4)
        batch = est.mle(C, p, search_radius=50.0)
        for m in range(4):
            np.testing.assert_allclose(est.mle(C[m], p, search_radius=50.0), batch[m])

    def test_zero_counts_returns_center(self):
        p = tcp_model(50.0)
        got = est.mle(np.zeros(4), p, search_radius=50.0, center=(1.0, -2.0))
        np.testing.assert_allclose(got, [1.0, -2.0])

    def test_estimate_stays_in_disk(self):
        p = tcp_model(50.0)
        # all photons in one peripheral exposure pulls the estimate far away
        C = np.array([[0, 0, 50, 0], [5, 0, 0, 0], [0, 3, 0, 0]])
        for R in (10.0, 30.0):
            got = est.mle(C, p, search_radius=R)
            self.assertTrue(np.all(np.hypot(got[:, 0], got[:, 1]) <= R * (1 + 1e-9)))

    def test_small_N_and_zero_center_counts(self):
        p = tcp_model(50.0)
        C = np.array([[1, 0, 0, 0], [2, 1, 1, 0], [0, 0, 0, 1], [3, 3, 3, 0]])
        got = est.mle(C, p, search_radius=50.0)
        self.assertTrue(np.all(np.isfinite(got)))
        np.testing.assert_allclose(got[3], [0.0, 0.0], atol=2e-3)   # symmetric counts -> centre

    def test_refinement_is_a_local_optimum_and_matches_scipy(self):
        p = tcp_model(50.0)
        C = np.random.default_rng(42).multinomial(100, p(np.array([4.0, -3.0])), size=5)
        a = est.mle(C, p, search_radius=50.0)
        b = est.mle(C, p, search_radius=50.0, refine="scipy")
        g = est.mle(C, p, search_radius=50.0, refine=False)
        np.testing.assert_allclose(a, b, atol=2e-3)
        for m in range(5):
            f0 = est.neg_loglike(a[m], C[m], p)
            self.assertLessEqual(f0, est.neg_loglike(g[m], C[m], p) + 1e-12)
            nb = a[m] + 1e-2 * np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
            self.assertTrue(np.all(est.neg_loglike(nb, C[m], p) >= f0 - 1e-9))

    def test_bad_args(self):
        p = tcp_model(50.0)
        with self.assertRaises(ValueError):
            est.mle([1, 2, 3, 4], p, search_radius=0.0)
        with self.assertRaises(ValueError):
            est.mle([1, 2, 3, 4], p, search_radius=10.0, refine="bogus")

    def test_nonpositive_tol_step_radius_raise_instead_of_hanging(self):
        # R1 defect: tol=0 made the pattern search loop forever
        p = tcp_model(50.0)
        c = [10, 10, 10, 1]
        for kw in ({"tol": 0.0}, {"tol": -1e-3}, {"tol": np.nan}, {"grid_step": 0.0},
                   {"grid_step": -1.0}, {"search_radius": -5.0}, {"search_radius": np.inf},
                   {"chunk": 0}, {"mem_budget": 0}):
            args = {"search_radius": 50.0}
            args.update(kw)
            with self.assertRaises(ValueError, msg=str(kw)):
                est.mle(c, p, **args)

    def test_counts_validation(self):
        p = tcp_model(50.0)
        with self.assertRaisesRegex(ValueError, ">= 0"):
            est.mle([-5, 10, 10, 0], p, 50.0)                # R1: silently accepted
        with self.assertRaisesRegex(ValueError, "finite"):
            est.mle([np.nan, 10, 10, 0], p, 50.0)
        with self.assertRaisesRegex(ValueError, "finite"):
            est.mle([[1, 2, 3, 4], [np.inf, 1, 1, 1]], p, 50.0)
        with self.assertRaisesRegex(ValueError, "K = 3 exposures but p_fn returns K = 4"):
            est.mle([10, 10, 10], p, 50.0)                   # R1: cryptic matmul error
        with self.assertRaises(ValueError):
            est.mle(np.ones((2, 2, 4)), p, 50.0)


class TestMLEMemory(unittest.TestCase):
    """R1 defect: the global-grid chunk did not scale with the grid size G."""

    def test_fine_grid_bounded_memory_and_same_result(self):
        import tracemalloc
        p = tcp_model(50.0, sbr=10.0)
        rng = np.random.default_rng(42)
        C = rng.multinomial(2000, p(np.array([3.0, -2.0])), size=200)
        coarse = est.mle(C, p, search_radius=50.0)                        # G = 1961
        tracemalloc.start()
        fine = est.mle(C, p, search_radius=50.0, grid_step=0.2)           # G = 196321
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        # measured: ~135 MiB with the fix (64 MiB block + grid/model temporaries);
        # ~310 MiB with the R1 fixed chunk (200 x 196321 x 8 B = 314 MB for one block)
        self.assertLess(peak, 220 * 2 ** 20, peak)
        np.testing.assert_allclose(fine, coarse, atol=2e-3)

    def test_chunking_does_not_change_the_result(self):
        p = tcp_model(50.0)
        rng = np.random.default_rng(42)
        C = rng.multinomial(300, p(np.array([4.0, 1.0])), size=37)
        a = est.mle(C, p, search_radius=30.0, refine=False)
        b = est.mle(C, p, search_radius=30.0, refine=False, mem_budget=1.0)   # chunk = 1
        c = est.mle(C, p, search_radius=30.0, refine=False, chunk=5)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(a, c)


class TestLMS(unittest.TestCase):

    def test_lms_tcp_equals_general_lms(self):
        C = np.random.default_rng(42).multinomial(200, [0.3, 0.3, 0.35, 0.05], size=6)
        for L, fwhm, sbr, rot in ((50.0, 300.0, None, np.pi / 2), (100.0, 360.0, None, 0.0),
                                  (100.0, 300.0, 10.0, np.pi / 2), (70.0, np.inf, 5.0, 1.0)):
            fw = 1e9 if np.isinf(fwhm) else fwhm
            p = tcp_model(L, fw, sbr, rotation=rot)
            a = est.lms(C, p, h=1e-3)
            b = est.lms_tcp(C, L, fwhm, sbr=sbr, rotation=rot)
            np.testing.assert_allclose(a, b, rtol=1e-5, atol=1e-6)

    def test_lms_unbiased_to_first_order_including_background(self):
        # expected counts at a small displacement -> LMS returns ~ the displacement
        for sbr in (None, 10.0):
            p = tcp_model(50.0, 300.0, sbr)
            r = np.array([0.4, -0.3])
            got = est.lms_tcp(p(r) * 1000.0, 50.0, 300.0, sbr=sbr)
            np.testing.assert_allclose(got, r, rtol=2e-2)
            # without the 1/s factor the estimate would shrink by s = SBR/(SBR+1)
        naive = est.lms_tcp(tcp_model(50.0, 300.0, 10.0)(r) * 1000.0, 50.0, 300.0, sbr=None)
        np.testing.assert_allclose(naive, r * 10.0 / 11.0, rtol=2e-2)

    def test_lms_zero_counts_and_single(self):
        p = tcp_model(50.0)
        np.testing.assert_allclose(est.lms(np.zeros(4), p, r_lin=(1.0, 2.0)), [1.0, 2.0])
        np.testing.assert_allclose(est.lms_tcp(np.zeros(4), 50.0, 300.0), [0.0, 0.0])
        self.assertEqual(est.lms_tcp([1, 2, 3, 4], 50.0, 300.0).shape, (2,))

    def test_lms_tcp_errors(self):
        with self.assertRaises(ValueError):
            est.lms_tcp([1, 2, 3], 50.0, 300.0)
        with self.assertRaises(ValueError):
            est.lms_tcp([1, 2, 3, 4], 400.0, 300.0)   # L^2 ln2 / fwhm^2 > 1

    def test_mlms(self):
        C = np.random.default_rng(42).multinomial(100, [0.3, 0.3, 0.3, 0.1], size=5)
        base = est.lms_tcp(C, 50.0, 300.0)
        p0 = C[:, 3] / C.sum(1)
        np.testing.assert_allclose(est.mlms_tcp(C, 50.0, 300.0),
                                   (1.27 + 3.8 * p0)[:, None] * base, rtol=1e-12)
        np.testing.assert_allclose(est.mlms_tcp(C, 50.0, 300.0, beta=(1.0,)), base, rtol=1e-12)
        self.assertEqual(est.mlms_tcp(C[0], 50.0, 300.0).shape, (2,))


if __name__ == "__main__":
    unittest.main()
