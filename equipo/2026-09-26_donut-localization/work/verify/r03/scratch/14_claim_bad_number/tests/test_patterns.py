import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np

from donutloc import patterns


def _tcp_reference(L):
    # copy of tests/test_acceptance.py::_tcp (not imported on purpose)
    ang = np.pi / 2 + 2 * np.pi * np.arange(3) / 3
    ring = np.stack([0.5 * L * np.cos(ang), 0.5 * L * np.sin(ang)], axis=1)
    return np.vstack([ring, [[0.0, 0.0]]])


class TestPatterns(unittest.TestCase):

    def test_tcp_matches_acceptance_convention_exactly(self):
        for L in (50.0, 100.0, 150.0, 37.3):
            np.testing.assert_array_equal(patterns.tcp_centers(L), _tcp_reference(L))

    def test_tcp_geometry(self):
        L = 100.0
        c = patterns.tcp_centers(L)
        self.assertEqual(c.shape, (4, 2))
        np.testing.assert_array_equal(c[3], [0.0, 0.0])
        np.testing.assert_allclose(np.hypot(c[:3, 0], c[:3, 1]), L / 2, rtol=1e-15)
        np.testing.assert_allclose(c[0], [0.0, L / 2], atol=1e-12)     # vertex up
        np.testing.assert_allclose(c[:3].sum(axis=0), 0.0, atol=1e-12)
        self.assertEqual(patterns.tcp_centers(L, center=False).shape, (3, 2))

    def test_polygon(self):
        c = patterns.polygon_centers(80.0, 6, rotation=0.0)
        self.assertEqual(c.shape, (7, 2))
        np.testing.assert_allclose(c[0], [40.0, 0.0], atol=1e-12)
        np.testing.assert_allclose(np.hypot(c[:6, 0], c[:6, 1]), 40.0, rtol=1e-15)
        np.testing.assert_array_equal(patterns.polygon_centers(50.0, 3), patterns.tcp_centers(50.0))
        with self.assertRaises(ValueError):
            patterns.polygon_centers(50.0, 0)
        with self.assertRaises(ValueError):
            patterns.polygon_centers(-1.0, 3)

    def test_perturb_random_direction_magnitude(self):
        c = patterns.tcp_centers(100.0)
        p = patterns.perturb_centers(c, 5.0, rng=42)
        np.testing.assert_allclose(np.linalg.norm(p - c, axis=1), 5.0, rtol=1e-12)
        self.assertFalse(np.allclose(p, c))
        # reproducible with seed; input not modified
        np.testing.assert_array_equal(p, patterns.perturb_centers(c, 5.0, rng=42))
        np.testing.assert_array_equal(c, patterns.tcp_centers(100.0))
        # per-center magnitudes
        mags = np.array([1.0, 2.0, 0.0, 3.0])
        p2 = patterns.perturb_centers(c, mags, rng=np.random.default_rng(42))
        np.testing.assert_allclose(np.linalg.norm(p2 - c, axis=1), mags, atol=1e-12)

    def test_perturb_directions_uniform(self):
        rng = np.random.default_rng(42)
        c = np.zeros((20000, 2))
        d = patterns.perturb_centers(c, 1.0, rng=rng)
        self.assertLess(abs(d[:, 0].mean()), 0.03)
        self.assertLess(abs(d[:, 1].mean()), 0.03)
        self.assertAlmostEqual((d[:, 0] ** 2).mean(), 0.5, delta=0.02)

    def test_perturb_explicit(self):
        c = patterns.tcp_centers(100.0)
        d = np.arange(8.0).reshape(4, 2)
        np.testing.assert_array_equal(patterns.perturb_centers(c, d), c + d)

    def test_perturb_errors(self):
        c = patterns.tcp_centers(100.0)
        with self.assertRaises(ValueError):
            patterns.perturb_centers(c, np.ones(3))
        with self.assertRaises(ValueError):
            patterns.perturb_centers(c, 1.0, mode="bogus")
        with self.assertRaises(ValueError):
            patterns.perturb_centers(np.ones(4), 1.0)


if __name__ == "__main__":
    unittest.main()
