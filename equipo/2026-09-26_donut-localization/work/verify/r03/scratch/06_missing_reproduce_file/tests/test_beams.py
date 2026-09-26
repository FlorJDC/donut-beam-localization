import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np
from scipy.optimize import minimize_scalar

from donutloc import beams


class TestBeams(unittest.TestCase):

    def test_ring_peak_is_one_at_ring_radius(self):
        for fwhm in (300.0, 360.0, 200.0):
            rr = beams.ring_radius(fwhm)
            self.assertAlmostEqual(float(beams.lg_donut(rr, 0.0, fwhm)), 1.0, places=14)
            # independent: numerical maximum of the radial profile
            res = minimize_scalar(lambda r: -beams.lg_donut(r, 0.0, fwhm),
                                  bounds=(1.0, 2 * fwhm), method="bounded",
                                  options={"xatol": 1e-9})
            self.assertAlmostEqual(res.x, rr, delta=1e-4)
            self.assertAlmostEqual(-res.fun, 1.0, places=10)

    def test_peak_to_peak_diameter_about_1p2_fwhm(self):
        fwhm = 300.0
        d = 2 * beams.ring_radius(fwhm)
        self.assertAlmostEqual(d / fwhm, 1.0 / np.sqrt(np.log(2)), places=12)
        self.assertAlmostEqual(d / fwhm, 1.2, delta=0.005)

    def test_eps_at_zero(self):
        for model in ("gaussian", "constant"):
            for e in (0.0, 0.002, 0.05, 0.15):
                self.assertEqual(float(beams.lg_donut(0, 0, eps=e, zero_model=model)), e)

    def test_eps_models_differ_away_from_zero(self):
        rr = beams.ring_radius(300.0)
        g = beams.lg_donut(rr, 0, eps=0.1, zero_model="gaussian")
        c = beams.lg_donut(rr, 0, eps=0.1, zero_model="constant")
        self.assertAlmostEqual(float(c), 1.1, places=12)
        self.assertAlmostEqual(float(g), 1.0 + 0.1 * np.exp(-1.0), places=12)  # a r_ring^2 = 1

    def test_quadratic_matches_donut_near_zero(self):
        fwhm = 300.0
        for r in (0.1, 1.0, 3.0):
            q = beams.quadratic(r, 0.0, fwhm)
            d = beams.lg_donut(r, 0.0, fwhm)
            a = 4 * np.log(2) / fwhm ** 2
            # relative difference is 1 - exp(-a r^2) ~ a r^2
            self.assertAlmostEqual(d / q, np.exp(-a * r * r), places=14)
            self.assertLess(abs(d / q - 1), 2 * a * r * r)

    def test_gaussian(self):
        fwhm = 300.0
        self.assertEqual(float(beams.gaussian(0, 0, fwhm)), 1.0)
        # half maximum at r = fwhm/2
        self.assertAlmostEqual(float(beams.gaussian(fwhm / 2, 0, fwhm)), 0.5, places=14)

    def test_vectorized_and_rotationally_symmetric(self):
        x = np.linspace(-200, 200, 7)[:, None]
        y = np.linspace(-150, 150, 5)[None, :]
        out = beams.lg_donut(x, y, 300.0, eps=0.01)
        self.assertEqual(out.shape, (7, 5))
        th = 0.7
        xr, yr = x * np.cos(th) - y * np.sin(th), x * np.sin(th) + y * np.cos(th)
        np.testing.assert_allclose(beams.lg_donut(xr, yr, 300.0, eps=0.01), out, rtol=1e-12)

    def test_make_beam(self):
        f = beams.make_beam("donut", fwhm=360.0, eps=0.03, zero_model="constant")
        x, y = np.array([0.0, 10.0, 150.0]), np.array([0.0, 5.0, -20.0])
        np.testing.assert_array_equal(f(x, y), beams.lg_donut(x, y, 360.0, 0.03, "constant"))
        self.assertEqual((f.kind, f.fwhm, f.eps, f.zero_model, f.power),
                         ("donut", 360.0, 0.03, "constant", 1.0))
        g = beams.make_beam("gaussian", fwhm=250.0)
        np.testing.assert_array_equal(g(x, y), beams.gaussian(x, y, 250.0))
        q = beams.make_beam("quadratic")
        np.testing.assert_array_equal(q(x, y), beams.quadratic(x, y, 300.0))

    def test_make_beam_power(self):
        f2 = beams.make_beam("donut", power=2)
        x, y = np.array([3.0, 100.0]), np.array([4.0, 0.0])
        np.testing.assert_allclose(f2(x, y), beams.lg_donut(x, y) ** 2, rtol=1e-15)
        f15 = beams.make_beam("donut", eps=0.01, power=1.5)
        np.testing.assert_allclose(f15(x, y), beams.lg_donut(x, y, eps=0.01) ** 1.5, rtol=1e-15)

    def test_errors(self):
        with self.assertRaises(ValueError):
            beams.lg_donut(0, 0, zero_model="bogus")
        with self.assertRaises(ValueError):
            beams.lg_donut(0, 0, eps=-0.1)
        with self.assertRaises(ValueError):
            beams.lg_donut(0, 0, fwhm=0.0)
        with self.assertRaises(ValueError):
            beams.make_beam("bessel")
        with self.assertRaises(ValueError):
            beams.make_beam("gaussian", eps=0.1)
        with self.assertRaises(ValueError):
            beams.make_beam("donut", zero_model="bogus")
        with self.assertRaises(ValueError):
            beams.make_beam("donut", power=0)
        with self.assertRaises(ValueError):
            beams.ring_radius(-1.0)


if __name__ == "__main__":
    unittest.main()
