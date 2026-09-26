# -*- coding: utf-8 -*-
"""Tests for donutloc.camera (ideal camera, Balzarotti2017 Eq. S59-S63)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np
from scipy.special import erf

from donutloc import camera as cam
from donutloc import fisher


class TestModel(unittest.TestCase):

    def test_shape_normalisation_and_pixel_values(self):
        p = cam.make_camera_model(100.0, 100.0, 9)
        self.assertEqual(p.K, 81)
        P = p(np.zeros((3, 4, 2)))
        self.assertEqual(P.shape, (3, 4, 81))
        np.testing.assert_allclose(P.sum(-1), 1.0, atol=1e-12)
        # centre pixel (index 40) against the erf formula, renormalised in the window
        q = lambda lo, hi: 0.5 * (erf(hi / (np.sqrt(2) * 100.0)) - erf(lo / (np.sqrt(2) * 100.0)))
        tot = q(-450.0, 450.0) ** 2
        self.assertAlmostEqual(float(p([0.0, 0.0])[40]), q(-50.0, 50.0) ** 2 / tot, places=12)
        # row-major: pixel index = iy * n + ix; emitter at +x shifts weight to ix = 5
        v = p([100.0, 0.0])
        self.assertEqual(int(np.argmax(v)), 4 * 9 + 5)

    def test_edges_centered(self):
        np.testing.assert_allclose(cam.pixel_edges(100.0, 3), [-150, -50, 50, 150])
        np.testing.assert_allclose(cam.pixel_edges(10.0, 2), [-10, 0, 10])

    def test_background_conventions(self):
        p = cam.make_camera_model(sbr=500.0, sbr_convention="per_pixel")
        q = cam.make_camera_model(sbr=500.0 / 81.0, sbr_convention="total")
        r = np.array([[0.0, 0.0], [37.0, -12.0]])
        np.testing.assert_allclose(p(r), q(r), rtol=1e-13)
        p0 = cam.make_camera_model()(r)
        # Eq. S62-S63
        np.testing.assert_allclose(p(r), 1.0 / (81 + 500) + 500.0 / (81 + 500) * p0, rtol=1e-13)
        np.testing.assert_allclose(cam.make_camera_model(sbr=0.0)(r), 1.0 / 81)

    def test_far_field_nan_and_errors(self):
        self.assertTrue(np.all(np.isnan(cam.make_camera_model()([1e5, 0.0]))))
        for kw in ({"sigma_psf": 0.0}, {"pixel": -1.0}, {"n_pix": 0}, {"sbr": -1.0},
                   {"sbr": 5.0, "sbr_convention": "bogus"}):
            with self.assertRaises(ValueError, msg=str(kw)):
                cam.make_camera_model(**kw)


class TestCRB(unittest.TestCase):

    def test_small_pixels_large_window_reach_ideal(self):
        got = cam.crb_camera(100.0, 400, pixel=10.0, n_pix=121)
        self.assertAlmostEqual(got / cam.crb_camera_ideal(100.0, 400), 1.0, delta=0.01)

    def test_balzarotti_p1_check(self):
        self.assertAlmostEqual(float(cam.crb_camera_ideal(100.0, 400)), 5.0, places=12)
        got = cam.crb_camera(100.0, 400)                     # a = 100 nm, 9 x 9, no background
        self.assertTrue(5.0 < got < 5.5, got)

    def test_pixelation_and_background_worsen(self):
        fine = cam.crb_camera(100.0, 400, pixel=10.0, n_pix=121)
        coarse = cam.crb_camera(100.0, 400, pixel=100.0, n_pix=9)
        self.assertGreater(coarse, fine * 1.02)
        vals = [cam.crb_camera(100.0, 400, sbr=s) for s in (None, 100.0, 10.0, 5.0, 1.0)]
        self.assertTrue(np.all(np.diff(vals) > 0), vals)

    def test_sbr_conventions_note_A_section9(self):
        a = cam.crb_camera(100.0, 600, sbr=500.0, sbr_convention="per_pixel")
        b = cam.crb_camera(100.0, 600, sbr=500.0 / 81.0)     # SBR_total = 6.17
        self.assertAlmostEqual(a / b, 1.0, delta=1e-12)
        self.assertAlmostEqual(500.0 / 81.0, 6.17, delta=0.01)

    def test_inverse_sqrt_N_and_array_r(self):
        a = cam.crb_camera(100.0, 100, sbr=10.0)
        b = cam.crb_camera(100.0, 400, sbr=10.0)
        self.assertAlmostEqual(a / b, 2.0, delta=1e-9)
        out = cam.crb_camera(100.0, 400, r=np.array([[0.0, 0.0], [50.0, 0.0]]))
        self.assertEqual(out.shape, (2,))
        p = cam.make_camera_model()
        self.assertAlmostEqual(out[1] / float(fisher.crb(p, [50.0, 0.0], 400)), 1.0, delta=1e-12)


if __name__ == "__main__":
    unittest.main()
