# -*- coding: utf-8 -*-
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import unittest

import numpy as np
from scipy.special import jv

from donutloc import vectorial as v
from donutloc import patterns, photons, beams

NB = dict(n=1.5)  # parameters of docs/literature/B_donut_optics.md (NA 1.4, 640 nm, F = 5/3)


class TestZeroDepth(unittest.TestCase):

    def test_correct_hand_perfect_zero(self):
        self.assertLess(v.zero_depth(handedness=+1, **NB), 1e-10)
        self.assertLess(v.zero_depth(handedness=+1), 1e-10)
        # generalized rule l*s = +1: charge -1 with s = -1 is also a perfect zero
        self.assertLess(v.zero_depth(charge=-1, handedness=-1, **NB), 1e-10)
        # and the zero holds off focus (J1, J2, J3 vanish on axis for every z)
        self.assertEqual(float(v.intensity(0.0, 0.0, z=300.0, **NB)), 0.0)

    def test_opposite_hand(self):
        self.assertAlmostEqual(v.zero_depth(handedness=-1, **NB), 0.869, delta=0.01)

    def test_linear(self):
        self.assertAlmostEqual(v.zero_depth(polarization="linear", **NB), 0.382, delta=0.01)

    def test_opposite_hand_bessel_orders(self):
        # on axis only Ez survives and equals -R(0)/sqrt2 with R(0) = int g sin(th) (J0(0) = 1)
        Ex, Ey, Ez = v.focal_field(np.array([0.0]), np.array([0.0]), handedness=-1, **NB)
        self.assertEqual(abs(Ex[0]), 0.0)
        self.assertEqual(abs(Ey[0]), 0.0)
        k, th, w, a = v._params(n=1.5)
        R0 = np.sum(w * a * np.sqrt(np.cos(th)) * np.sin(th) ** 2)
        self.assertAlmostEqual(abs(Ez[0]) / (R0 / np.sqrt(2.0)), 1.0, places=12)
        # transverse part is pure J1: |Ex|^2+|Ey|^2 = (|I_A|^2 + |I_B'|^2)/4
        rho = 120.0
        u = k * rho * np.sin(th)
        g = w * a * np.sqrt(np.cos(th)) * np.sin(th)
        IA = np.sum(g * (1 + np.cos(th)) * jv(1, u))
        IB = np.sum(g * (1 - np.cos(th)) * jv(1, u))
        Ex, Ey, Ez = v.focal_field(rho, 0.7, handedness=-1, **NB)
        tr = abs(Ex) ** 2 + abs(Ey) ** 2
        self.assertAlmostEqual(tr / (0.25 * (IA ** 2 + IB ** 2)), 1.0, places=12)
        IC = np.sum(g * np.sin(th) * jv(0, u))
        self.assertAlmostEqual(abs(Ez) / (abs(IC) / np.sqrt(2.0)), 1.0, places=12)


class TestShape(unittest.TestCase):

    def test_peak_to_peak(self):
        self.assertAlmostEqual(v.peak_to_peak_diameter(filling=100.0, **NB) / 380.0, 1.0, delta=0.03)
        self.assertAlmostEqual(v.peak_to_peak_diameter(filling=2.0, **NB) / 392.0, 1.0, delta=0.03)

    def test_curvature(self):
        c = v.zero_curvature(**NB)
        self.assertAlmostEqual(c / 7.0e-5, 1.0, delta=0.05)
        # equals the analytic small-rho limit (k^2/16)|int g (1+cos) sin|^2 / I_max
        k, th, w, a = v._params(n=1.5)
        g = w * a * np.sqrt(np.cos(th)) * np.sin(th)
        Imax = v._max_search(dict(NB))[0]
        c_an = k ** 2 / 16.0 * np.sum(g * (1 + np.cos(th)) * np.sin(th)) ** 2 / Imax
        self.assertAlmostEqual(c / c_an, 1.0, places=5)

    def test_rotational_symmetry_circular(self):
        rho = np.array([10.0, 80.0, 190.0, 400.0])
        for h in (+1, -1):
            I = [v.intensity(rho * np.cos(p), rho * np.sin(p), handedness=h, **NB) for p in (0., 1., 2.)]
            for Ii in I[1:]:
                np.testing.assert_allclose(Ii, I[0], rtol=1e-10, atol=0)
            np.testing.assert_allclose(v.radial_profile(rho, handedness=h, **NB), I[0], rtol=1e-12)

    def test_linear_not_symmetric_and_radial_profile_refuses(self):
        a = v.intensity(150.0, 0.0, polarization="linear", **NB)
        b = v.intensity(0.0, 150.0, polarization="linear", **NB)
        self.assertGreater(abs(a - b) / max(a, b), 1e-3)
        with self.assertRaises(ValueError):
            v.radial_profile(np.array([1.0]), polarization="linear")

    def test_bad_arguments(self):
        with self.assertRaises(ValueError):
            v.focal_field(1.0, 0.0, handedness=0)
        with self.assertRaises(ValueError):
            v.focal_field(1.0, 0.0, polarization="radial")
        with self.assertRaises(ValueError):
            v.focal_field(1.0, 0.0, NA=1.6, n=1.5)
        with self.assertRaises(ValueError):
            v.focal_field(1.0, 0.0, method="quad")
        with self.assertRaises(TypeError):
            v.intensity(1.0, 0.0, fwhm=3.0)
        with self.assertRaises(ValueError):
            v.lg_equivalent_fwhm({}, match="area")


class TestAgainst2D(unittest.TestCase):

    def test_bessel_equals_direct_2d(self):
        rho = np.array([0.0, 40.0, 150.0, 330.0])
        phi = np.array([0.3, 1.1, -2.0, 2.9])
        cases = [dict(polarization="circular", handedness=+1),
                 dict(polarization="circular", handedness=-1),
                 dict(polarization="linear", pol_angle=0.4)]
        for c in cases:
            for z in (0.0, 250.0):
                A = v.focal_field(rho, phi, z=z, n_theta=201, **dict(c, **NB))
                B = v.focal_field(rho, phi, z=z, n_theta=201, method="2d", n_phi=200, **dict(c, **NB))
                sc = max(np.abs(x).max() for x in A)
                for x, y in zip(A, B):
                    self.assertLess(np.abs(x - y).max() / sc, 1e-6, msg=str((c, z)))


class TestEnergy(unittest.TestCase):

    def test_flux_same_for_both_hands(self):
        # int |E|^2 2 pi rho drho at z = 0 is fixed by Parseval, (2pi/k^2) int a^2 sin(th) dth,
        # independent of s (|E_0|^2 = a^2 for both hands): the total (incl. Ez) must agree.
        r = np.arange(0.0, 2000.0, 2.0)
        k, th, w, a = v._params(n_theta=401)
        flux = {}
        for h in (+1, -1):
            I = v.radial_profile(r, handedness=h, n_theta=401)
            _trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy>=2 renamed trapz
            flux[h] = _trapz(I * 2 * np.pi * r, r)
        self.assertLess(abs(flux[1] / flux[-1] - 1.0), 1e-3)
        pars = 2 * np.pi / k ** 2 * np.sum(w * a ** 2 * np.sin(th))
        self.assertAlmostEqual(flux[1] / pars, 1.0, delta=0.05)  # truncation tail ~ 1/R


class TestBeam(unittest.TestCase):

    def test_beam_peak_and_model(self):
        b = v.make_vectorial_beam(rho_max=700.0, **NB)
        r = np.linspace(0, 600, 6001)
        vals = b(r, 0 * r)
        self.assertAlmostEqual(vals.max(), 1.0, delta=1e-4)  # 1 nm table step
        self.assertEqual(float(b(0.0, 0.0)), 0.0)
        self.assertEqual(float(b(2000.0, 0.0)), 0.0)
        # linear in rho^2 below the first node: I/rho^2 -> curvature
        c = v.zero_curvature(**NB)
        self.assertAlmostEqual(float(b(1e-3, 0.0)) / 1e-6 / c, 1.0, delta=1e-3)
        p = photons.make_model(patterns.tcp_centers(100.0), b)
        pts = np.array([[5.0, -3.0], [30.0, 20.0], [0.2, 0.0]])
        np.testing.assert_allclose(p(pts).sum(axis=-1), 1.0, rtol=1e-12)
        be = v.make_vectorial_beam(rho_max=100.0, eps=0.01, **NB)
        self.assertAlmostEqual(float(be(0.0, 0.0)), 0.01, places=12)
        with self.assertRaises(ValueError):
            v.make_vectorial_beam(eps=-1.0)

    def test_linear_beam_grid(self):
        b = v.make_vectorial_beam(rho_max=300.0, grid_step=5.0, polarization="linear", n_theta=201, **NB)
        Imax, rpk, ppk = v._max_search(dict(polarization="linear", n_theta=201, **NB))
        self.assertAlmostEqual(float(b(rpk * np.cos(ppk), rpk * np.sin(ppk))), 1.0, delta=2e-3)
        self.assertAlmostEqual(float(b(0.0, 0.0)), v.zero_depth(polarization="linear", n_theta=201, **NB),
                               places=10)
        p = photons.make_model(patterns.tcp_centers(50.0), b)
        self.assertAlmostEqual(float(p(np.array([3.0, 4.0])).sum()), 1.0, places=12)

    def test_lg_equivalent_and_crb(self):
        c = v.zero_curvature(**NB)
        fc = v.lg_equivalent_fwhm(NB, match="curvature")
        self.assertAlmostEqual(4 * np.e * np.log(2) / fc ** 2 / c, 1.0, places=10)
        fd = v.lg_equivalent_fwhm(NB, match="diameter")
        self.assertAlmostEqual(2 * beams.ring_radius(fd) / v.peak_to_peak_diameter(**NB), 1.0, places=10)
        out = v.compare_crb_vectorial_vs_lg([50.0], N=100, **NB)
        # near the zero the curvature-matched LG must reproduce the vectorial CRB at small L
        self.assertAlmostEqual(out["crb_vectorial"][0] / out["crb_lg_curvature"][0], 1.0, delta=0.01)


class TestLinearBeamHarmonic(unittest.TestCase):
    """R2 refuted: bilinear Cartesian table biased the linear-pol CRB (+7.4 % at (7,3)).
    Reference (n = 1.518, L = 50, N = 100, mode='exact'): CRB(7,3) = 62.36, crb_limit = 60.77."""

    @classmethod
    def setUpClass(cls):
        from donutloc import fisher
        cls.fisher = fisher
        o = dict(polarization="linear", n_theta=401)
        cls.bh = staticmethod(v.make_vectorial_beam(rho_max=120.0, **o))
        cls.be = staticmethod(v.make_vectorial_beam(mode="exact", **o))
        cls.c = patterns.tcp_centers(50.0)

    def test_harmonic_matches_exact_intensity(self):
        rng = np.random.default_rng(0)
        pts = rng.uniform(-80, 80, size=(50, 2))
        np.testing.assert_allclose(self.bh(pts[:, 0], pts[:, 1]), self.be(pts[:, 0], pts[:, 1]),
                                   rtol=1e-6)
        self.assertEqual(float(self.bh(200.0, 0.0)), 0.0)
        self.assertEqual(self.bh.linear_method, "harmonic")

    def test_crb_within_half_percent_of_exact(self):
        f = self.fisher
        ph = photons.make_model(self.c, self.bh)
        pe = photons.make_model(self.c, self.be)
        r = np.array([7.0, 3.0])
        a, b = float(f.crb(ph, r, 100)), float(f.crb(pe, r, 100))
        la, lb = f.crb_limit(ph, 100), f.crb_limit(pe, 100)
        self.assertAlmostEqual(a / b, 1.0, delta=5e-3)
        self.assertAlmostEqual(la / lb, 1.0, delta=5e-3)
        self.assertAlmostEqual(b, 62.36, delta=0.1)
        self.assertAlmostEqual(lb, 60.77, delta=0.1)

    def test_cartesian_warns_and_is_biased(self):
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bc = v.make_vectorial_beam(rho_max=120.0, grid_step=5.0, polarization="linear",
                                       n_theta=401, linear_method="cartesian")
        self.assertTrue(any(issubclass(x.category, UserWarning) for x in w))
        f = self.fisher
        a = float(f.crb(photons.make_model(self.c, bc), np.array([7.0, 3.0]), 100))
        b = float(f.crb(photons.make_model(self.c, self.be), np.array([7.0, 3.0]), 100))
        self.assertGreater(abs(a / b - 1.0), 0.02)
        with self.assertRaises(ValueError):
            v.make_vectorial_beam(rho_max=50.0, polarization="linear", n_theta=101,
                                  linear_method="bogus")


if __name__ == "__main__":
    unittest.main()
