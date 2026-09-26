# -*- coding: utf-8 -*-
"""Verifier r05: independent numeric checks of the r4 demo keys and a few text formulas
(own Fisher matrix by central finite differences; does NOT import donutloc)."""
import numpy as np

FWHM = 300.0
A = 4 * np.log(2) / FWHM ** 2


def ring(Lv, rot=np.pi / 2):
    ang = rot + 2 * np.pi * np.arange(3) / 3
    b = np.zeros((4, 2))
    b[:3, 0] = Lv / 2 * np.cos(ang)
    b[:3, 1] = Lv / 2 * np.sin(ang)
    return b


def beam(u, eps=0.0, model="constant"):
    if model == "constant":
        return np.e * A * u * np.exp(-A * u) + eps
    return (np.e * A * u + eps) * np.exp(-A * u)


def p_of(r, b, eps=0.0, model="constant", sbr=None, bg_abs=None):
    u = np.sum((np.asarray(r)[None, :] - b) ** 2, axis=1)
    lam = beam(u, eps, model)
    if bg_abs is not None:
        lam = lam + bg_abs
    p = lam / lam.sum()
    if sbr is not None:
        s = sbr / (sbr + 1)
        p = s * p + (1 - s) / 4
    return p


def crb(r, b, N, h=1e-4, **kw):
    r = np.asarray(r, float)
    p = p_of(r, b, **kw)
    g = np.zeros((4, 2))
    for k in range(2):
        e = np.zeros(2); e[k] = h
        g[:, k] = (p_of(r + e, b, **kw) - p_of(r - e, b, **kw)) / (2 * h)
    m = p > 0
    F = N * (g[m].T / p[m]) @ g[m]
    return np.sqrt(np.trace(np.linalg.inv(F)) / 2)


def s31(Lv, N, sbr):
    x = Lv ** 2 * np.log(2) / FWHM ** 2
    return Lv / (2 * np.sqrt(2 * N)) / (1 - x) * np.sqrt((1 + 1 / sbr) * (1 + 3 / (4 * sbr)))


Lv, N, eps = 100.0, 100, 0.05
b = ring(Lv)
x = Lv ** 2 * np.log(2) / FWHM ** 2
sbr_eps = 3 * np.e * x * np.exp(-x) / (4 * eps)
# also from the definition: sum of LG at centre / (4 eps)
sum_lg = np.sum(beam(np.sum(b ** 2, 1), 0.0))
print("SBR_eps closed %.6f   from sum I_LG(0)/(4 eps) %.6f" % (sbr_eps, sum_lg / (4 * eps)))
c_const = crb([0, 0], b, N, eps=eps, model="constant")
c_gauss = crb([0, 0], b, N, eps=eps, model="gaussian")
ref = s31(Lv, N, sbr_eps)
print("constant pedestal CRB %.6f / S31(SBR_eps) %.6f = %.7f" % (c_const, ref, c_const / ref))
print("gaussian pedestal CRB %.6f / S31(SBR_eps) %.6f = %.7f" % (c_gauss, ref, c_gauss / ref))
# point (r=0) vs small r for the pedestal (continuity)
for rr in (1e-3, 1e-2):
    print("  r=%g: const %.6f gauss %.6f" % (rr, crb([rr, 0], b, N, eps=eps, model="constant"),
                                            crb([rr, 0], b, N, eps=eps, model="gaussian")))
# degradation table checks vs S27
s27 = Lv / (2 * np.sqrt(2 * N)) / (1 - x)
print("degradation L100 eps0.05: const %.5f gauss %.5f (registry 1.30023 / 1.30715)"
      % (c_const / s27, c_gauss / s27))
# combined background, 'beam' convention: 1+1/SBR_eff = (1+1/SBR)(1+1/SBR_eps)
sbr = 10.0
lam_b = (sum_lg + 4 * eps) / (4 * sbr)
c_comb = crb([0, 0], b, N, eps=eps, model="constant", bg_abs=lam_b)
sbr_eff = 1 / ((1 + 1 / sbr) * (1 + 1 / sbr_eps) - 1)
print("const pedestal + SBR=10 (beam conv.): numeric %.6f  S31(SBR_eff) %.6f" % (c_comb, s31(Lv, N, sbr_eff)))

# LMS shrink factor (Balzarotti Eq. S50 with rotation pi/2, own implementation)
L2 = 50.0
b2 = ring(L2)
g2 = 1 - L2 ** 2 * np.log(2) / FWHM ** 2
for r0 in ([5, 0], [0, 3], [4, -2], [20, 7]):
    pb = p_of(r0, b2, sbr=10.0)
    lms_raw = -(1 / g2) * pb[:3] @ b2[:3]
    lms_cor = lms_raw / (10 / 11)
    p0 = p_of(r0, b2)
    lms_nobg = -(1 / g2) * p0[:3] @ b2[:3]
    print("LMS at %s: raw/corrected = %.7f ; raw/background-free LMS = %s"
          % (r0, lms_raw[0] / lms_cor[0] if lms_cor[0] else np.nan, lms_raw / lms_nobg))

# small checks of text statements
L3 = 50.0
x3 = L3 ** 2 * np.log(2) / FWHM ** 2
g3 = 1 - x3
alpha = 8 * 100 * g3 ** 2 / L3 ** 2
beta = 16 * 100 * np.exp(x3) / (3 * L3 ** 2)
print("L=50 limit axes par %.4f perp %.4f  sigma_lim %.5f" % ((alpha + beta) ** -.5, alpha ** -.5,
      np.sqrt((1 / alpha + 1 / (alpha + beta)) / 2)))
num_lim = np.mean([crb([1e-3 * np.cos(t), 1e-3 * np.sin(t)], ring(L3), 100)
                   for t in np.linspace(0, 2 * np.pi, 12, endpoint=False)])
print("numeric limit (12 dirs, r=1e-3) %.5f ; point value %.5f" % (num_lim, crb([0, 0], ring(L3), 100)))
# ratios wrong-hand / correct etc. from Table I
T = {50: (1.6, 147.9, 60.77), 100: (3.323, 89.59, 38.06), 150: (5.331, 84.2, 35.72)}
for k, (c, o, l) in T.items():
    print("L=%d opposite/correct %.1f linear/correct %.1f" % (k, o / c, l / c))
# r_c at SBR 10, 5 as fraction of L
for s in (5, 10, 20):
    print("r_c/L at SBR=%d: %.3f (L=50) %.3f (L=100)" % (s, np.sqrt(3) / 4 * np.exp(-x3 / 2) / np.sqrt(s),
          np.sqrt(3) / 4 * np.exp(-x / 2) / np.sqrt(s)))
# honest bias vs floor
print("honest centre excess over floor: %.3f SE ; ratio %.3f" % ((0.1797 - 0.1643) / 0.005103, 0.1797 / 0.1643))
print("honest Lq excess over floor: %.3f SE ; ratio %.3f" % ((0.2759 - 0.2309) / 0.0083, 0.2759 / 0.2309))
# MC weighted slope using the population ratios and the MC SE weights
se = np.array([0.03262, 0.08025, 0.1583]); dl = np.array([2., 5., 10.]); pr = np.array([0.7471, 0.7549, 0.7835])
w = 1 / se ** 2
print("population ratios with MC weights -> slope %.4f" % (np.sum(w * dl * pr * dl) / np.sum(w * dl ** 2)))
print("Masullo L=150 fwhm360: %.4f" % (150 / (2 * np.sqrt(1000)) / (1 - 150 ** 2 * np.log(2) / 360 ** 2)
                                       * np.sqrt(1.2 * 1.15)))
