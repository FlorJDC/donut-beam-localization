# -*- coding: utf-8 -*-
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 5))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from scipy.stats import norm
from donutloc import camera as cam, fisher

p = cam.make_camera_model(100., 100., 3)
q = p(np.array([60., -130.]))
# independent pixel probabilities via norm.cdf, row-major y outer
e = cam.pixel_edges(100., 3)
qx = np.diff(norm.cdf(e, loc=60., scale=100.)); qy = np.diff(norm.cdf(e, loc=-130., scale=100.))
ref = np.outer(qy, qx).ravel(); ref /= ref.sum()
print("pixel order/values max err:", np.max(abs(q - ref)))
# asymmetric check: emitter at +x -> larger p in last column
print("row-major check (x=+150): p reshaped\n", p(np.array([150., 0.])).reshape(3, 3).round(3))
# ideal limit & N scaling
print("small pix big window / ideal:", cam.crb_camera(100., 400, pixel=5., n_pix=241) / 5.0)
print("Balzarotti p1 9x9 a=100 N=400:", cam.crb_camera(100., 400))
# Independent analytic Fisher for pixelated camera without bg (1D product structure) at r=0
def fisher1d(a, n, s, x0=0.):
    e = a * (np.arange(n + 1) - n / 2)
    Phi = norm.cdf(e, x0, s); phi = norm.pdf(e, x0, s)
    q = np.diff(Phi); dq = -np.diff(phi)
    return q, dq
qx, dqx = fisher1d(100., 9, 100.)
# full 2D with renormalization: p = qx qy / (Sx Sy); window renorm
Sx = qx.sum()
px = qx / Sx; dpx = dqx / Sx - qx * dqx.sum() / Sx ** 2
Fxx = 400 * np.sum(dpx ** 2 / px)  # separable: p_ij = px_i py_j -> F_xx = N sum dpx^2/px
print("independent separable CRB (renormalized):", 1 / np.sqrt(Fxx))
# window edge: emitter near window edge -> renormalized model
for x in (0., 300., 450., 600., 2000., 1e5):
    print("r=(%g,0): crb %s" % (x, cam.crb_camera(100., 400, r=(x, 0.))))
# conventions
a = cam.crb_camera(100., 600, sbr=500., sbr_convention="per_pixel"); b = cam.crb_camera(100., 600, sbr=500. / 81, sbr_convention="total")
print("SBR_c=500 per pixel vs total 6.17:", a, b)
# sbr=0
print("sbr=0 crb:", cam.crb_camera(100., 400, sbr=0.))
# p_min exclusion effect: small-pixel far tails
print("n_pix=1 crb:", cam.crb_camera(100., 400, n_pix=1))
