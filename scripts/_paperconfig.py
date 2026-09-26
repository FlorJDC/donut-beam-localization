# -*- coding: utf-8 -*-
"""Canonical parameters shared by compute_paper_numbers.py and every fig_*.py (nm units)."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
SEED = 42
FWHM = 300.0                     # LG donut size parameter (Balzarotti Eq. S17)
FWHM_MASULLO = 360.0             # reproduces Masullo table
WAVELENGTH, NA, N_MEDIUM, FILLING = 640.0, 1.4, 1.518, 5.0 / 3.0
L_REF, N_REF = 50.0, 100
SBR_MLE = 10.0                   # mle_efficiency_center convention
L_FIT = (5.0, 10.0, 20.0, 40.0)  # crb_exponent_L: L << fwhm
N_FIT = (100, 300, 1000, 3000, 10000)
L_LIST = (50.0, 100.0, 150.0)
N_REP_MLE = 10000
SIGMA_PSF = 100.0                # camera: ideal sigma_PSF/sqrt(N)
CAM_PIXEL, CAM_NPIX = 100.0, 9
ITER = dict(n_iter=4, L0=150.0, L_min=25.0, photon_split="equal", recenter=True, rule="fixed")
ITER_N_TOTAL = 1000              # iterative_sigma_nm / camera_sigma_nm at the same photons
ITER_N_LIST = (250, 500, 1000, 2000, 4000, 8000)
ITER_N_REP = 10000
N_BOOT = 2000
EPS_LIST = (0.002, 0.01, 0.05, 0.15)
MIS = dict(L=100.0, N=500, sbr=10, n_patterns=400, n_rep=200)
MIS_DELTAS = (0.0, 2.0, 5.0, 10.0)
DATA = os.path.join(ROOT, "data"); MC = os.path.join(DATA, "mc"); FIGDIR = os.path.join(ROOT, "paper", "figures")
# MLE search-disk radii (fig 5 and compute_paper_numbers must agree): centre points use radius
# L (r1 convention of mle_efficiency_center), off-centre / x-sweep points use 2L.
MLE_RADIUS_CENTRE_OVER_L = 1.0
MLE_RADIUS_SWEEP_OVER_L = 2.0
N_REP_BIAS0 = 40000              # MLE bias near the centre without background (mle_nobg_bias_x_r2_nm)
EPS_TRANSITION = 0.002           # eps of zero_depth_transition_scale_nm (fig 7c, V10)
# noise-free (population) naive-MLE bias under misalignment (inbox r4)
MIS_POP_N_PATTERNS = 4000
MIS_POP_DELTAS = (2.0, 5.0, 10.0)
QUICK_DIRNAME = "quick"          # --quick outputs go to paper/figures/quick, data/quick, ...
