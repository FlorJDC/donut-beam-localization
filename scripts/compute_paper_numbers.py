# -*- coding: utf-8 -*-
"""Compute every number the paper cites and write ``data/paper_numbers.json`` (the ONLY writer)
plus the LaTeX macros ``paper/generated/numbers.tex``.

Format of ``data/paper_numbers.json``::

    {key: {"value": float|int, "unit": str, "script": str, "description": str[, "se": float]}}

``se`` is the Monte Carlo standard error (bootstrap for sigmas) of MC numbers.  All parameters
come from ``scripts/_paperconfig.py`` (nm, seed 42); only the public API of ``donutloc`` is used.
Iterative-MINFLUX numbers are read from ``data/iterative_sweep.json``
(``scripts/run_iterative_sweep.py``); if that file is missing the sweep is run first.

``numbers.tex`` defines, for every key ``k``, the macro ``\\csname pnum@k\\endcsname`` and the
user commands ``\\pnum{k}`` (formatted value) and ``\\pnumse{k}`` (formatted SE, MC numbers).

Usage:  python scripts/compute_paper_numbers.py [--quick]
``--quick`` reduces the MC repetitions (development only) and writes to the separate files
``data/quick/paper_numbers.json`` (with a top-level ``"quick": true``) and
``paper/generated/quick/numbers.tex``; the final files are never touched by a quick run.  Shared points use the same parameters and seed as the figure scripts, so the figures
print the same values (``NUMBER <key> = ...``).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _paperconfig as C  # noqa: E402

SCRIPT = "scripts/compute_paper_numbers.py"
NUMBERS_PATH = os.path.join(C.DATA, "paper_numbers.json")
TEX_PATH = os.path.join(C.ROOT, "paper", "generated", "numbers.tex")
SWEEP_PATH = os.path.join(C.DATA, "iterative_sweep.json")
QUICK_NUMBERS_PATH = os.path.join(C.DATA, C.QUICK_DIRNAME, "paper_numbers.json")
QUICK_TEX_PATH = os.path.join(C.ROOT, "paper", "generated", C.QUICK_DIRNAME, "numbers.tex")


# ---------------------------------------------------------------------------- tex formatting
# Number formatting rules for paper/generated/numbers.tex (tested in tests/test_paper_tooling.py):
#
#   * ints (and bools) print as ints: ``3``.
#   * a value WITHOUT a standard error prints with ``SIG_FIGS`` = 4 significant figures and keeps
#     its trailing zeros (``1.600``, ``1.000``, ``0.07500``); 1e4 <= |v| < 1e5 prints as an
#     integer (``12346``); an exactly integral float with 100 <= |v| < 1e5 (a set parameter such
#     as fwhm = 300.0 nm) prints as an integer (``300``, not ``300.0``); |v| < 1e-3 or
#     |v| >= 1e5 prints in scientific notation (``7.040\times10^{-5}``).
#   * a value WITH a standard error: the SE is rounded to 2 significant figures when its two
#     leading digits are 10-24, otherwise to 1 (``0.0101 -> 0.010``, ``0.005131 -> 0.005``), and
#     the value is rounded to the same decimal place (``0.99123 +- 0.005131 -> 0.991 +- 0.005``).
#     ``\pnum{k}`` and ``\pnumse{k}`` of such a key both use this pair, so a value is never quoted
#     with more digits than its SE supports.  If |v| is in the scientific range, value and SE
#     share the value's power of ten.  A zero or non-finite SE falls back to the value-only rule.
SIG_FIGS = 4


def _sci_range(a):
    return a != 0.0 and (a < 1e-3 or a >= 1e5)


def _sci(mant, exp):
    return r"\ensuremath{%s\times10^{%d}}" % (mant, exp)


def _fixed(v, decimals):
    """``v`` rounded to ``decimals`` decimal places (negative: tens, hundreds, ...), as text."""
    if decimals <= 0:
        return "%d" % int(round(v, decimals))
    return "%.*f" % (decimals, v)


def fmt_value(v, sig=SIG_FIGS):
    """Format a value without SE (see the rules above)."""
    if isinstance(v, (bool, np.bool_)):
        return "%d" % int(v)
    if isinstance(v, (int, np.integer)):
        return "%d" % int(v)
    v = float(v)
    if not np.isfinite(v):
        return r"\ensuremath{\infty}" if v > 0 else r"\ensuremath{-\infty}"
    if v == 0.0:
        return "0"
    a = abs(v)
    if v.is_integer() and 100 <= a < 1e5:
        return "%d" % int(v)          # exact parameters such as fwhm = 300.0 nm print as '300'
    if _sci_range(a):
        m, e = ("%.*e" % (sig - 1, v)).split("e")
        return _sci(m, int(e))
    s = "%#.*g" % (sig, v)           # '#' keeps trailing zeros
    if "e" in s:                      # 1e4 <= |v| < 1e5 with sig=4
        return "%.0f" % v
    return s.rstrip(".")               # '1000.' -> '1000'


def se_decimals(se):
    """Decimal place to which an SE is rounded: 2 significant figures if its two leading digits
    are 10-24, else 1.  Returns the number of decimals (negative for tens, hundreds, ...)."""
    se = abs(float(se))
    e = int(np.floor(np.log10(se)))
    lead = se / 10.0 ** (e - 1)       # in [10, 100)
    nsig = 2 if round(lead, 6) < 25 else 1
    return nsig - 1 - e


def fmt_value_se(v, se):
    """Return ``(value_text, se_text)`` for a value with a standard error (rules above)."""
    if isinstance(v, (bool, np.bool_, int, np.integer)) or se is None:
        return fmt_value(v), (None if se is None else fmt_value(se))
    v, se = float(v), abs(float(se))
    if not (np.isfinite(v) and np.isfinite(se)) or se == 0.0:
        return fmt_value(v), fmt_value(se)
    if _sci_range(abs(v)):
        exp = int(np.floor(np.log10(abs(v))))
        scale = 10.0 ** exp
        d = max(se_decimals(se / scale), 0)
        return _sci("%.*f" % (d, v / scale), exp), _sci("%.*f" % (d, se / scale), exp)
    d = se_decimals(se)
    return _fixed(v, d), _fixed(se, d)


_fmt = fmt_value   # backwards-compatible name


def numbers_tex(numbers):
    """Text of paper/generated/numbers.tex for the dict ``numbers`` (deterministic)."""
    lines = [
        "% Generated by scripts/compute_paper_numbers.py from data/paper_numbers.json.",
        "% DO NOT EDIT: rerun the script. Use \\pnum{key} (value) and \\pnumse{key} (MC SE).",
        "\\providecommand{\\pnum}[1]{\\ifcsname pnum@#1\\endcsname\\csname pnum@#1\\endcsname"
        "\\else\\textbf{??\\detokenize{#1}??}\\fi}",
        "\\providecommand{\\pnumse}[1]{\\ifcsname pnumse@#1\\endcsname\\csname pnumse@#1\\endcsname"
        "\\else\\textbf{??se:\\detokenize{#1}??}\\fi}",
    ]
    for k in sorted(numbers):
        e = numbers[k]
        if not isinstance(e, dict):          # metadata such as "quick": true
            continue
        se = e.get("se")
        val, set_ = fmt_value_se(e["value"], se)
        lines.append("\\expandafter\\def\\csname pnum@%s\\endcsname{%s}" % (k, val))
        if se is not None:
            lines.append("\\expandafter\\def\\csname pnumse@%s\\endcsname{%s}" % (k, set_))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------- registry
class Registry(dict):
    def add(self, key, value, unit, description, se=None, script=SCRIPT):
        if key in self:
            raise KeyError("duplicate key %s" % key)
        if isinstance(value, (np.integer,)):
            value = int(value)
        elif isinstance(value, (bool, np.bool_)):
            value = int(bool(value))
        elif not isinstance(value, int):
            value = float(value)
        e = {"value": value, "unit": unit, "script": script, "description": description}
        if se is not None:
            e["se"] = float(se)
        self[key] = e
        print("NUMBER %s = %s%s %s" % (key, value, "" if se is None else " +- %.3g" % se, unit),
              flush=True)


def _tag(v):
    return ("%g" % v).replace(".", "p")


def _slope(x, y):
    return float(np.polyfit(np.log(np.asarray(x, float)), np.log(np.asarray(y, float)), 1)[0])


# ---------------------------------------------------------------------------- sections
def crb_section(R):
    from donutloc import beams, closed_forms, fisher, patterns, photons
    lg = beams.make_beam("donut", fwhm=C.FWHM)

    def model(L, sbr=None, fwhm=C.FWHM):
        b = lg if fwhm == C.FWHM else beams.make_beam("donut", fwhm=fwhm)
        return photons.make_model(patterns.tcp_centers(L), b, sbr=sbr)

    L, N = C.L_REF, C.N_REF
    R.add("fwhm_nm", C.FWHM, "nm", "LG donut size parameter fwhm (Balzarotti Eq. S17); ring "
          "peak-to-peak diameter fwhm/sqrt(ln2)")
    R.add("fwhm_masullo_nm", C.FWHM_MASULLO, "nm", "fwhm that reproduces the Masullo table")
    R.add("lg_ring_diameter_nm", 2 * beams.ring_radius(C.FWHM), "nm",
          "LG ring peak-to-peak diameter fwhm/sqrt(ln2) at fwhm=300 (also where the eps=0 centre "
          "CRB diverges)")
    lim = fisher.crb_limit(model(L), N)
    R.add("crb_center_lg_L50_N100_nm", lim, "nm", "centre CRB of the TCP, r->0 limit (acceptance "
          "definition, fisher.crb_limit, 12 directions at 1e-3 nm), LG fwhm=300, L=50, N=100, no "
          "background")
    R.add("crb_center_lg_L50_N100_closed_nm", closed_forms.crb_tcp_center_limit(L, N, C.FWHM),
          "nm", "same, closed form sigma_lim^2 = L^2/(8Ng^2)(3g^2+e^x)/(3g^2+2e^x)")
    R.add("crb_center_point_S27_L50_N100_nm", closed_forms.crb_tcp_center_point(L, N, C.FWHM),
          "nm", "point value at r=0 exactly (Balzarotti Eq. S27, central p=0 term excluded), L=50, "
          "N=100")
    for Lr in (5.0, 50.0, 100.0, 150.0):
        R.add("limit_to_point_ratio_L%d" % Lr, closed_forms.limit_to_point_ratio(Lr, C.FWHM), "",
              "sigma_lim/sigma_S27 at fwhm=300, L=%g nm" % Lr)
    R.add("limit_to_point_ratio_quadratic", 2.0 / np.sqrt(5.0), "",
          "sigma_lim/sigma_S27 for a quadratic zero (fwhm->inf): 2/sqrt5 exactly")
    R.add("limit_to_point_ratio_quadratic_numeric",
          fisher.crb_limit(model(L, fwhm=1e9), N) / closed_forms.crb_tcp_center_point(L, N),
          "", "numeric crb_limit / S27 with fwhm=1e9 nm (quadratic case check)")
    # scaling exponents
    cL = [fisher.crb_limit(model(l), N) for l in C.L_FIT]
    R.add("crb_exponent_L", _slope(C.L_FIT, cL), "", "log-log slope of the centre CRB (limit) vs "
          "L over L_FIT=%s nm (L<<fwhm), N=100, fwhm=300" % (list(C.L_FIT),))
    cN = fisher.crb(model(L), [0.0, 0.0], np.array(C.N_FIT, float), zero_policy="limit")
    R.add("crb_exponent_N", _slope(C.N_FIT, cN), "", "log-log slope of the centre CRB (limit) vs "
          "N over N_FIT=%s, L=50" % (list(C.N_FIT),))
    # limit ellipse (V4)
    for (Lr, fw, tag) in ((100.0, np.inf, "L100_N100_quadratic"), (50.0, C.FWHM, "L50_N100")):
        par, perp = closed_forms.crb_tcp_center_limit_axes(Lr, N, fw)
        R.add("crb_limit_axis_par_%s_nm" % tag, par, "nm", "r->0 limit error ellipse: sigma along "
              "the approach direction (V4), %s" % tag)
        R.add("crb_limit_axis_perp_%s_nm" % tag, perp, "nm", "r->0 limit error ellipse: sigma "
              "perpendicular to the approach direction (= S27 value), %s" % tag)
    R.add("crb_limit_axis_ratio_quadratic", np.sqrt(3.0 / 5.0), "",
          "sigma_par/sigma_perp of the limit ellipse, quadratic zero: sqrt(3/5)")
    # background (S31)
    R.add("crb_center_sbr10_L50_N100_nm", closed_forms.crb_tcp_center_point(L, N, C.FWHM, 10.0),
          "nm", "centre CRB with SBR=10 (Balzarotti Eq. S31; point value = limit), L=50, N=100")
    R.add("crb_center_sbr10_limit_L50_N100_nm", fisher.crb_limit(model(L, 10.0), N), "nm",
          "same, numeric r->0 limit (continuous with background)")
    # photons for 5 nm (W9)
    for s in (np.inf, 50.0, 20.0, 10.0, 5.0):
        c = closed_forms.crb_tcp_center_point(L, 100, C.FWHM, s)
        R.add("photons_for_5nm_L50_sbr%s" % ("inf" if np.isinf(s) else "%d" % s),
              100.0 * (c / 5.0) ** 2, "photons", "photons for a 5 nm centre CRB (point value, Eq. "
              "S31/S27), L=50, fwhm=300, SBR=%s" % s)
    R.add("photons_for_5nm_L50_limit_nobg", 100.0 * (lim / 5.0) ** 2, "photons",
          "photons for 5 nm with the r->0 limit, no background, L=50")
    # Masullo table (V12)
    for fw in (C.FWHM_MASULLO, C.FWHM):
        for Lr in C.L_LIST:
            R.add("crb_center_sbr5_N500_fwhm%d_L%d_nm" % (fw, Lr),
                  closed_forms.crb_tcp_center_point(Lr, 500, fw, 5.0), "nm",
                  "Eq. S31 centre CRB, N=500, SBR=5, fwhm=%g, L=%g (Masullo table check)" % (fw, Lr))
    R.add("crb_1d_quadratic_L50_N100_nm", closed_forms.crb_1d_center(L, N, kind="quadratic"), "nm",
          "1D two-exposure CRB L/(4 sqrt N) (Eq. S22c)")
    for c in (1, 2, 3):
        R.add("crb_center_multiphoton_c%d_L100_N100_quadratic_nm" % c,
              closed_forms.crb_tcp_center_point(100.0, N, np.inf, power=c), "nm",
              "centre CRB with multiphoton exponent c=%d (sigma ~ 1/c), L=100, quadratic" % c)


def camera_section(R):
    from donutloc import camera
    R.add("camera_sigma_nm", C.SIGMA_PSF / np.sqrt(C.ITER_N_TOTAL), "nm",
          "IDEAL camera sigma_PSF/sqrt(N) with sigma_PSF=100 nm, N=1000 (same photons as "
          "iterative_sigma_nm; declared ideal, no pixelation/background)")
    R.add("camera_ideal_N400_nm", C.SIGMA_PSF / np.sqrt(400.0), "nm",
          "ideal camera sigma_PSF/sqrt(N), N=400 (Balzarotti p. 1: 5 nm)")
    R.add("camera_pixelated_9x9_N400_nm", camera.crb_camera(C.SIGMA_PSF, 400, C.CAM_PIXEL,
          C.CAM_NPIX), "nm", "pixelated camera CRB, 9x9 pixels of 100 nm, sigma_PSF=100, N=400, "
          "no background")
    R.add("camera_pixelated_9x9_N1000_nm", camera.crb_camera(C.SIGMA_PSF, C.ITER_N_TOTAL,
          C.CAM_PIXEL, C.CAM_NPIX), "nm", "same, N=1000")
    c600 = camera.crb_camera(C.SIGMA_PSF, 600, C.CAM_PIXEL, C.CAM_NPIX, sbr=500.0,
                             sbr_convention="per_pixel")
    R.add("camera_perpixel_sbr500_N600_nm", c600, "nm", "pixelated 9x9 camera, SBR_c=500 defined "
          "per pixel (total signal / background per pixel; SBR_total=500/81), N=600")
    R.add("camera_photons_for_5nm_perpixel_sbr500", 600.0 * (c600 / 5.0) ** 2, "photons",
          "photons for 5 nm, same convention")
    t600 = camera.crb_camera(C.SIGMA_PSF, 600, C.CAM_PIXEL, C.CAM_NPIX, sbr=500.0,
                             sbr_convention="total")
    R.add("camera_total_sbr500_N600_nm", t600, "nm", "same camera if SBR=500 were the total SBR")
    R.add("camera_photons_for_5nm_total_sbr500", 600.0 * (t600 / 5.0) ** 2, "photons",
          "photons for 5 nm if SBR=500 were total")
    R.add("camera_sbr_perpixel_to_total", 500.0 / C.CAM_NPIX ** 2, "",
          "SBR_total equivalent to SBR_c=500 per pixel on 9x9 (500/81)")


def vectorial_section(R):
    from donutloc import beams, fisher, patterns, photons, vectorial
    VO = dict(wavelength=C.WAVELENGTH, NA=C.NA, n=C.N_MEDIUM, filling=C.FILLING)
    R.add("vectorial_zero_depth_correct", vectorial.zero_depth(**VO), "",
          "I(0)/I_max of the Richards-Wolf vortex donut, circular pol. correct hand (l s=+1), "
          "lambda=640, NA=1.4, n=1.518, F=5/3 (exact zero; machine precision)")
    R.add("vectorial_zero_depth_wrong_handedness", vectorial.zero_depth(handedness=-1, **VO), "",
          "same, opposite handedness (on-axis Ez != 0)")
    R.add("vectorial_zero_depth_linear", vectorial.zero_depth(polarization="linear", **VO), "",
          "same, linear polarization (I_max searched in 2D)")
    R.add("vectorial_peak_to_peak_diameter_nm", vectorial.peak_to_peak_diameter(**VO), "nm",
          "ring peak-to-peak diameter, correct hand")
    R.add("vectorial_zero_curvature_nm2", vectorial.zero_curvature(**VO), "nm^-2",
          "curvature c of I/I_max ~ c rho^2 at the zero, correct hand")
    fc = vectorial.lg_equivalent_fwhm(VO, match="curvature")
    R.add("vectorial_lg_fwhm_curvature_nm", fc, "nm", "LG fwhm with the same zero curvature")
    R.add("vectorial_lg_fwhm_diameter_nm", vectorial.lg_equivalent_fwhm(VO, match="diameter"),
          "nm", "LG fwhm with the same peak-to-peak diameter")
    cases = (("correct", {}), ("opposite", dict(handedness=-1)),
             ("linear", dict(polarization="linear", pol_angle=0.0)))
    vb = {name: vectorial.make_vectorial_beam(mode="exact", **dict(VO, **o)) for name, o in cases}
    lgs = {"lg300": beams.make_beam("donut", fwhm=C.FWHM), "lgcurv": beams.make_beam("donut",
                                                                                     fwhm=fc)}
    for Lr in C.L_LIST:
        tag = "L%d_N%d" % (Lr, C.N_REF)
        cen = patterns.tcp_centers(Lr)
        val = {}
        for name, b in list(vb.items()) + list(lgs.items()):
            val[name] = fisher.crb_limit(photons.make_model(cen, b), C.N_REF)
        for name in ("correct", "opposite", "linear"):
            R.add("crb_center_vec_%s_%s_nm" % (name, tag), val[name], "nm",
                  "centre CRB (r->0 limit), vectorial donut %s (mode='exact'), %s" % (name, tag))
        R.add("crb_center_lg300_%s_nm" % tag, val["lg300"], "nm", "centre CRB, LG fwhm=300, %s"
              % tag)
        R.add("crb_center_lgcurv_%s_nm" % tag, val["lgcurv"], "nm",
              "centre CRB, curvature-matched LG (fwhm=%.2f), %s" % (fc, tag))
        R.add("crb_vec_over_lgcurv_%s" % tag, val["correct"] / val["lgcurv"], "",
              "vectorial (correct) / curvature-matched LG centre CRB, %s" % tag)
        R.add("crb_vec_over_lg300_%s" % tag, val["correct"] / val["lg300"], "",
              "vectorial (correct) / LG300 centre CRB, %s" % tag)


def estimator_section(R, quick):
    from donutloc import beams, closed_forms, estimators, fisher, montecarlo, patterns, photons
    L, N = C.L_REF, C.N_REF
    beam = beams.make_beam("donut", fwhm=C.FWHM)
    p_bg = photons.make_model(patterns.tcp_centers(L), beam, sbr=C.SBR_MLE)
    p_0 = photons.make_model(patterns.tcp_centers(L), beam)
    reps = 1000 if quick else C.N_REP_MLE
    rad = C.MLE_RADIUS_CENTRE_OVER_L * L       # MLE search disk radius (r1 convention, fig 5)

    def sig(mc, r_true):
        e = mc["estimates"] - np.asarray(r_true, float)
        return (montecarlo.sigma_of_errors(e),
                montecarlo.bootstrap_sigma_se(e, n_boot=C.N_BOOT, seed=C.SEED))

    mc = montecarlo.run_mc(lambda c: estimators.mle(c, p_bg, search_radius=rad), p_bg,
                           np.zeros(2), N, reps, seed=C.SEED, n_boot=0)
    s, se = sig(mc, [0, 0])
    crb10 = fisher.crb_limit(p_bg, N)
    common = "L=50, N=100, fwhm=300, MLE in disk radius L, %d reps, seed 42, bootstrap SE" % reps
    R.add("mle_sigma_center_sbr10_nm", s, "nm", "MLE sigma at the TCP centre, SBR=10; " + common,
          se=se)
    R.add("mle_efficiency_center", s / crb10, "", "sigma_MLE/CRB at the centre with SBR=10 "
          "(regular model); " + common, se=se / crb10)
    mc0 = montecarlo.run_mc(lambda c: estimators.mle(c, p_0, search_radius=rad), p_0,
                            np.zeros(2), N, reps, seed=C.SEED, n_boot=0)
    s0, se0 = sig(mc0, [0, 0])
    lim = fisher.crb_limit(p_0, N)
    s27 = closed_forms.crb_tcp_center_point(L, N, C.FWHM)
    R.add("mle_nobg_sigma_center_nm", s0, "nm", "MLE sigma at the centre WITHOUT background; "
          + common, se=se0)
    R.add("mle_nobg_sigma_over_crb_limit_center", s0 / lim, "", "superefficiency: sigma_MLE / "
          "CRB r->0 limit at the centre, no background (non-regular model); " + common,
          se=se0 / lim)
    R.add("mle_nobg_sigma_over_s27_center", s0 / s27, "", "sigma_MLE / S27 point value, centre, "
          "no background; " + common, se=se0 / s27)
    rb = C.N_REP_BIAS0 // 10 if quick else C.N_REP_BIAS0
    r2 = np.array([2.0, 0.0])
    mcb = montecarlo.run_mc(lambda c: estimators.mle(c, p_0, search_radius=rad), p_0, r2, N, rb,
                            seed=C.SEED, n_boot=0)
    e = mcb["estimates"] - r2
    R.add("mle_nobg_bias_x_r2_nm", float(e[:, 0].mean()), "nm", "MLE bias along x at r=(2,0), no "
          "background (toward the centre); L=50, N=100, %d reps, seed 42" % rb,
          se=float(e[:, 0].std(ddof=1) / np.sqrt(rb)))
    # off-centre MLE with SBR=10 (fig 5 x-sweep point x0=15: search disk radius 2L)
    r15 = np.array([15.0, 0.0])
    for n15 in (N, 1000):
        mc15 = montecarlo.run_mc(lambda c: estimators.mle(c, p_bg, search_radius=C.MLE_RADIUS_SWEEP_OVER_L * L), p_bg,
                                 r15, n15, reps, seed=C.SEED, n_boot=0)
        s15, se15 = sig(mc15, r15)
        crb15 = float(fisher.crb(p_bg, r15, n15))
        R.add("mle_sbr10_sigma_over_crb_x15_N%d" % n15, s15 / crb15, "", "sigma_MLE/CRB at "
              "r=(15,0) nm with SBR=10, N=%d, L=50, MLE disk radius 2L (fig 5 sweep), %d reps, "
              "seed 42, bootstrap SE (outlier tail at low N)" % (n15, reps), se=se15 / crb15)
    mcl = montecarlo.run_mc(lambda c: estimators.lms_tcp(c, L, C.FWHM), p_0, np.zeros(2), N, reps,
                            seed=C.SEED, n_boot=0)
    sl, sel = sig(mcl, [0, 0])
    R.add("lms_nobg_sigma_over_s27_center", sl / s27, "", "LMS (Eq. S49-S50) sigma / S27 at the "
          "centre, no background (analytically 1); L=50, N=100, %d reps, seed 42" % reps,
          se=sel / s27)
    R.add("lms_nobg_sigma_over_crb_limit_center", sl / lim, "", "LMS sigma / CRB limit at the "
          "centre (analytically 1/rho = 1.123)", se=sel / lim)
    # LMS with background but without the 1/s factor: noise-free shrink factor (= s = 10/11)
    r5 = np.array([5.0, 0.0])
    cnt = N * np.asarray(p_bg(r5), float)
    R.add("lms_sbr10_shrink_factor", estimators.lms_tcp(cnt, L, C.FWHM)[0]
          / estimators.lms_tcp(cnt, L, C.FWHM, sbr=C.SBR_MLE)[0], "",
          "LMS on the expected counts at r=(5,0) with SBR=10: estimate without / with the 1/s "
          "factor (= s = SBR/(SBR+1) = 10/11 at any r), L=50, fwhm=300")


def eps_section(R):
    from donutloc import closed_forms, experiments
    for Lr in (50.0, 100.0):
        s27 = closed_forms.crb_tcp_center_point(Lr, C.N_REF, C.FWHM)
        for e in C.EPS_LIST:
            for zm in ("constant", "gaussian"):
                R.add("crb_eps_degradation_L%d_eps%s_%s" % (Lr, _tag(e), zm),
                      closed_forms.crb_tcp_center_eps(Lr, C.N_REF, C.FWHM, e, zero_model=zm) / s27,
                      "", "centre CRB(eps)/CRB_S27, %s pedestal, L=%g, fwhm=300, no background "
                      "(V9)" % (zm, Lr))
    for e in C.EPS_LIST:
        o = experiments.optimal_L(e, N=C.N_REF, fwhm=C.FWHM, zero_model="gaussian")
        t = _tag(e)
        R.add("L_opt_eps%s_nm" % t, o["L_opt"], "nm", "L minimizing the centre CRB, gaussian "
              "pedestal eps=%g, N=100, fwhm=300 (experiments.optimal_L)" % e)
        R.add("crb_opt_eps%s_nm" % t, o["crb_opt"], "nm", "centre CRB at L_opt, eps=%g, N=100" % e)
        R.add("L_opt_over_fwhm_sqrt_eps_eps%s" % t, o["L_opt"] / (C.FWHM * np.sqrt(e)), "",
              "L_opt/(fwhm sqrt eps), eps=%g (0.78 only for eps<~0.01)" % e)
    e5, L5 = 0.05, 100.0
    se5 = closed_forms.sbr_eff_constant_pedestal(L5, C.FWHM, e5)
    s31 = closed_forms.crb_tcp_center_point(L5, C.N_REF, C.FWHM, sbr=se5)
    R.add("eps_constant_pedestal_sbr_eff_L100_eps0p05", se5, "", "effective SBR of a constant "
          "pedestal eps=0.05 at the centre, SBR_eps = 3 e x e^-x/(4 eps), L=100, fwhm=300")
    for zm in ("constant", "gaussian"):
        R.add("eps_%s_pedestal_over_s31_L100_eps0p05" % zm,
              closed_forms.crb_tcp_center_eps(L5, C.N_REF, C.FWHM, e5, zero_model=zm) / s31, "",
              "centre CRB with a %s pedestal eps=0.05 / Eq. S31 at SBR=SBR_eps (1 exactly for the "
              "constant pedestal; not for the gaussian one), L=100, N=100, fwhm=300" % zm)
    a = 4.0 * np.log(2.0) / C.FWHM ** 2
    R.add("zero_depth_transition_scale_nm", np.sqrt(C.EPS_TRANSITION / (np.e * a)), "nm",
          "sqrt(eps/(e a)) for eps=%g: radial scale where the pedestal and the LG zero are "
          "comparable (V10; a scale, not the argmin)" % C.EPS_TRANSITION)
    R.add("eps0_crb_divergence_L_nm", C.FWHM / np.sqrt(np.log(2.0)), "nm",
          "with eps=0 the centre CRB diverges at L = fwhm/sqrt(ln2) (ring diameter); monotonic "
          "only for smaller L")


def misalignment_section(R, quick):
    from donutloc import experiments
    mis = dict(C.MIS)
    if quick:
        mis.update(n_patterns=10, n_rep=100)
    m = experiments.misalignment_study(C.MIS_DELTAS, seed=C.SEED, **mis)
    P, Rr = mis["n_patterns"], mis["n_rep"]
    base = "L=%g, N=%d, SBR=%g, %d patterns x %d reps, seed 42" % (mis["L"], mis["N"], mis["sbr"],
                                                                   P, Rr)
    tags = ["center", "Lq"]
    for i, d in enumerate(m["delta"]):
        dt = _tag(d)
        for q, tag in enumerate(tags):
            where = "at the centre" if q == 0 else "at (L/4, 0)"
            for est in ("naive", "honest"):
                r = m[est]
                R.add("misalignment_%s_bias_abs_%s_d%s_nm" % (est, tag, dt),
                      r["bias_abs_mean"][i, q], "nm", "mean |bias| of the %s MLE %s, delta=%g; %s"
                      % (est, where, d, base), se=r["bias_abs_se"][i, q])
                R.add("misalignment_%s_sigma_%s_d%s_nm" % (est, tag, dt), r["sigma_mean"][i, q],
                      "nm", "mean sigma of the %s MLE %s, delta=%g (SE between patterns); %s"
                      % (est, where, d, base), se=r["sigma_se"][i, q])
                R.add("misalignment_%s_rmse_%s_d%s_nm" % (est, tag, dt), r["rmse"][i, q], "nm",
                      "rmse of the %s MLE %s, delta=%g; %s" % (est, where, d, base))
            R.add("misalignment_crb_honest_%s_d%s_nm" % (tag, dt), m["crb_honest"][i, q], "nm",
                  "mean honest CRB %s, delta=%g" % (where, d))
            R.add("misalignment_mc_floor_%s_d%s_nm" % (tag, dt), m["bias_noise_floor"][i, q],
                  "nm", "MC noise floor of |bias| sigma sqrt(pi/(2R)) %s, delta=%g" % (where, d))
    R.add("misalignment_delta0_identical", int(m["identical"][0]), "",
          "1 if naive and honest estimates coincide bit for bit at delta=0")
    dl = np.asarray(m["delta"])
    for q, tag in enumerate(tags):
        b = m["naive"]["bias_abs_mean"][:, q]
        se = m["naive"]["bias_abs_se"][:, q]
        k = dl >= 2.0
        w = 1.0 / se[k] ** 2
        s = np.sum(w * dl[k] * b[k]) / np.sum(w * dl[k] ** 2)
        R.add("misalignment_naive_bias_over_delta_%s" % tag, s, "", "weighted LS slope of the "
              "naive |bias| vs delta through the origin, delta in %s (paper numbers; fig 8 fits "
              "a denser delta grid under misalignment_naive_bias_slope_*); %s"
              % ([float(v) for v in dl[k]], base), se=float(1.0 / np.sqrt(np.sum(w * dl[k] ** 2))))


def misalignment_population_section(R, quick):
    """Noise-free (population) naive-MLE bias under misalignment (inbox r4): MLE on expected
    counts, mean over MIS_POP_N_PATTERNS random patterns (seed 42), centre and (L/4, 0)."""
    from donutloc import experiments
    P = 200 if quick else C.MIS_POP_N_PATTERNS
    mis = C.MIS
    m = experiments.misalignment_population_bias(C.MIS_POP_DELTAS, L=mis["L"], N=mis["N"],
                                                 sbr=mis["sbr"], n_patterns=P, seed=C.SEED)
    base = ("noise-free: naive MLE (ideal-TCP model, disk radius 0.75 L) on the EXPECTED counts "
            "of each perturbed TCP (no Poisson noise, no MC floor), mean over %d random "
            "misalignment patterns (4 zeros each moved by delta in a random direction, seed 42); "
            "L=%g, SBR=%g; SE = pattern sampling only" % (P, mis["L"], mis["sbr"]))
    R.add("misalignment_pop_n_patterns", int(P), "", "number of random misalignment patterns "
          "of the noise-free population estimate (misalignment_naive_bias_over_delta_pop_*)")
    for q, (tag, where) in enumerate((("center", "at the centre"), ("Lq", "at (L/4, 0)"))):
        R.add("misalignment_naive_bias_over_delta_pop_%s" % tag, m["slope"][q], "",
              "population slope of the naive |bias| vs delta %s: unweighted least-squares fit "
              "through the origin sum(delta b)/sum(delta^2) over delta=%s nm; %s"
              % (where, [float(d) for d in C.MIS_POP_DELTAS], base), se=m["slope_se"][q])
        for i, d in enumerate(m["delta"]):
            R.add("misalignment_naive_bias_over_delta_pop_%s_d%s" % (tag, _tag(d)),
                  m["ratio"][i, q], "", "population mean naive |bias| / delta %s, delta=%g nm "
                  "(grows with delta: nonlinear); %s" % (where, d, base), se=m["ratio_se"][i, q])


def iterative_section(R, quick):
    if quick and not os.path.exists(SWEEP_PATH):
        import run_iterative_sweep
        d = run_iterative_sweep.run_sweep(quick=True)
    elif not os.path.exists(SWEEP_PATH):
        import run_iterative_sweep
        d = run_iterative_sweep.run_sweep(quick=False)
        with open(SWEEP_PATH, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(d, fh, indent=1, sort_keys=True)
            fh.write("\n")
    else:
        with open(SWEEP_PATH, encoding="utf-8") as fh:
            d = json.load(fh)
    src = "scripts/run_iterative_sweep.py"
    p = d["params"]
    base = ("iterative MINFLUX, 4 iterations L=150->25 nm, equal split, re-centring, MLE, no "
            "background, %d reps, seed 42 (data/iterative_sweep.json)" % p["n_rep"])
    j = d["N"].index(C.ITER_N_TOTAL)
    R.add("iterative_sigma_nm", d["sigma"][j], "nm", "final sigma at N_total=1000; " + base,
          se=d["sigma_se"][j], script=src)
    R.add("iterative_ratio_to_camera_N1000", d["ratio_to_camera"][j], "",
          "iterative sigma / ideal camera sigma_PSF/sqrt(N), N=1000", se=d["per_N"][j]["ratio_se"],
          script=src)
    for i, n in enumerate(d["N"]):
        R.add("iterative_sigma_N%d_nm" % n, d["sigma"][i], "nm", "final sigma at N_total=%d; %s"
              % (n, base), se=d["sigma_se"][i], script=src)
    for s in ("slope_all", "slope_N_ge_500"):
        e = d[s]
        R.add("iterative_%s" % s, e["value"], "", "log-log slope of sigma vs N_total over N=%s "
              "(bootstrap SE)" % e["N_used"], se=e["se"], script=src)
        R.add("iterative_%s_ci95_lo" % s, e["ci95"][0], "", "bootstrap 95%% CI low, %s" % s,
              script=src)
        R.add("iterative_%s_ci95_hi" % s, e["ci95"][1], "", "bootstrap 95%% CI high, %s" % s,
              script=src)
    R.add("iterative_ratio_min", d["ratio_range"][0], "", "min over N of sigma/camera", script=src)
    R.add("iterative_ratio_max", d["ratio_range"][1], "", "max over N of sigma/camera", script=src)
    ex = d["extras_N_ref"]
    R.add("iterative_sbr10_sigma_nm", ex["sbr10"]["sigma"], "nm", "same protocol with SBR=10, "
          "N=1000", se=ex["sbr10"]["sigma_se"], script=src)
    R.add("iterative_adaptive_sigma_nm", ex["adaptive"]["sigma"], "nm", "adaptive L rule "
          "(kappa=6, L=%s), N=1000" % [round(v, 2) for v in ex["adaptive"]["L"]],
          se=ex["adaptive"]["sigma_se"], script=src)
    R.add("iterative_no_recenter_sigma_nm_ARTEFACT", ex["no_recenter"]["sigma"], "nm",
          "recenter=False, N=1000: ARTEFACT of the MLE search disk (radius 0.75 L_k < emitter "
          "spread 37.5 nm), not physics", se=ex["no_recenter"]["sigma_se"], script=src)
    R.add("iterative_crb_all_photons_L25_N1000_nm", d["per_N"][j]["crb_all_photons_Lmin"], "nm",
          "centre CRB if all 1000 photons were spent at L=25", script=src)


def adaptive_coverage_section(R, quick):
    from donutloc import experiments
    reps = 2000 if quick else C.ITER_N_REP
    Nk = np.array([C.ITER_N_TOTAL // C.ITER["n_iter"]] * C.ITER["n_iter"])
    Ls = experiments.l_schedule(C.ITER["n_iter"], C.ITER["L0"], C.ITER["L_min"], rule="adaptive",
                                N_k=Nk)
    it0 = experiments.iterative_minflux(int(Nk[0]), L_schedule=[C.ITER["L0"]], n_rep=reps,
                                        seed=C.SEED)
    err = np.hypot(*(it0["estimates"] - it0["r_true"]).T)
    f = float(np.mean(err > Ls[1] / 2.0))
    R.add("adaptive_frac_outside_next_radius", f, "", "fraction of emitters whose iteration-0 "
          "error exceeds L_1/2=%.1f nm (adaptive rule, kappa=6, L0=150, 250 photons, %d reps, "
          "seed 42): coverage ~97%%, not 99%%" % (Ls[1] / 2, reps), se=np.sqrt(f * (1 - f) / reps))
    R.add("adaptive_iter0_sigma_nm", it0["sigma"], "nm", "real per-axis error of iteration 0 "
          "(L=150, 250 photons, emitters uniform in L0/4)", se=it0["sigma_se"])
    R.add("adaptive_iter0_crb_center_nm", it0["crb_center_iter"][0], "nm",
          "centre CRB of iteration 0 used by the adaptive rule (underestimates the real error)")


def physical_background_section(R, quick):
    """Fixed SBR (Eq. S30, our default) vs a physical background constant per exposure (Eq. S28)
    matched at the TCP centre (SBR_c = 10): fig 5 x-sweep points and the iterative protocol."""
    from donutloc import beams, estimators, experiments, fisher, montecarlo, patterns, photons
    L, N, sbr0 = C.L_REF, C.N_REF, C.SBR_MLE
    beam = beams.make_beam("donut", fwhm=C.FWHM)
    cen = patterns.tcp_centers(L)

    def bg_for(centers):
        # SBR_c = sum_i I_i(0) / (K b)  (Eq. S29 at the centre)
        return float(photons.intensities(np.zeros(2), centers, beam).sum()) / (centers.shape[0] * sbr0)

    bg = bg_for(cen)
    models = {"fixed": photons.make_model(cen, beam, sbr=sbr0),
              "phys": photons.make_model(cen, beam, bg_per_exposure=bg)}
    reps = 1000 if quick else C.N_REP_MLE
    rad = C.MLE_RADIUS_SWEEP_OVER_L * L
    base = ("L=50, N=100, fwhm=300, %d reps, seed 42, bootstrap SE; 'phys' = background constant "
            "per exposure with SBR=10 at the centre (Eq. S28), 'fixed' = SBR=10 at every position "
            "(Eq. S30, fig 5); MLE with the true model (disk radius 2L), LMS/mLMS with the 1/s "
            "factor for SBR=10" % reps)
    for x in C.BGPHYS_X:
        r = np.array([x, 0.0])
        xt = "x%d" % x
        R.add("bgphys_sbr_%s" % xt, float(photons.sbr_at(r, cen, beam, bg)), "",
              "SBR at r=(%g,0) with the physical background matched to SBR=10 at the centre "
              "(L=50, fwhm=300)" % x)
        for mk, p_fn in models.items():
            crb = float(fisher.crb(p_fn, r, N))
            R.add("bgphys_%s_crb_%s_nm" % (mk, xt), crb, "nm", "CRB at r=(%g,0); %s" % (x, base))
            fns = {"mle": lambda c, p=p_fn: estimators.mle(c, p, search_radius=rad),
                   "lms": lambda c: estimators.lms_tcp(c, L, C.FWHM, sbr=sbr0),
                   "mlms": lambda c: estimators.mlms_tcp(c, L, C.FWHM, sbr=sbr0)}
            for ek, fn in fns.items():
                mc = montecarlo.run_mc(fn, p_fn, r, N, reps, seed=C.SEED, n_boot=0)
                e = mc["estimates"] - r
                s = montecarlo.sigma_of_errors(e)
                se = montecarlo.bootstrap_sigma_se(e, n_boot=C.N_BOOT, seed=C.SEED)
                R.add("bgphys_%s_%s_bias_x_%s_nm" % (mk, ek, xt), float(e[:, 0].mean()), "nm",
                      "%s bias along x at r=(%g,0), %s background; %s" % (ek, x, mk, base),
                      se=float(e[:, 0].std(ddof=1) / np.sqrt(e.shape[0])))
                R.add("bgphys_%s_%s_sigma_over_crb_%s" % (mk, ek, xt), s / crb, "",
                      "%s sigma/CRB at r=(%g,0), %s background; %s" % (ek, x, mk, base),
                      se=se / crb)
    # iterative protocol (fig 6) with a physical background matched at L0 or at L_min
    ireps = 2000 if quick else C.ITER_N_REP
    for Lm in C.BGPHYS_ITER_MATCH_L:
        b_it = bg_for(patterns.tcp_centers(Lm))
        kw = dict(C.ITER)
        res = experiments.iterative_minflux(C.ITER_N_TOTAL, n_rep=ireps, seed=C.SEED,
                                            sigma_psf=C.SIGMA_PSF, bg_per_exposure=b_it, **kw)
        err = res["estimates"] - res["r_true"]
        se = montecarlo.bootstrap_sigma_se(err, n_boot=C.N_BOOT, seed=C.SEED)
        t = "L%d" % Lm
        ib = ("iterative protocol of fig 6 (4 iterations L=150->25, equal split, re-centring, "
              "MLE), N_tot=1000, %d reps, seed 42, with a background constant per exposure fixed "
              "so that SBR=10 at the centre of the L=%g pattern" % (ireps, Lm))
        R.add("bgphys_iter_match%s_sigma_nm" % t, float(res["sigma"]), "nm",
              "final sigma; " + ib, se=se)
        R.add("bgphys_iter_match%s_sbr_first" % t, float(res["sbr_center"][0]), "",
              "centre SBR of the first iteration (L=150); " + ib)
        R.add("bgphys_iter_match%s_sbr_last" % t, float(res["sbr_center"][-1]), "",
              "centre SBR of the last iteration (L=25); " + ib)


def misspecified_estimator_section(R):
    """Noise-free (population) bias of estimators that ignore a constant zero pedestal or
    misjudge the SBR: MLE / LMS applied to the EXPECTED counts N p_true(r) (as in
    misalignment_population_bias).  L=50, fwhm=300, true SBR=10 (fixed-SBR convention)."""
    from donutloc import beams, estimators, fisher, patterns, photons
    L, N, sbr0 = C.L_REF, C.N_REF, C.SBR_MLE
    cen = patterns.tcp_centers(L)
    rad = 0.75 * L
    b0 = beams.make_beam("donut", fwhm=C.FWHM)
    p_naive = photons.make_model(cen, b0, sbr=sbr0)
    base = ("noise-free: estimator applied to the expected counts N p_true(r) (no Poisson noise), "
            "L=50, fwhm=300, true SBR=10 at every position, MLE disk radius 0.75 L, LMS Eq. S50 "
            "with the 1/s factor of the assumed SBR")
    honest_max = 0.0
    for eps in C.NAIVE_EPS:
        bt = beams.make_beam("donut", fwhm=C.FWHM, eps=eps, zero_model="constant")
        p_true = photons.make_model(cen, bt, sbr=sbr0)
        et = _tag(eps)
        for x in C.NAIVE_X:
            r = np.array([x, 0.0])
            cnt = float(N) * np.asarray(p_true(r), float)
            cnt0 = float(N) * np.asarray(p_naive(r), float)
            xt = "x%d" % x
            b = float(np.hypot(*(estimators.mle(cnt, p_naive, search_radius=rad) - r)))
            crb = float(fisher.crb(p_true, r, N))
            R.add("naive_eps%s_mle_bias_abs_%s_nm" % (et, xt), b, "nm",
                  "|bias| (vector) of the MLE that ignores a constant pedestal eps=%g (model "
                  "eps=0) at r=(%g,0); %s" % (eps, x, base))
            R.add("naive_eps%s_crb_%s_nm" % (et, xt), crb, "nm",
                  "CRB (N=100) of the true model with constant pedestal eps=%g at r=(%g,0), "
                  "L=50, SBR=10" % (eps, x))
            R.add("naive_eps%s_mle_N_bias_eq_crb_%s" % (et, xt), float(N) * (crb / b) ** 2, "",
                  "photon number at which the CRB (~N^-1/2) falls to the N-independent |bias| "
                  "of the naive MLE, N (CRB_100/|bias|)^2, eps=%g, r=(%g,0)" % (eps, x))
            R.add("naive_eps%s_lms_extra_bias_x_%s_nm" % (et, xt),
                  float(estimators.lms_tcp(cnt, L, C.FWHM, sbr=sbr0)[0]
                        - estimators.lms_tcp(cnt0, L, C.FWHM, sbr=sbr0)[0]), "nm",
                  "change of the LMS estimate along x caused by an ignored constant pedestal "
                  "eps=%g at r=(%g,0) (LMS with eps minus LMS without eps; the linearization "
                  "bias of fig 5 is removed); %s" % (eps, x, base))
            h = estimators.mle(cnt, p_true, search_radius=rad)
            honest_max = max(honest_max, float(np.hypot(*(h - r))))
    for sa in C.NAIVE_SBR_ASSUMED:
        p_ass = photons.make_model(cen, b0, sbr=None if np.isinf(sa) else sa)
        st = "inf" if np.isinf(sa) else _tag(sa)
        for x in C.NAIVE_X:
            r = np.array([x, 0.0])
            cnt = float(N) * np.asarray(p_naive(r), float)
            R.add("naive_sbr%s_mle_bias_abs_x%d_nm" % (st, x),
                  float(np.hypot(*(estimators.mle(cnt, p_ass, search_radius=rad) - r))), "nm",
                  "|bias| (vector) of the MLE that assumes SBR=%g when the true SBR is 10 (eps=0) "
                  "at r=(%g,0); %s" % (sa, x, base))
    R.add("naive_honest_mle_max_abs_bias_nm", honest_max, "nm",
          "largest |bias| of the MLE with the correct pedestal model over the eps and x above "
          "(numerical zero of the noise-free test); " + base)


def simuflux_convention_section(R):
    """Centre CRB of the TCP with the SimuFLUX default donut (fwhm=310, L=75, N=100, no
    background): point value vs r->0 limit, to compare with their Fig. 2e (~2.8 nm)."""
    from donutloc import beams, closed_forms, fisher, patterns, photons
    L, fw, N = 75.0, 310.0, 100
    p0 = photons.make_model(patterns.tcp_centers(L), beams.make_beam("donut", fwhm=fw))
    R.add("simuflux_tcp_crb_point_L75_fwhm310_nm", closed_forms.crb_tcp_center_point(L, N, fw),
          "nm", "TCP centre CRB, point value Eq. S27 (centre exposure excluded), L=75, fwhm=310 "
          "(SimuFLUX default), N=100, no background")
    R.add("simuflux_tcp_crb_limit_L75_fwhm310_nm", float(fisher.crb_limit(p0, N)), "nm",
          "TCP centre CRB, r->0 limit, L=75, fwhm=310, N=100, no background")


def compute(quick=False):
    t0 = time.time()
    R = Registry()
    for name, fn in (("crb", lambda: crb_section(R)), ("camera", lambda: camera_section(R)),
                     ("eps", lambda: eps_section(R)),
                     ("iterative", lambda: iterative_section(R, quick)),
                     ("estimators", lambda: estimator_section(R, quick)),
                     ("adaptive", lambda: adaptive_coverage_section(R, quick)),
                     ("vectorial", lambda: vectorial_section(R)),
                     ("misalignment", lambda: misalignment_section(R, quick)),
                     ("misalignment_pop", lambda: misalignment_population_section(R, quick)),
                     ("bgphys", lambda: physical_background_section(R, quick)),
                     ("misspecified", lambda: misspecified_estimator_section(R)),
                     ("simuflux", lambda: simuflux_convention_section(R))):
        t = time.time()
        fn()
        print("# section %s done in %.1f s" % (name, time.time() - t), flush=True)
    print("# total %.1f s, %d keys" % (time.time() - t0, len(R)))
    return dict(R)


def write(numbers, numbers_path=NUMBERS_PATH, tex_path=TEX_PATH):
    os.makedirs(os.path.dirname(numbers_path), exist_ok=True)
    with open(numbers_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(numbers, fh, indent=1, sort_keys=True)
        fh.write("\n")
    os.makedirs(os.path.dirname(tex_path), exist_ok=True)
    with open(tex_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(numbers_tex(numbers))
    print("wrote %s and %s" % (numbers_path, tex_path))


def main(argv=None):
    ap = argparse.ArgumentParser(description="compute data/paper_numbers.json")
    ap.add_argument("--quick", action="store_true", help="reduced MC (development only); "
                    "writes data/quick/paper_numbers.json and paper/generated/quick/numbers.tex "
                    "(with 'quick': true) and never touches the final files")
    a = ap.parse_args(argv)
    numbers = compute(quick=a.quick)
    if a.quick:
        numbers["quick"] = True
        write(numbers, QUICK_NUMBERS_PATH, QUICK_TEX_PATH)
    else:
        write(numbers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
