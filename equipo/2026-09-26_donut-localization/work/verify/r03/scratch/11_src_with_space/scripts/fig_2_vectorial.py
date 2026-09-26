# -*- coding: utf-8 -*-
"""Figure 2 -- vectorial (Richards-Wolf) donut vs scalar LG donut.

(a-c) Focal-plane intensity (normalized to the ring maximum) of the vortex donut with circular
      polarization of the correct hand, the opposite hand, and linear (x) polarization.
(d)   Profiles through the zero and zero depths I(0)/I_max.
(e)   Centre CRB (r -> 0 limit, N = 100, no background) of the TCP vs L for the vectorial
      beams and the LG donut (fwhm 300 nm and curvature-matched 327.14 nm).
(f)   Relative difference of the correct-hand vectorial and LG 327 CRBs with respect to LG 300.

All vectorial beams are evaluated exactly (make_vectorial_beam(mode="exact")); no Cartesian
interpolation grid is used for the linear case.

Usage: python scripts/fig_2_vectorial.py [--quick] [--no-cache]
"""
import numpy as np

import _paperconfig as C
import _paperstyle as S
from donutloc import beams, fisher, patterns, photons, vectorial

VOPT = dict(wavelength=C.WAVELENGTH, NA=C.NA, n=C.N_MEDIUM, filling=C.FILLING)
CASES = (("correct", dict(polarization="circular", handedness=+1)),
         ("opposite", dict(polarization="circular", handedness=-1)),
         ("linear", dict(polarization="linear", pol_angle=0.0)))


def compute(quick):
    n_map = 61 if quick else 121
    g = np.linspace(-500.0, 500.0, n_map)
    X, Y = np.meshgrid(g, g)
    xp = np.linspace(-500.0, 500.0, 201 if quick else 401)
    L_arr = np.array([10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 125.0, 150.0, 200.0, 250.0])
    out = dict(g=g, xp=xp, L=L_arr)
    f_curv = vectorial.lg_equivalent_fwhm(dict(VOPT), match="curvature")
    out["f_curv"] = f_curv
    for name, o in CASES:
        vb = vectorial.make_vectorial_beam(mode="exact", **dict(VOPT, **o))
        out["map_" + name] = vb(X, Y)
        out["px_" + name] = vb(xp, 0.0 * xp)
        out["py_" + name] = vb(0.0 * xp, xp)
        out["depth_" + name] = vectorial.zero_depth(**dict(VOPT, **o))
        out["crb_" + name] = np.array([fisher.crb_limit(photons.make_model(patterns.tcp_centers(L), vb), C.N_REF)
                                       for L in L_arr])
    for key, f in (("lg300", C.FWHM), ("lgcurv", f_curv)):
        b = beams.make_beam("donut", fwhm=f)
        out["crb_" + key] = np.array([fisher.crb_limit(photons.make_model(patterns.tcp_centers(L), b), C.N_REF)
                                      for L in L_arr])
    return out


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    d = S.cached("fig2_vectorial", lambda: compute(args.quick), quick=args.quick, force=args.no_cache)
    g, xp, L = d["g"], d["xp"], d["L"]
    f_curv = float(d["f_curv"])

    fig, ax = S.figure("double", aspect=0.64, nrows=2, ncols=3)
    titles = {"correct": "circular, correct hand", "opposite": "circular, opposite hand",
              "linear": "linear (x)"}
    im = None
    for i, (name, _) in enumerate(CASES):
        a = ax[0, i]
        im = a.imshow(d["map_" + name], extent=(g[0], g[-1], g[0], g[-1]), origin="lower",
                      cmap=S.CMAP_SEQ, vmin=0, vmax=1, interpolation="bilinear")
        a.set_title(titles[name], fontsize=7.5, pad=3)
        a.set_xlabel("x (nm)")
        if i == 0:
            a.set_ylabel("y (nm)")
        else:
            a.set_yticklabels([])
        a.set_xticks([-400, 0, 400])
        a.set_yticks([-400, 0, 400])
        a.text(0.04, 0.05, r"$I(0)/I_{\max} = %.3f$" % float(d["depth_" + name]), transform=a.transAxes,
               color="w", fontsize=6.5)
        if name == "linear":
            a.annotate("", xy=(430, 400), xytext=(250, 400),
                       arrowprops=dict(arrowstyle="<->", color="w", lw=0.8))
            a.text(340, 440, "E", color="w", ha="center", fontsize=6.5)
        S.panel_label(a, "abc"[i])
    cb = fig.colorbar(im, ax=ax[0, :], shrink=0.92, pad=0.01, aspect=18)
    cb.set_label(r"$I/I_{\max}$", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    # (d) profiles through the zero
    a = ax[1, 0]
    a.plot(xp, beams.lg_donut(xp, 0.0 * xp, fwhm=C.FWHM), color="0.55", lw=0.9, label="LG 300 nm")
    a.plot(xp, d["px_correct"], color=S.COLORS["vec"], label="correct hand")
    a.plot(xp, d["px_opposite"], color=S.COLORS["opp"], label="opposite hand")
    a.plot(xp, d["px_linear"], color=S.COLORS["lin"], label="linear, along x")
    a.plot(xp, d["py_linear"], color=S.COLORS["lin"], ls="--", lw=1.0, label="linear, along y")
    a.set_xlim(-500, 500)
    a.set_ylim(0, 1.5)
    a.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    a.set_xlabel("x or y (nm)")
    a.set_ylabel(r"$I/I_{\max}$")
    a.legend(loc="upper center", ncol=2, fontsize=5.8, columnspacing=0.8, handlelength=1.3)
    S.panel_label(a, "d")

    # (e) centre CRB vs L
    a = ax[1, 1]
    a.loglog(L, d["crb_opposite"], "o-", color=S.COLORS["opp"], ms=2.5, label="opposite hand")
    a.loglog(L, d["crb_linear"], "s-", color=S.COLORS["lin"], ms=2.5, label="linear (x)")
    a.loglog(L, d["crb_lg300"], color="0.55", lw=0.9, label="LG 300 nm")
    a.loglog(L, d["crb_lgcurv"], color=S.COLORS["lg"], ls="--", lw=1.0, label="LG %.1f nm" % f_curv)
    a.loglog(L, d["crb_correct"], "^", color=S.COLORS["vec"], ms=3, label="correct hand")
    a.set_xlabel("L (nm)")
    a.set_ylabel(r"centre CRB, $r\to0$ (nm)")
    a.set_xlim(8, 300)
    a.set_ylim(0.2, 1000)
    a.legend(loc="center left", bbox_to_anchor=(0.0, 0.45), fontsize=5.8, handlelength=1.6)
    a.text(0.97, 0.04, "N = %d, no background" % C.N_REF, transform=a.transAxes, ha="right", fontsize=6)
    S.panel_label(a, "e")

    # (f) relative difference to LG 300
    a = ax[1, 2]
    rv = 100.0 * (d["crb_correct"] / d["crb_lg300"] - 1.0)
    rc = 100.0 * (d["crb_lgcurv"] / d["crb_lg300"] - 1.0)
    rvc = 100.0 * (d["crb_correct"] / d["crb_lgcurv"] - 1.0)
    a.plot(L, rv, "^-", color=S.COLORS["vec"], ms=3, label="vectorial / LG 300")
    a.plot(L, rc, "--", color=S.COLORS["lg"], lw=1.0, label="LG %.1f / LG 300" % f_curv)
    a.plot(L, rvc, "o:", color="0.3", ms=2.5, lw=0.8, label="vectorial / LG %.1f" % f_curv)
    a.axhline(0, color="0.6", lw=0.5)
    a.set_xlabel("L (nm)")
    a.set_ylabel("CRB difference (%)")
    a.set_xlim(0, 260)
    a.legend(loc="lower left", fontsize=5.8, handlelength=1.6)
    S.panel_label(a, "f")

    S.savefig(fig, "fig2_vectorial")

    S.report("vectorial_zero_depth_correct", float(d["depth_correct"]), "")
    S.report("vectorial_zero_depth_wrong_handedness", float(d["depth_opposite"]), "")
    S.report("vectorial_zero_depth_linear", float(d["depth_linear"]), "")
    S.report("vectorial_lg_fwhm_curvature_nm", f_curv, "nm")
    for Lr in C.L_LIST:
        j = int(np.argmin(np.abs(L - Lr)))
        tag = "L%d_N%d" % (Lr, C.N_REF)
        S.report("crb_center_vec_correct_%s_nm" % tag, float(d["crb_correct"][j]), "nm")
        S.report("crb_center_vec_opposite_%s_nm" % tag, float(d["crb_opposite"][j]), "nm")
        S.report("crb_center_vec_linear_%s_nm" % tag, float(d["crb_linear"][j]), "nm")
        S.report("crb_center_lg300_%s_nm" % tag, float(d["crb_lg300"][j]), "nm")
        S.report("crb_center_lgcurv_%s_nm" % tag, float(d["crb_lgcurv"][j]), "nm")
        S.report("crb_vec_over_lgcurv_%s" % tag, float(d["crb_correct"][j] / d["crb_lgcurv"][j]), "")
        S.report("crb_vec_over_lg300_%s" % tag, float(d["crb_correct"][j] / d["crb_lg300"][j]), "")
    S.write_summary("fig2")


if __name__ == "__main__":
    main()
