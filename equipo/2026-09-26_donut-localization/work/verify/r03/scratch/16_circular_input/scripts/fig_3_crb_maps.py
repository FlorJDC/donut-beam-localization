# -*- coding: utf-8 -*-
"""Figure 3 -- CRB maps of the TCP and the discontinuity at the centre.

(a-d) fisher.crb_map (zero_policy="limit") over the field of view |x|, |y| <= L for L = 50 and
      100 nm, N = 100, fwhm = 300 nm, without background and with SBR = 10.
(e)   Cut along x (and y) through the centre, L = 50 nm: without background the r -> 0 limit
      (1.605 nm) differs from the point value at r = 0 of Balzarotti Eq. S27 (1.802 nm), which
      excludes the p = 0 central exposure; with SBR = 10 the CRB is continuous (Eq. S31).
(f)   Ratio limit/point vs L (closed_forms.limit_to_point_ratio) with the quadratic-zero
      asymptote 2/sqrt(5), and numerical fisher.crb_limit / fisher.crb(r = 0) checks.

Usage: python scripts/fig_3_crb_maps.py [--quick] [--no-cache]
"""
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter

import _paperconfig as C
import _paperstyle as S
from donutloc import beams, closed_forms, fisher, patterns, photons

L_MAPS = (50.0, 100.0)
SBRS = (None, 10.0)


def compute(quick):
    b = beams.make_beam("donut", fwhm=C.FWHM)
    n = 61 if quick else 161
    out = {}
    for L in L_MAPS:
        xs = np.linspace(-L, L, n)                  # odd n: the origin is a grid node
        out["xs_%d" % L] = xs
        for s in SBRS:
            p = photons.make_model(patterns.tcp_centers(L), b, sbr=s)
            out["map_%d_%s" % (L, s)] = fisher.crb_map(p, xs, xs, C.N_REF, zero_policy="limit")
    # cut through the centre, L = 50 (even number of points: r = 0 excluded)
    L = C.L_REF
    xc = np.linspace(-L, L, 400)
    zeros = np.zeros_like(xc)
    for s in SBRS:
        p = photons.make_model(patterns.tcp_centers(L), b, sbr=s)
        out["cutx_%s" % s] = fisher.crb(p, np.stack([xc, zeros], -1), C.N_REF, zero_policy="point")
        out["cuty_%s" % s] = fisher.crb(p, np.stack([zeros, xc], -1), C.N_REF, zero_policy="point")
        out["point0_%s" % s] = float(fisher.crb(p, np.array([0.0, 0.0]), C.N_REF, zero_policy="point"))
        out["limit0_%s" % s] = fisher.crb_limit(p, C.N_REF)
    out["xc"] = xc
    # numerical limit/point ratio
    Lchk = np.array([5.0, 25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0])
    rnum = []
    for Lk in Lchk:
        p = photons.make_model(patterns.tcp_centers(Lk), b)
        rnum.append(fisher.crb_limit(p, C.N_REF) / float(fisher.crb(p, np.array([0.0, 0.0]), C.N_REF,
                                                                       zero_policy="point")))
    out["L_chk"] = Lchk
    out["ratio_num"] = np.array(rnum)
    return out


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    d = S.cached("fig3_crb_maps", lambda: compute(args.quick), quick=args.quick, force=args.no_cache)

    fig, ax = S.figure("double", aspect=0.66, nrows=2, ncols=3,
                       gridspec_kw=dict(width_ratios=[1.0, 1.0, 1.25]))
    letters = iter("abcd")
    for i, L in enumerate(L_MAPS):
        xs = d["xs_%d" % L]
        c0 = float(d["map_%d_None" % L][xs.size // 2, xs.size // 2])
        vmin = min(float(d["map_%d_%s" % (L, s)].min()) for s in SBRS)
        vmax = max(float(d["map_%d_%s" % (L, s)].max()) for s in SBRS)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        im = None
        cen = patterns.tcp_centers(L)
        for j, s in enumerate(SBRS):
            a = ax[i, j]
            m = d["map_%d_%s" % (L, s)]
            im = a.imshow(m, extent=(xs[0], xs[-1], xs[0], xs[-1]), origin="lower", cmap=S.CMAP_SEQ,
                          norm=norm, interpolation="bilinear")
            levels = [lv for lv in (2, 3, 5, 10, 20, 50, 100) if vmin < lv < vmax]
            cs = a.contour(xs, xs, m, levels=levels, colors="w", linewidths=0.4, alpha=0.8)
            a.clabel(cs, fmt="%g", fontsize=5, inline_spacing=1)
            t = np.linspace(0, 2 * np.pi, 200)
            a.plot(0.5 * L * np.cos(t), 0.5 * L * np.sin(t), color="w", lw=0.5, ls="--")
            a.plot(cen[:3, 0], cen[:3, 1], "x", color=S.COLORS["vec"], ms=4, mew=1.0)
            a.plot(0, 0, "+", color=S.COLORS["vec"], ms=5, mew=1.0)
            a.set_title("L = %d nm, %s" % (L, "no background" if s is None else "SBR = %g" % s),
                        fontsize=7, pad=3)
            a.set_xlabel("x (nm)")
            if j == 0:
                a.set_ylabel("y (nm)")
            else:
                a.set_yticklabels([])
            S.panel_label(a, next(letters))
        cb = fig.colorbar(im, ax=list(ax[i, :2]), shrink=0.95, pad=0.02, aspect=14)
        cb.set_label("CRB (nm)", fontsize=7)
        cb.ax.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
        cb.ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
        cb.ax.yaxis.set_minor_formatter(NullFormatter())
        cb.ax.tick_params(labelsize=6)

    # (e) cut through the centre
    a = ax[0, 2]
    xc = d["xc"]
    h1, = a.plot(xc, d["cutx_None"], color=S.COLORS["lg"], label="no bkg.")
    a.plot(xc, d["cuty_None"], color=S.COLORS["lg"], ls="--", lw=0.9)
    h2, = a.plot(xc, d["cutx_10.0"], color=S.COLORS["sbr"], label="SBR = 10")
    a.plot(xc, d["cuty_10.0"], color=S.COLORS["sbr"], ls="--", lw=0.9)
    # point values from the closed forms (Eq. S27 / S31); the numerical r = 0 values agree to ~2e-9
    p0 = closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM)
    s0 = closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM, sbr=C.SBR_MLE)
    l0 = float(d["limit0_None"])
    h3, = a.plot(0, p0, "o", mfc="w", mec="k", ms=4, mew=0.9, zorder=5, label="r = 0 (S27)")
    h4, = a.plot(0, l0, "o", color="k", ms=3.5, zorder=5, label=r"$r\to0$ limit")
    a.plot(0, s0, "s", mfc="w", mec=S.COLORS["sbr"], ms=3.5, mew=0.9, zorder=5)
    a.set_xlim(-50, 50)
    a.set_ylim(1.5, 2.7)
    a.set_xlabel("position through the centre (nm)")
    a.set_ylabel("CRB (nm)")
    hx, = a.plot([], [], color="0.3", label="along x")
    hy, = a.plot([], [], color="0.3", ls="--", lw=0.9, label="along y")
    leg = a.legend(handles=[h1, h2, hx, hy], loc="upper left", fontsize=5.8, handlelength=1.5,
                   borderaxespad=0.3)
    a.add_artist(leg)
    a.legend(handles=[h3, h4], loc="lower right", fontsize=5.8, handlelength=1.0, borderaxespad=0.3)
    a.text(0.03, 0.04, "L = %d nm\nN = %d" % (C.L_REF, C.N_REF), transform=a.transAxes, fontsize=5.8,
           va="bottom")
    S.panel_label(a, "e")

    # (f) limit / point ratio vs L
    a = ax[1, 2]
    Lg = np.linspace(1.0, 355.0, 400)
    a.plot(Lg, closed_forms.limit_to_point_ratio(Lg, C.FWHM), color=S.COLORS["lg"],
           label="closed form, fwhm = %d nm" % C.FWHM)
    a.plot(Lg, closed_forms.limit_to_point_ratio(Lg, C.FWHM_MASULLO), color=S.COLORS["lms"], ls="--", lw=1.0,
           label="closed form, fwhm = %d nm" % C.FWHM_MASULLO)
    a.plot(d["L_chk"], d["ratio_num"], "o", mfc="w", mec="k", ms=3.5, mew=0.8, label="numerical (fwhm = 300 nm)")
    a.axhline(2 / np.sqrt(5), color="0.4", lw=0.6, ls=":")
    a.text(150, 2 / np.sqrt(5) + 0.006, r"$2/\sqrt{5}$ (quadratic zero)", fontsize=6, color="0.3")
    a.set_xlim(0, 360)
    a.set_ylim(0.64, 0.93)
    a.set_xlabel("L (nm)")
    a.set_ylabel(r"$\sigma_{r\to0}\;/\;\sigma_{r=0}$ (S27)")
    a.legend(loc="lower left", fontsize=5.8)
    S.panel_label(a, "f")

    S.savefig(fig, "fig3_crb_maps")

    S.report("crb_center_lg_L50_N100_nm", float(d["limit0_None"]), "nm")
    p27 = closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM)
    s31 = closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM, sbr=C.SBR_MLE)
    S.report("crb_center_point_S27_L50_N100_nm", float(p27), "nm")
    S.report("crb_center_sbr10_L50_N100_nm", float(s31), "nm")
    S.report("fig3_point0_numeric_rel_dev_S27", abs(float(d["point0_None"]) / p27 - 1.0), "")
    S.report("fig3_point0_numeric_rel_dev_S31_sbr10", abs(float(d["point0_10.0"]) / s31 - 1.0), "")
    S.report("crb_center_sbr10_limit_L50_N100_nm", float(d["limit0_10.0"]), "nm")
    for L in L_MAPS:
        xs = d["xs_%d" % L]
        for s in SBRS:
            m = d["map_%d_%s" % (L, s)]
            tag = "L%d_N%d_%s" % (L, C.N_REF, "nobkg" if s is None else "sbr%g" % s)
            S.report("fig3_map_center_%s_nm" % tag, float(m[xs.size // 2, xs.size // 2]), "nm")
            S.report("fig3_map_min_%s_nm" % tag, float(m.min()), "nm")
            inside = np.hypot(*np.meshgrid(xs, xs)) <= 0.5 * L
            S.report("fig3_map_max_inside_tcp_%s_nm" % tag, float(m[inside].max()), "nm")
    for Lr in (5.0, 50.0, 100.0, 150.0):
        S.report("limit_to_point_ratio_L%d" % Lr, float(closed_forms.limit_to_point_ratio(Lr, C.FWHM)), "")
    S.report("limit_to_point_ratio_quadratic", 2 / np.sqrt(5), "")
    dev = np.max(np.abs(d["ratio_num"] / closed_forms.limit_to_point_ratio(d["L_chk"], C.FWHM) - 1))
    S.report("fig3_ratio_numeric_vs_closed_max_rel_dev", float(dev), "")
    S.write_summary("fig3")


if __name__ == "__main__":
    main()
