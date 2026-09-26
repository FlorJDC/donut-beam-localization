# -*- coding: utf-8 -*-
"""Figure 1 -- beam model and TCP geometry.

(a) Normalized radial profiles: scalar LG01 donut (fwhm = 300 nm, Balzarotti Eq. S17) vs the
    vectorial Richards-Wolf donut with the correct circular handedness (lambda = 640 nm,
    NA = 1.4, n = 1.518, filling 5/3), and the curvature-matched LG (fwhm = 327.14 nm).
(b) Zoom of the zero on a log radius axis: I/rho^2 -> curvature c at the zero (quadratic zeros).
(c) TCP geometry: 3 donut zeros on a circle of diameter L plus the centre; map of one exposure.

Usage: python scripts/fig_1_schematic.py [--quick] [--no-cache]
"""
import numpy as np

import _paperconfig as C
import _paperstyle as S
from donutloc import beams, patterns, vectorial

VOPT = dict(wavelength=C.WAVELENGTH, NA=C.NA, n=C.N_MEDIUM, filling=C.FILLING)


def compute():
    rho = np.linspace(0.0, 700.0, 701)
    rho_z = np.logspace(0.0, np.log10(300.0), 120)
    vb = vectorial.make_vectorial_beam(mode="exact", polarization="circular", handedness=+1, **VOPT)
    f_curv = vectorial.lg_equivalent_fwhm(dict(VOPT), match="curvature")
    f_diam = vectorial.lg_equivalent_fwhm(dict(VOPT), match="diameter")
    return dict(
        rho=rho, rho_z=rho_z,
        I_vec=vb(rho, 0.0 * rho), I_vec_z=vb(rho_z, 0.0 * rho_z),
        I_vec0=float(vb(np.array(0.0), np.array(0.0))),
        f_curv=f_curv, f_diam=f_diam,
        d_pp=vectorial.peak_to_peak_diameter(**VOPT),
        curv=vectorial.zero_curvature(**VOPT),
    )


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    d = S.cached("fig1_profiles", compute, quick=args.quick, force=args.no_cache)
    rho, rho_z = d["rho"], d["rho_z"]
    f_curv = float(d["f_curv"])
    lg = lambda r, f: beams.lg_donut(r, 0.0 * r, fwhm=f)  # noqa: E731

    fig, ax = S.figure("double", aspect=0.34, ncols=3, gridspec_kw=dict(width_ratios=[1.15, 1.0, 1.0]))

    # (a) radial profiles
    a = ax[0]
    a.plot(rho, lg(rho, C.FWHM), color=S.COLORS["lg"], label="LG %.0f nm" % C.FWHM)
    a.plot(rho, lg(rho, f_curv), color=S.COLORS["lg"], ls="--", lw=1.0,
           label="LG %.1f nm\n(curv. matched)" % f_curv)
    a.plot(rho, d["I_vec"], color=S.COLORS["vec"], label="vectorial\n(circ., correct hand)")
    a.axvline(beams.ring_radius(C.FWHM), color="0.6", lw=0.5, ls=":")
    a.axvline(0.5 * float(d["d_pp"]), color=S.COLORS["vec"], lw=0.5, ls=":")
    a.set_xlim(0, 700)
    a.set_ylim(0, 1.08)
    a.set_xlabel(r"radius $\rho$ (nm)")
    a.set_ylabel("intensity / ring peak")
    a.legend(loc="upper left", bbox_to_anchor=(0.55, 1.0), fontsize=5.8, borderaxespad=0.2, handlelength=1.4)
    S.panel_label(a, "a")

    # (b) zoom of the zero: I / rho^2 on a log radius axis (curvature at the zero)
    b = ax[1]
    c = float(d["curv"])
    b.semilogx(rho_z, lg(rho_z, C.FWHM) / rho_z ** 2 * 1e5, color=S.COLORS["lg"], label="LG 300 nm")
    b.semilogx(rho_z, lg(rho_z, f_curv) / rho_z ** 2 * 1e5, color=S.COLORS["lg"], ls="--", lw=1.0,
               label="LG %.1f nm" % f_curv)
    b.semilogx(rho_z, d["I_vec_z"] / rho_z ** 2 * 1e5, color=S.COLORS["vec"], label="vectorial")
    b.axhline(c * 1e5, color="0.4", lw=0.6, ls=":")
    b.set_xlim(1, 300)
    b.set_ylim(0, 10)
    b.set_xlabel(r"radius $\rho$ (nm)")
    b.set_ylabel(r"$I/\rho^2$ ($10^{-5}$ nm$^{-2}$)")
    b.legend(loc="lower left", fontsize=6.3)
    b.text(1.3, c * 1e5 + 0.25, r"vectorial: $c = %.3f\times10^{-5}$ nm$^{-2}$" % (c * 1e5),
           fontsize=6, va="bottom", color="0.3")
    S.panel_label(b, "b")

    # (c) TCP geometry
    e = ax[2]
    Lg = 100.0
    cen = patterns.tcp_centers(Lg)
    g = np.linspace(-200.0, 200.0, 401)
    X, Y = np.meshgrid(g, g)
    Imap = beams.lg_donut(X - cen[0, 0], Y - cen[0, 1], fwhm=C.FWHM)
    im = e.imshow(Imap, extent=(g[0], g[-1], g[0], g[-1]), origin="lower", cmap=S.CMAP_SEQ,
                  vmin=0, vmax=1, interpolation="bilinear")
    t = np.linspace(0, 2 * np.pi, 200)
    e.plot(0.5 * Lg * np.cos(t), 0.5 * Lg * np.sin(t), color="w", lw=0.7, ls="--")
    for k, (x, y) in enumerate(cen):
        e.plot(x, y, marker="x" if k < 3 else "+", color="w", ms=5, mew=1.1)
        off = (0, 8) if k == 0 else ((-8, -6) if k == 1 else ((8, -6) if k == 2 else (0, -8)))
        e.annotate(str(k), (x, y), xytext=off, textcoords="offset points", color="w",
                   fontsize=7, ha="center", va="center")
    e.annotate("", xy=(-0.5 * Lg, -95), xytext=(0.5 * Lg, -95),
               arrowprops=dict(arrowstyle="<->", color="w", lw=0.7))
    e.text(0, -108, "L = %.0f nm" % Lg, color="w", ha="center", va="top", fontsize=6.5)
    e.set_xlabel("x (nm)")
    e.set_ylabel("y (nm)")
    e.set_xticks([-200, -100, 0, 100, 200])
    e.set_yticks([-200, -100, 0, 100, 200])
    cb = fig.colorbar(im, ax=e, shrink=0.9, pad=0.02)
    cb.set_label("intensity of exposure 0", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    S.panel_label(e, "c", x=-0.22)

    S.savefig(fig, "fig1_schematic")

    S.report("fwhm_nm", C.FWHM, "nm")
    S.report("lg_ring_diameter_nm", 2 * beams.ring_radius(C.FWHM), "nm")
    S.report("vectorial_peak_to_peak_diameter_nm", float(d["d_pp"]), "nm")
    S.report("vectorial_zero_curvature_nm2", c, "nm^-2")
    S.report("vectorial_lg_fwhm_curvature_nm", f_curv, "nm")
    S.report("vectorial_lg_fwhm_diameter_nm", float(d["f_diam"]), "nm")
    rz = d["rho_z"]
    dev = np.abs(d["I_vec_z"] / lg(rz, f_curv) - 1.0)[rz <= 100.0].max()
    S.report("fig1_vec_vs_lgcurv_max_rel_dev_rho_le_100nm", float(dev), "")
    S.report("fig1_tcp_L_nm", Lg, "nm")
    S.write_summary("fig1")


if __name__ == "__main__":
    main()
