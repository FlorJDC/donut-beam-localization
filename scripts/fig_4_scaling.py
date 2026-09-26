# -*- coding: utf-8 -*-
"""Figure 4 -- scaling of the centre CRB with L, N and SBR, and comparison with a camera.

(a) Centre CRB vs L (N = 100, fwhm = 300 nm, no background): r -> 0 limit (numerical
    fisher.crb_limit and closed form) and point value of Eq. S27; log-log fit over L_FIT
    (L << fwhm). With eps = 0 the CRB increases monotonically only for L < 360 nm (ring
    diameter fwhm/sqrt(ln2)), diverges there and decreases beyond.
(b) Centre CRB vs N (L = 50 nm): limit, S27, and Eq. S31 with SBR = 10; fit over N_FIT.
(c) Centre CRB vs SBR (L = 50 nm, N = 100): Eq. S31 (point = limit for finite SBR) tends to the
    S27 point value, not to the background-free r -> 0 limit (the two limits do not commute).
(d) MINFLUX (L = 50 nm) vs camera localization: ideal sigma_PSF/sqrt(N) (sigma_PSF = 100 nm)
    and the pixelated 9 x 9 camera (100 nm pixels) without background and with SBR_c = 500
    defined per pixel (total SBR = 500/81); photons needed for 5 nm.

Usage: python scripts/fig_4_scaling.py [--quick] [--no-cache]
"""
import numpy as np

import _paperconfig as C
import _paperstyle as S
from donutloc import beams, camera, closed_forms, fisher, patterns, photons

SBR_CAM = 500.0          # Balzarotti SBR_c: total signal / background PER PIXEL
TARGET = 5.0             # nm


def _lim(L, N, sbr=None):
    b = beams.make_beam("donut", fwhm=C.FWHM)
    return fisher.crb_limit(photons.make_model(patterns.tcp_centers(L), b, sbr=sbr), N)


def _pt(L, N, sbr=None):
    b = beams.make_beam("donut", fwhm=C.FWHM)
    p = photons.make_model(patterns.tcp_centers(L), b, sbr=sbr)
    return float(fisher.crb(p, np.array([0.0, 0.0]), N, zero_policy="point"))


def compute(quick):
    out = {}
    # (a) vs L
    L_lo = np.logspace(0.0, np.log10(355.0), 25 if quick else 60)
    L_hi = np.linspace(365.0, 700.0, 10 if quick else 30)
    out["L_lo"], out["L_hi"] = L_lo, L_hi
    out["lim_L_lo"] = np.array([_lim(L, C.N_REF) for L in L_lo])
    out["lim_L_hi"] = np.array([_lim(L, C.N_REF) for L in L_hi])
    out["pt_L_lo"] = np.array([_pt(L, C.N_REF) for L in L_lo])
    out["lim_L_fit"] = np.array([_lim(L, C.N_REF) for L in C.L_FIT])
    # (b) vs N
    N_arr = np.logspace(0, 5, 21)
    out["N"] = N_arr
    out["lim_N"] = np.array([_lim(C.L_REF, N) for N in N_arr])
    out["pt_N"] = np.array([_pt(C.L_REF, N) for N in N_arr])
    out["sbr10_N"] = np.array([_pt(C.L_REF, N, C.SBR_MLE) for N in N_arr])
    out["lim_N_fit"] = np.array([_lim(C.L_REF, N) for N in C.N_FIT])
    # (c) vs SBR
    sbr = np.logspace(-0.5, 4.5, 41)
    out["sbr"] = sbr
    out["pt_sbr"] = np.array([_pt(C.L_REF, C.N_REF, s) for s in sbr])
    out["lim_sbr"] = np.array([_lim(C.L_REF, C.N_REF, s) for s in sbr])
    # (d) camera
    out["cam_pix_N"] = np.array([camera.crb_camera(C.SIGMA_PSF, N, pixel=C.CAM_PIXEL, n_pix=C.CAM_NPIX)
                                 for N in N_arr])
    out["cam_pixbg_N"] = np.array([camera.crb_camera(C.SIGMA_PSF, N, pixel=C.CAM_PIXEL, n_pix=C.CAM_NPIX,
                                                     sbr=SBR_CAM, sbr_convention="per_pixel")
                                   for N in N_arr])
    return out


def _n_for(sigma_at_N, N, target=TARGET):
    """Photons for ``target`` nm, exact for sigma proportional to 1/sqrt(N) (fixed SBR)."""
    return N * (sigma_at_N / target) ** 2


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    d = S.cached("fig4_scaling", lambda: compute(args.quick), quick=args.quick, force=args.no_cache)
    ring_d = 2 * beams.ring_radius(C.FWHM)

    slope_L, icpt_L = np.polyfit(np.log(C.L_FIT), np.log(d["lim_L_fit"]), 1)
    slope_N, icpt_N = np.polyfit(np.log(C.N_FIT), np.log(d["lim_N_fit"]), 1)

    fig, ax = S.figure("double", aspect=0.72, nrows=2, ncols=2)

    # (a) vs L
    a = ax[0, 0]
    L_lo, L_hi = d["L_lo"], d["L_hi"]
    a.axvspan(min(C.L_FIT), max(C.L_FIT), color="0.9", lw=0)
    a.loglog(L_lo, d["pt_L_lo"], color=S.COLORS["lms"], ls="--", lw=1.0, label=r"$r=0$ point (Eq. S27)")
    a.loglog(L_lo, d["lim_L_lo"], color=S.COLORS["lg"], label=r"$r\to0$ limit")
    a.loglog(L_hi, d["lim_L_hi"], color=S.COLORS["lg"], ls=":", lw=1.0)
    Lc = L_lo[L_lo < 355]
    a.loglog(Lc, closed_forms.crb_tcp_center_limit(Lc, C.N_REF, C.FWHM), "o", mfc="none",
             mec=S.COLORS["lg"], ms=2.2, mew=0.5, markevery=4, label="closed form")
    a.loglog(np.array(C.L_FIT), np.exp(icpt_L) * np.array(C.L_FIT) ** slope_L, color="k", lw=0.6)
    a.axvline(ring_d, color="0.5", lw=0.6, ls="--")
    a.text(ring_d * 0.93, 0.25, r"$L=\mathrm{fwhm}/\sqrt{\ln 2}$" + "\n= %.0f nm" % ring_d, fontsize=5.8,
           ha="right", va="bottom", color="0.3")
    a.text(14, 0.1, "fit %g-%g nm:\nslope %.3f" % (min(C.L_FIT), max(C.L_FIT), slope_L), fontsize=5.8,
           ha="left", va="top")
    a.set_xlim(1, 700)
    a.set_ylim(0.02, 300)
    a.set_xlabel("L (nm)")
    a.set_ylabel("centre CRB (nm)")
    a.legend(loc="upper left", fontsize=5.8)
    a.text(0.03, 0.62, "N = %d, no bkg." % C.N_REF, transform=a.transAxes, ha="left", fontsize=5.8)
    S.panel_label(a, "a")

    # (b) vs N
    a = ax[0, 1]
    N = d["N"]
    a.loglog(N, d["sbr10_N"], color=S.COLORS["sbr"], label="SBR = 10 (Eq. S31)")
    a.loglog(N, d["pt_N"], color=S.COLORS["lms"], ls="--", lw=1.0, label=r"$r=0$ point (Eq. S27)")
    a.loglog(N, d["lim_N"], color=S.COLORS["lg"], label=r"$r\to0$ limit")
    Nf = np.array(C.N_FIT, dtype=float)
    a.loglog(Nf, d["lim_N_fit"], "o", mfc="w", mec="k", ms=3, mew=0.7, label="fit points, slope %.4f" % slope_N)
    a.set_xlabel("detected photons N")
    a.set_ylabel("centre CRB (nm)")
    a.set_xlim(10, 1e5)
    a.legend(loc="upper right", fontsize=5.8)
    a.text(0.03, 0.04, "L = %d nm" % C.L_REF, transform=a.transAxes, fontsize=5.8)
    S.panel_label(a, "b")

    # (c) vs SBR
    a = ax[1, 0]
    sbr = d["sbr"]
    p27 = closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM)
    l0 = closed_forms.crb_tcp_center_limit(C.L_REF, C.N_REF, C.FWHM)
    a.semilogx(sbr, closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM, sbr=sbr), color=S.COLORS["sbr"],
               label="Eq. S31")
    a.semilogx(sbr, d["pt_sbr"], "o", mfc="none", mec="k", ms=2.5, mew=0.5, markevery=2,
               label="numerical, r = 0")
    a.semilogx(sbr, d["lim_sbr"], "x", color="0.4", ms=2.5, mew=0.5, markevery=[i for i in range(1, sbr.size, 2)],
               label=r"numerical, $r\to0$")
    a.axhline(p27, color=S.COLORS["lms"], ls="--", lw=0.8)
    a.axhline(l0, color=S.COLORS["lg"], ls="-", lw=0.8)
    a.text(sbr[-1] * 0.8, p27 + 0.04, "S27 (no bkg., r = 0)", ha="right", va="bottom", fontsize=5.8,
           color=S.COLORS["lms"])
    a.text(sbr[-1] * 0.8, l0 - 0.04, r"no bkg., $r\to0$", ha="right", va="top", fontsize=5.8, color=S.COLORS["lg"])
    s10 = closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM, sbr=C.SBR_MLE)
    a.plot(C.SBR_MLE, s10, "s", color=S.COLORS["sbr"], ms=3.5)
    a.annotate("SBR = 10: %.3f nm" % s10, (C.SBR_MLE, s10), xytext=(8, 10), textcoords="offset points",
               fontsize=5.8)
    a.set_xlim(sbr[0], sbr[-1])
    a.set_ylim(1.4, 5.0)
    a.set_xlabel("SBR (total signal / total background)")
    a.set_ylabel("centre CRB (nm)")
    a.legend(loc="upper right", fontsize=5.8)
    a.text(0.3, 0.93, "L = %d nm, N = %d" % (C.L_REF, C.N_REF), transform=a.transAxes, fontsize=5.8,
           va="top")
    S.panel_label(a, "c")

    # (d) camera vs MINFLUX
    a = ax[1, 1]
    cam_ideal = camera.crb_camera_ideal(C.SIGMA_PSF, N)
    a.loglog(N, cam_ideal, color=S.COLORS["cam"], label=r"camera, ideal $\sigma_{PSF}/\sqrt{N}$")
    # open squares: the pixelated no-background curve lies about 4 % above the ideal line
    a.loglog(N, d["cam_pix_N"], color=S.COLORS["cam"], ls="--", lw=0.7, marker="s", ms=2.4, mfc="w",
             mew=0.5, markevery=2, label="camera 9x9 px, no bkg.")
    a.loglog(N, d["cam_pixbg_N"], color=S.COLORS["cam"], ls=":", lw=1.0,
             label=r"camera 9x9 px, SBR$_c$ = 500 per px")
    a.loglog(N, d["sbr10_N"], color=S.COLORS["sbr"], label="MINFLUX, SBR = 10")
    a.loglog(N, d["pt_N"], color=S.COLORS["lms"], ls="--", lw=1.0, label="MINFLUX, S27")
    a.loglog(N, d["lim_N"], color=S.COLORS["lg"], label=r"MINFLUX, $r\to0$")
    a.axhline(TARGET, color="0.6", lw=0.5)
    n5 = {
        "cam_ideal": _n_for(float(camera.crb_camera_ideal(C.SIGMA_PSF, 400.0)), 400.0),
        "cam_pix": _n_for(float(np.interp(np.log(400.0), np.log(N), d["cam_pix_N"])), 400.0),
        "cam_pixbg": _n_for(camera.crb_camera(C.SIGMA_PSF, 600.0, pixel=C.CAM_PIXEL, n_pix=C.CAM_NPIX,
                                              sbr=SBR_CAM, sbr_convention="per_pixel"), 600.0),
        "mf_s27": _n_for(closed_forms.crb_tcp_center_point(C.L_REF, C.N_REF, C.FWHM), C.N_REF),
        "mf_lim": _n_for(_lim(C.L_REF, C.N_REF), C.N_REF),
        "mf_sbr10": _n_for(s10, C.N_REF),
    }
    # exact pixelated value at N = 400 (the interpolated one is only for plotting)
    n5["cam_pix"] = _n_for(camera.crb_camera(C.SIGMA_PSF, 400.0, pixel=C.CAM_PIXEL, n_pix=C.CAM_NPIX), 400.0)
    for k, col in (("cam_ideal", S.COLORS["cam"]), ("cam_pixbg", S.COLORS["cam"]), ("mf_lim", S.COLORS["lg"]),
                   ("mf_s27", S.COLORS["lms"]), ("mf_sbr10", S.COLORS["sbr"])):
        a.plot(n5[k], TARGET, "v", color=col, ms=3.0, mec="w", mew=0.3, zorder=5)
    a.text(n5["mf_lim"] * 0.8, TARGET * 0.72, "%.1f" % n5["mf_lim"], fontsize=5.5, ha="center", va="top",
           color=S.COLORS["lg"])
    a.text(n5["mf_s27"] * 1.3, TARGET * 0.72, "%.0f" % n5["mf_s27"], fontsize=5.5, ha="center", va="top",
           color=S.COLORS["lms"])
    a.text(n5["mf_sbr10"] * 1.25, TARGET * 1.12, "%.1f" % n5["mf_sbr10"], fontsize=5.5, ha="center", va="bottom",
           color=S.COLORS["sbr"])
    a.text(n5["cam_ideal"] * 0.85, TARGET * 0.72, "%.0f" % n5["cam_ideal"], fontsize=5.5, ha="center", va="top")
    a.text(n5["cam_pixbg"] * 1.2, TARGET * 1.12, "%.0f" % n5["cam_pixbg"], fontsize=5.5, ha="center", va="bottom")
    a.text(1.3, TARGET * 1.1, "5 nm", fontsize=5.5, color="0.4", va="bottom")
    a.set_xlim(1, 1e5)
    a.set_ylim(0.004, 150)
    a.set_xlabel("detected photons N")
    a.set_ylabel("localization precision (nm)")
    a.legend(loc="lower left", fontsize=5.5, ncol=2, columnspacing=0.8, handlelength=1.6)
    S.panel_label(a, "d")

    S.savefig(fig, "fig4_scaling")

    S.report("crb_exponent_L", float(slope_L), "")
    S.report("crb_exponent_N", float(slope_N), "")
    S.report("crb_center_lg_L50_N100_nm", float(_lim(C.L_REF, C.N_REF)), "nm")
    S.report("crb_center_point_S27_L50_N100_nm", float(p27), "nm")
    S.report("crb_center_sbr10_L50_N100_nm", float(s10), "nm")
    S.report("lg_ring_diameter_nm", float(ring_d), "nm")
    S.report("fig4_crb_limit_L355_N100_nm", float(d["lim_L_lo"][-1]), "nm")
    S.report("fig4_crb_limit_L365_N100_nm", float(d["lim_L_hi"][0]), "nm")
    S.report("fig4_crb_limit_L700_N100_nm", float(d["lim_L_hi"][-1]), "nm")
    S.report("camera_sigma_nm", float(camera.crb_camera_ideal(C.SIGMA_PSF, C.ITER_N_TOTAL)), "nm")
    S.report("camera_pixelated_9x9_N400_nm", float(camera.crb_camera(C.SIGMA_PSF, 400.0, pixel=C.CAM_PIXEL,
                                                             n_pix=C.CAM_NPIX)), "nm")
    S.report("camera_perpixel_sbr500_N600_nm", float(camera.crb_camera(C.SIGMA_PSF, 600.0, pixel=C.CAM_PIXEL,
                                                                     n_pix=C.CAM_NPIX, sbr=SBR_CAM,
                                                                     sbr_convention="per_pixel")), "nm")
    S.report("camera_sbr_perpixel_to_total", SBR_CAM / C.CAM_NPIX ** 2, "")
    S.report("fig4_photons_for_5nm_camera_ideal", float(n5["cam_ideal"]), "photons")
    S.report("fig4_photons_for_5nm_camera_pixelated_9x9", float(n5["cam_pix"]), "photons")
    S.report("camera_photons_for_5nm_perpixel_sbr500", float(n5["cam_pixbg"]), "photons")
    S.report("photons_for_5nm_L50_sbrinf", float(n5["mf_s27"]), "photons")
    S.report("photons_for_5nm_L50_limit_nobg", float(n5["mf_lim"]), "photons")
    S.report("photons_for_5nm_L50_sbr10", float(n5["mf_sbr10"]), "photons")
    S.write_summary("fig4")


if __name__ == "__main__":
    main()
