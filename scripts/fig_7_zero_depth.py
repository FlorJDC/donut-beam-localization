# -*- coding: utf-8 -*-
"""Figure 7 -- finite depth of the donut zero (residual intensity eps at the centre).

LG donut (fwhm = 300 nm) with a Gaussian pedestal of relative height eps (zero_model =
"gaussian"), TCP centre, N = 100 photons, no background (scripts/_paperconfig.py).

(a) centre CRB versus L for eps = 0 and EPS_LIST (experiments.eps_L_sweep), with the optimum
    L_opt of each eps marked (experiments.optimal_L).  With eps = 0 the centre value is the
    r -> 0 limit; it grows monotonically with L only for L < fwhm/sqrt(ln 2) = 360 nm.
(b) L_opt / (fwhm sqrt(eps)) versus eps: the scaling L_opt ~ 0.78 fwhm sqrt(eps) holds only for
    eps <~ 0.01.
(c) eps = 0.002, L = 100 nm: CRB versus distance from the centre along three directions; the
    minimum is not at the centre.  sqrt(eps/(e a)) = 4.9 nm (a = 4 ln2/fwhm^2) is marked as the
    transition SCALE, not as the argmin (which depends on L and on the direction).

Deterministic (no Monte Carlo); results cached in data/mc/fig7_*.npz like the other figures.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _paperconfig as C  # noqa: E402
import _paperstyle as S  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

from donutloc import beams, experiments, fisher, patterns, photons  # noqa: E402

N, FWHM = C.N_REF, C.FWHM
ZM = "gaussian"
L_RING = FWHM / np.sqrt(np.log(2.0))          # ring diameter where eps=0 diverges (360.7 nm)
EPS_C, L_C = C.EPS_TRANSITION, 100.0                    # panel (c), V10


def compute(quick):
    eps_all = np.array((0.0,) + tuple(C.EPS_LIST))
    Ls = np.geomspace(2.0, 340.0, 60 if quick else 160)
    sw = experiments.eps_L_sweep(eps_all, Ls, N=N, fwhm=FWHM, zero_model=ZM)
    lopt = np.array([experiments.optimal_L(e, N=N, fwhm=FWHM, zero_model=ZM)["L_opt"]
                     for e in C.EPS_LIST])
    copt = np.array([experiments.optimal_L(e, N=N, fwhm=FWHM, zero_model=ZM)["crb_opt"]
                     for e in C.EPS_LIST])
    eps_d = np.geomspace(1e-4, 0.3, 12 if quick else 30)
    lopt_d = np.array([experiments.optimal_L(e, N=N, fwhm=FWHM, zero_model=ZM)["L_opt"]
                       for e in eps_d])
    # (c): CRB versus |r| along 0, 45 and 90 degrees, eps = 0.002, L = 100
    beam = beams.make_beam("donut", fwhm=FWHM, eps=EPS_C, zero_model=ZM)
    p_fn = photons.make_model(patterns.tcp_centers(L_C), beam)
    rr = np.linspace(0.0, 20.0, 41 if quick else 201)
    ang = np.deg2rad([0.0, 45.0, 90.0])
    pts = rr[None, :, None] * np.stack([np.cos(ang), np.sin(ang)], -1)[:, None, :]
    crb_r = np.asarray(fisher.crb(p_fn, pts, N), float)
    return {"eps": eps_all, "L": Ls, "crb": sw["crb_center"], "eps_opt": np.array(C.EPS_LIST),
            "L_opt": lopt, "crb_opt": copt, "eps_dense": eps_d, "L_opt_dense": lopt_d,
            "r": rr, "ang_deg": np.array([0.0, 45.0, 90.0]), "crb_r": crb_r}


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    t0 = time.time()
    d = S.cached("fig7_zero_depth", lambda: compute(args.quick), args.quick, args.no_cache)

    a = 4.0 * np.log(2.0) / FWHM ** 2
    scale = np.sqrt(EPS_C / (np.e * a))
    ratio = d["L_opt"] / (FWHM * np.sqrt(d["eps_opt"]))
    for e, lo, co, q in zip(d["eps_opt"], d["L_opt"], d["crb_opt"], ratio):
        tag = ("%g" % e).replace(".", "p")
        S.report("L_opt_eps%s_nm" % tag, lo, "nm")
        S.report("crb_opt_eps%s_nm" % tag, co, "nm")
        S.report("L_opt_over_fwhm_sqrt_eps_eps%s" % tag, q, "")
    S.report("L_opt_over_fwhm_sqrt_eps_small_eps_limit", float(
        d["L_opt_dense"][0] / (FWHM * np.sqrt(d["eps_dense"][0]))), "")
    S.report("fig7_eps_dense_min", float(d["eps_dense"][0]), "")
    S.report("eps0_crb_divergence_L_nm", L_RING, "nm")
    S.report("fig7_crb_eps0_L50_nm", experiments.crb_center(50.0, N, FWHM, 0.0, zero_model=ZM),
             "nm")
    for k, ang in enumerate(d["ang_deg"]):
        i = int(np.argmin(d["crb_r"][k]))
        S.report("fig7_eps0p002_L100_crb_min_dir%d_nm" % ang, float(d["crb_r"][k, i]), "nm")
        S.report("fig7_eps0p002_L100_argmin_dir%d_nm" % ang, float(d["r"][i]), "nm")
    S.report("fig7_eps0p002_L100_crb_center_nm", float(d["crb_r"][0, 0]), "nm")
    S.report("zero_depth_transition_scale_nm", float(scale), "nm")

    # ------------------------------------------------------------------ figure
    fig, ax = S.figure("double", aspect=0.34, nrows=1, ncols=3)
    cols = ["#000000"] + [S.CYCLE[i] for i in (0, 2, 4, 1)]
    for i, e in enumerate(d["eps"]):
        lab = r"$\epsilon=0$ ($r\to0$ limit)" if e == 0 else r"$\epsilon=%g$" % e
        ax[0].plot(d["L"], d["crb"][i], color=cols[i], label=lab, lw=1.1,
                   ls="--" if e == 0 else "-")
    for i in range(d["eps_opt"].size):
        ax[0].plot(d["L_opt"][i], d["crb_opt"][i], marker="*", ms=7, color=cols[i + 1],
                   mec="k", mew=0.4, ls="none", zorder=5)
    ax[0].plot([], [], marker="*", ms=7, color="0.7", mec="k", mew=0.4, ls="none",
               label=r"$L_\mathrm{opt}$")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"TCP diameter $L$ (nm)")
    ax[0].set_ylabel(r"centre CRB $\sigma_\mathrm{CRB}$ (nm), $N=%d$" % N)
    ax[0].legend(loc="upper center", ncol=2, fontsize=6.2, handlelength=1.5, columnspacing=0.8)
    ax[0].set_ylim(0.03, 4000)
    ax[0].xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
    ax[0].yaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))

    rd = d["L_opt_dense"] / (FWHM * np.sqrt(d["eps_dense"]))
    ax[1].plot(d["eps_dense"], rd, color=S.COLORS["lg"], lw=1.1)
    for i in range(d["eps_opt"].size):
        ax[1].plot(d["eps_opt"][i], ratio[i], marker="*", ms=7, color=cols[i + 1], mec="k",
                   mew=0.4, ls="none", zorder=5)
    ax[1].axhline(0.78, color="0.4", ls=":", lw=0.8)
    ax[1].text(0.25, 0.777, "0.78", fontsize=6.5, color="0.3", va="top", ha="right")
    ax[1].axvspan(1e-4, 0.01, color="0.93", lw=0, zorder=0)
    ax[1].text(1e-3, 0.25, "0.78 scaling\nvalid here",
               transform=ax[1].get_xaxis_transform(), ha="center", va="bottom", fontsize=6.2,
               color="0.3")
    ax[1].set_xscale("log")
    ax[1].set_xlim(1e-4, 0.3)
    ax[1].set_xlabel(r"residual zero $\epsilon = I(0)/I_\mathrm{max}$")
    ax[1].set_ylabel(r"$L_\mathrm{opt}\,/\,(\mathrm{fwhm}\sqrt{\epsilon})$")

    lss = ["-", "--", ":"]
    for k, ang in enumerate(d["ang_deg"]):
        ax[2].plot(d["r"], d["crb_r"][k], color=S.COLORS["lg"], ls=lss[k],
                   label=r"along %d$^\circ$" % ang)
    ax[2].axvline(scale, color="0.5", lw=0.7, ls="-.")
    ax[2].text(scale, 0.97, r" $\sqrt{\epsilon/(e\,a)}$ = %.1f nm" % scale, fontsize=6.2,
               transform=ax[2].get_xaxis_transform(), va="top", color="0.3")
    ax[2].set_xlabel(r"distance from centre $|\mathbf{r}|$ (nm)")
    ax[2].set_ylabel(r"$\sigma_\mathrm{CRB}$ (nm)")
    ax[2].set_title(r"$\epsilon=%g$, $L=%g$ nm, $N=%d$" % (EPS_C, L_C, N), fontsize=7)
    ax[2].legend(loc="lower right", fontsize=6.2)
    for axx, l in zip(ax, "abc"):
        S.panel_label(axx, l)
    S.savefig(fig, "fig7_zero_depth")
    S.write_summary("fig7")
    print("fig7: total %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
