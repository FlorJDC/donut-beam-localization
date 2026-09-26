# -*- coding: utf-8 -*-
"""Figure 5 -- position estimators for the TCP (MLE, LMS, mLMS), Monte Carlo.

L = 50 nm, N = 100 photons, fwhm = 300 nm, seed 42 (scripts/_paperconfig.py).

(a) bias along x and (b) sigma/CRB versus the true position x0 on the x axis (0..L, i.e. inside
    and outside the TCP circle of radius L/2) for MLE, LMS (Balzarotti2017 Eq. S49-S50) and
    mLMS (Eq. S51) with SBR = 10.
(c) MLE without background at the TCP centre: sigma / CRB_lim (r -> 0 limit) and sigma / S27
    (point value) versus N -- the MLE is superefficient (model not regular, p_centre ~ r^2);
    LMS without background shown for reference (sigma_LMS = S27 exactly).
(d) MLE without background, N = 100: bias along x versus x0 near the centre (bias toward the
    centre; r = (2, 0) is the point quoted in the text).

Only the public donutloc API is used. Heavy Monte Carlo is cached in data/mc/fig5_*.npz;
``--quick`` uses a separate reduced cache, ``--no-cache`` recomputes.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _paperconfig as C  # noqa: E402
import _paperstyle as S  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

from donutloc import beams, closed_forms, estimators, fisher, montecarlo, patterns, photons  # noqa: E402

L, N, FWHM, SBR = C.L_REF, C.N_REF, C.FWHM, C.SBR_MLE
# MLE search disk (_paperconfig.MLE_RADIUS_*_OVER_L): radius 2L for the x-sweep (true positions reach x0 = L); the centre points use
# radius L, the r1 convention for mle_efficiency_center (worker B, verified); the same point must
# be computed with the same radius and seed in compute_paper_numbers.py to give the same value.
R_SWEEP = C.MLE_RADIUS_SWEEP_OVER_L * L
R_CENTRE = C.MLE_RADIUS_CENTRE_OVER_L * L


def sigma_se(err):
    """Bootstrap SE of sigma = sqrt((var_x+var_y)/2) if available (W1, r3), else Gaussian."""
    f = getattr(montecarlo, "bootstrap_sigma_se", None)
    if f is not None:
        return float(f(err, n_boot=C.N_BOOT, seed=C.SEED)), 1
    s = np.sqrt(0.5 * err.var(axis=0, ddof=1).sum())
    return float(s / (2.0 * np.sqrt(err.shape[0]))), 0


def stats(est, r_true):
    e = np.asarray(est, float) - np.asarray(r_true, float)
    R = e.shape[0]
    var = e.var(axis=0, ddof=1)
    se, boot = sigma_se(e)
    return dict(bias=e.mean(axis=0), bias_se=np.sqrt(var / R),
                sigma=float(np.sqrt(0.5 * var.sum())), sigma_se=se, boot=boot)


def model(sbr):
    beam = beams.make_beam("donut", fwhm=FWHM)
    return photons.make_model(patterns.tcp_centers(L), beam, sbr=sbr)


def run_sweep(quick):
    """(a,b): x-sweep with SBR=10, the three estimators on the SAME counts (run_mc seed 42)."""
    R = 1000 if quick else C.N_REP_MLE
    xs = np.linspace(0.0, L, 11)
    p_fn = model(SBR)
    est_fns = {
        "mle": lambda c: estimators.mle(c, p_fn, search_radius=R_SWEEP),
        "lms": lambda c: estimators.lms_tcp(c, L, FWHM, sbr=SBR),
        "mlms": lambda c: estimators.mlms_tcp(c, L, FWHM, sbr=SBR),
    }
    out = {"x": xs, "n_rep": R}
    for k in est_fns:
        for f in ("bias", "bias_se"):
            out["%s_%s" % (k, f)] = np.zeros((xs.size, 2))
        for f in ("sigma", "sigma_se"):
            out["%s_%s" % (k, f)] = np.zeros(xs.size)
    out["crb"] = np.zeros(xs.size)
    boot = 1
    for i, x in enumerate(xs):
        r = np.array([x, 0.0])
        out["crb"][i] = float(fisher.crb(p_fn, r, N))
        for k, fn in est_fns.items():
            mc = montecarlo.run_mc(fn, p_fn, r, N, R, seed=C.SEED)
            st = stats(mc["estimates"], r)
            boot = min(boot, st["boot"])
            for f in ("bias", "bias_se", "sigma", "sigma_se"):
                out["%s_%s" % (k, f)][i] = st[f]
    out["boot"] = boot
    return out


def run_centre(quick):
    """Centre point with SBR=10 and N_REP_MLE (= mle_efficiency_center) and the no-background
    N-scan of the MLE and LMS at the centre."""
    R = 1000 if quick else C.N_REP_MLE
    Ns = np.array(C.N_FIT)
    p_bg, p_0 = model(SBR), model(None)
    out = {"n_rep": R, "N": Ns}
    mc = montecarlo.run_mc(lambda c: estimators.mle(c, p_bg, search_radius=R_CENTRE), p_bg,
                           np.zeros(2), N, R, seed=C.SEED)
    st = stats(mc["estimates"], np.zeros(2))
    out["sbr_sigma"], out["sbr_sigma_se"] = st["sigma"], st["sigma_se"]
    out["sbr_crb"] = float(fisher.crb_limit(p_bg, N))
    boot = st["boot"]
    for k in ("mle", "lms"):
        out[k + "_sigma"] = np.zeros(Ns.size)
        out[k + "_sigma_se"] = np.zeros(Ns.size)
    out["crb_lim"] = np.zeros(Ns.size)
    out["s27"] = np.zeros(Ns.size)
    for i, n in enumerate(Ns):
        out["crb_lim"][i] = float(fisher.crb_limit(p_0, int(n)))
        out["s27"][i] = float(closed_forms.crb_tcp_center_point(L, int(n), FWHM))
        fns = {"mle": lambda c: estimators.mle(c, p_0, search_radius=R_CENTRE),
               "lms": lambda c: estimators.lms_tcp(c, L, FWHM)}
        for k, fn in fns.items():
            mc = montecarlo.run_mc(fn, p_0, np.zeros(2), int(n), R, seed=C.SEED)
            st = stats(mc["estimates"], np.zeros(2))
            boot = min(boot, st["boot"])
            out[k + "_sigma"][i], out[k + "_sigma_se"][i] = st["sigma"], st["sigma_se"]
    out["boot"] = boot
    return out


def run_bias0(quick):
    """(d): MLE without background, N=100, bias along x near the centre (many repetitions so
    that SE(bias) <= 0.01 nm at x0 = 2)."""
    R = C.N_REP_BIAS0 // 10 if quick else C.N_REP_BIAS0
    xs = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0])
    p_0 = model(None)
    out = {"x": xs, "n_rep": R, "bias": np.zeros((xs.size, 2)), "bias_se": np.zeros((xs.size, 2)),
           "sigma": np.zeros(xs.size), "sigma_se": np.zeros(xs.size), "crb": np.zeros(xs.size)}
    boot = 1
    for i, x in enumerate(xs):
        r = np.array([x, 0.0])
        mc = montecarlo.run_mc(lambda c: estimators.mle(c, p_0, search_radius=R_CENTRE), p_0, r,
                               N, R, seed=C.SEED)
        st = stats(mc["estimates"], r)
        boot = min(boot, st["boot"])
        for f in ("bias", "bias_se", "sigma", "sigma_se"):
            out[f][i] = st[f]
        out["crb"][i] = float(fisher.crb(p_0, r, N))
    out["boot"] = boot
    return out


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    t0 = time.time()
    sw = S.cached("fig5_sweep", lambda: run_sweep(args.quick), args.quick, args.no_cache)
    ce = S.cached("fig5_centre", lambda: run_centre(args.quick), args.quick, args.no_cache)
    b0 = S.cached("fig5_bias0", lambda: run_bias0(args.quick), args.quick, args.no_cache)
    print("fig5: Monte Carlo ready in %.1f s (bootstrap SE: %s)"
          % (time.time() - t0, "yes" if min(sw["boot"], ce["boot"], b0["boot"]) else
             "NO -- Gaussian sigma/(2 sqrt R) fallback"))

    # ------------------------------------------------------------------ numbers
    eff = ce["sbr_sigma"] / ce["sbr_crb"]
    S.report("mle_efficiency_center", eff, "", se=ce["sbr_sigma_se"] / ce["sbr_crb"])
    S.report("fig5_mle_sigma_center_sbr10_nm", ce["sbr_sigma"], "nm", se=ce["sbr_sigma_se"])
    S.report("fig5_crb_center_sbr10_nm", ce["sbr_crb"], "nm")
    S.report("fig5_n_rep_center", int(ce["n_rep"]))
    for i, n in enumerate(ce["N"]):
        S.report("fig5_mle_nobg_sigma_over_crblim_N%d" % n, ce["mle_sigma"][i] / ce["crb_lim"][i],
                 "", se=ce["mle_sigma_se"][i] / ce["crb_lim"][i])
        S.report("fig5_mle_nobg_sigma_over_s27_N%d" % n, ce["mle_sigma"][i] / ce["s27"][i],
                 "", se=ce["mle_sigma_se"][i] / ce["s27"][i])
        S.report("fig5_lms_nobg_sigma_over_s27_N%d" % n, ce["lms_sigma"][i] / ce["s27"][i],
                 "", se=ce["lms_sigma_se"][i] / ce["s27"][i])
    j2 = int(np.argmin(np.abs(b0["x"] - 2.0)))
    S.report("mle_nobg_bias_x_r2_nm", b0["bias"][j2, 0], "nm", se=b0["bias_se"][j2, 0])
    S.report("fig5_n_rep_bias0", int(b0["n_rep"]))
    for k in ("mle", "lms", "mlms"):
        for i, x in enumerate(sw["x"]):
            if x in (0.0, 25.0, 50.0):
                S.report("fig5_%s_sbr10_bias_x_x%d_nm" % (k, x), sw[k + "_bias"][i, 0], "nm",
                         se=sw[k + "_bias_se"][i, 0])
                S.report("fig5_%s_sbr10_sigma_over_crb_x%d" % (k, x),
                         sw[k + "_sigma"][i] / sw["crb"][i], "",
                         se=sw[k + "_sigma_se"][i] / sw["crb"][i])
    S.report("fig5_mle_sbr10_eff_max_dev_inside", float(np.max(np.abs(
        sw["mle_sigma"][sw["x"] <= L / 2] / sw["crb"][sw["x"] <= L / 2] - 1))), "")

    # ------------------------------------------------------------------ figure
    fig, axs = S.figure("double", aspect=0.62, nrows=2, ncols=2)
    ax = axs.ravel()
    lab = {"mle": "MLE", "lms": "LMS", "mlms": "mLMS"}
    mk = {"mle": "o", "lms": "s", "mlms": "^"}
    x = sw["x"]
    for a in ax[:2]:
        a.axvspan(0, L / 2, color="0.93", zorder=0, lw=0)
    ax[1].text(L / 4, 0.03, "inside TCP\n($|x_0|<L/2$)", transform=ax[1].get_xaxis_transform(),
               ha="center", va="bottom", fontsize=6.5, color="0.35")
    for k in ("mle", "lms", "mlms"):
        ax[0].errorbar(x, sw[k + "_bias"][:, 0], yerr=sw[k + "_bias_se"][:, 0], marker=mk[k],
                       color=S.COLORS[k], label=lab[k], ms=3.5, capsize=0)
        ax[1].errorbar(x, sw[k + "_sigma"] / sw["crb"], yerr=sw[k + "_sigma_se"] / sw["crb"],
                       marker=mk[k], color=S.COLORS[k], label=lab[k], ms=3.5, capsize=0)
    ax[0].axhline(0, color="0.5", lw=0.6, ls=":")
    ax[0].set_xlabel(r"true position $x_0$ (nm), $y_0=0$")
    ax[0].set_ylabel(r"bias $\langle\hat{x}\rangle - x_0$ (nm)")
    ax[0].legend(loc="lower left", title="SBR = %g" % SBR, title_fontsize=7)
    ax[1].axhline(1, color=S.COLORS["crb"], lw=0.8, ls="--")
    ax[1].text(L * 0.7, 0.97, "CRB", ha="center", va="top", fontsize=7)
    ax[1].set_yscale("log")
    ax[1].set_yticks([0.1, 0.2, 0.5, 1.0, 2.0])
    ax[1].yaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
    ax[1].yaxis.set_minor_formatter(mticker.NullFormatter())
    ax[1].set_xlabel(r"true position $x_0$ (nm), $y_0=0$")
    ax[1].set_ylabel(r"$\sigma\,/\,\sigma_\mathrm{CRB}$")
    for a in ax[:2]:
        a.set_xlim(-1, L + 1)

    Ns = ce["N"]
    ax[2].errorbar(Ns, ce["mle_sigma"] / ce["crb_lim"], yerr=ce["mle_sigma_se"] / ce["crb_lim"],
                   marker="o", color=S.COLORS["mle"], label=r"MLE, $\sigma/\sigma_\mathrm{lim}$",
                   capsize=0)
    ax[2].errorbar(Ns, ce["mle_sigma"] / ce["s27"], yerr=ce["mle_sigma_se"] / ce["s27"],
                   marker="o", mfc="white", color=S.COLORS["mle"], ls="--",
                   label=r"MLE, $\sigma/\sigma_\mathrm{S27}$", capsize=0)
    ax[2].errorbar(Ns, ce["lms_sigma"] / ce["crb_lim"], yerr=ce["lms_sigma_se"] / ce["crb_lim"],
                   marker="s", color=S.COLORS["lms"], label=r"LMS, $\sigma/\sigma_\mathrm{lim}$",
                   capsize=0)
    ax[2].axhline(1, color=S.COLORS["crb"], lw=0.8, ls="--")
    ax[2].axhline(float(np.mean(ce["s27"] / ce["crb_lim"])), color="0.5", lw=0.8, ls=":")
    ax[2].text(Ns[0], 1.0, r"  $\sigma_\mathrm{lim}$ ($r\to0$)", va="bottom", fontsize=6.5)
    ax[2].text(Ns[0], float(np.mean(ce["s27"] / ce["crb_lim"])), r"  S27 (point value)",
               va="bottom", fontsize=6.5, color="0.35")
    ax[2].set_xscale("log")
    ax[2].set_ylim(0.6, 1.2)
    ax[2].set_xlabel(r"photons $N$")
    ax[2].set_ylabel(r"$\sigma$ / bound  (centre, no bkg)")
    ax[2].legend(loc="lower center", ncol=3, fontsize=6.5, columnspacing=1.0, handlelength=1.4)

    ax[3].errorbar(b0["x"], b0["bias"][:, 0], yerr=b0["bias_se"][:, 0], marker="o",
                   color=S.COLORS["mle"], capsize=0, label="MLE, no bkg, $N=%d$" % N)
    ax[3].axhline(0, color="0.5", lw=0.6, ls=":")
    ax[3].annotate("%.2f nm" % b0["bias"][j2, 0], (b0["x"][j2], b0["bias"][j2, 0]),
                   xytext=(4.6, b0["bias"][j2, 0] - 0.03), textcoords="data", va="center", fontsize=6.5,
                   arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3"))
    ax[3].set_xlabel(r"true position $x_0$ (nm), $y_0=0$")
    ax[3].set_ylabel(r"bias $\langle\hat{x}\rangle - x_0$ (nm)")
    ax[3].set_ylim(float(b0["bias"][:, 0].min()) - 0.08, None)
    ax[3].legend(loc="upper left")
    for a, l in zip(ax, "abcd"):
        S.panel_label(a, l)
    S.savefig(fig, "fig5_estimators")
    S.write_summary("fig5")
    print("fig5: total %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
