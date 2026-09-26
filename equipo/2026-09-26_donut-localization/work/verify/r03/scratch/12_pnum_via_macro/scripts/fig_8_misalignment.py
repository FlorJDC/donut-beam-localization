# -*- coding: utf-8 -*-
"""Figure 8 -- TCP misalignment: naive versus honest MLE.

experiments.misalignment_study with _paperconfig.MIS (L = 100 nm, N = 500, SBR = 10,
n_patterns = 50 random misaligned TCPs x n_rep = 200 repetitions, seed 42): every one of the
4 zeros is displaced by delta in an independent random direction; the photons are simulated with
the true (misaligned) TCP and estimated by the MLE with the true model ("honest") and with the
ideal TCP ("naive").  True positions: the TCP centre and (L/4, 0).  The deltas include
MIS_DELTAS; the study re-creates its random generator for every delta (common random numbers),
so the value at each delta does not depend on which other deltas are computed.

(a) mean |bias| versus delta, with the Monte Carlo noise floor of an unbiased estimator
    (sigma sqrt(pi/(2R))) and the fitted line |bias| ~ slope * delta of the naive MLE (~0.75).
(b) sigma versus delta and the honest CRB (r -> r_true limit, averaged over patterns).
(c) rmse versus delta.

Error bars: bias_abs_se and sigma_se, both the standard error over the 50 patterns (the
between-pattern sigma_se of experiments.misalignment_study, r3 fix); rmse has no error bar.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _paperconfig as C  # noqa: E402
import _paperstyle as S  # noqa: E402

import numpy as np  # noqa: E402

from donutloc import experiments  # noqa: E402

DELTAS = tuple(sorted(set(C.MIS_DELTAS) | {1.0, 3.5, 7.5, 15.0}))
KEYS = ("bias_abs_mean", "bias_abs_se", "bias_chi2", "sigma_mean", "sigma_se", "sigma_se_within",
        "rmse")


def compute(quick):
    mis = dict(C.MIS)
    deltas = DELTAS
    if quick:
        mis.update(n_patterns=10, n_rep=100)
        deltas = tuple(C.MIS_DELTAS)
    r = experiments.misalignment_study(deltas, seed=C.SEED, **mis)
    out = {"delta": r["delta"], "positions": r["positions"], "crb_honest": r["crb_honest"],
           "floor": r["bias_noise_floor"], "identical": r["identical"].astype(int),
           "n_patterns": mis["n_patterns"], "n_rep": mis["n_rep"]}
    for m in ("honest", "naive"):
        for k in KEYS:
            out["%s_%s" % (m, k)] = r[m][k]
    return out


def fit_slope(delta, b, se, dmin=2.0):
    """Weighted least squares through the origin, |bias| = s * delta, for delta >= dmin."""
    m = delta >= dmin
    w = 1.0 / se[m] ** 2
    s = np.sum(w * delta[m] * b[m]) / np.sum(w * delta[m] ** 2)
    return float(s), float(1.0 / np.sqrt(np.sum(w * delta[m] ** 2)))


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    t0 = time.time()
    d = S.cached("fig8_misalignment", lambda: compute(args.quick), args.quick, args.no_cache)
    dl = d["delta"]
    P, R = int(d["n_patterns"]), int(d["n_rep"])
    print("fig8: n_patterns = %d, n_rep = %d, deltas = %s" % (P, R, list(dl)))

    # ------------------------------------------------------------------ numbers
    pos_tag = ["center", "Lq"]
    slopes = []
    canon = np.isin(dl, np.array(C.MIS_DELTAS))
    for q, tag in enumerate(pos_tag):
        s, se = fit_slope(dl, d["naive_bias_abs_mean"][:, q], d["naive_bias_abs_se"][:, q])
        slopes.append(s)
        S.report("fig8_naive_bias_slope_all_deltas_%s" % tag, s, "", se=se)
        s, se = fit_slope(dl[canon], d["naive_bias_abs_mean"][canon, q],
                          d["naive_bias_abs_se"][canon, q])
        S.report("misalignment_naive_bias_slope_%s" % tag, s, "", se=se)
    for i, dd in enumerate(dl):
        dtag = ("%g" % dd).replace(".", "p")
        for q, tag in enumerate(pos_tag):
            for m in ("naive", "honest"):
                S.report("misalignment_%s_bias_abs_%s_d%s_nm" % (m, tag, dtag),
                         d[m + "_bias_abs_mean"][i, q], "nm", se=d[m + "_bias_abs_se"][i, q])
                S.report("misalignment_%s_sigma_%s_d%s_nm" % (m, tag, dtag),
                         d[m + "_sigma_mean"][i, q], "nm", se=d[m + "_sigma_se"][i, q])
                S.report("misalignment_%s_rmse_%s_d%s_nm" % (m, tag, dtag), d[m + "_rmse"][i, q],
                         "nm")
            S.report("misalignment_crb_honest_%s_d%s_nm" % (tag, dtag), d["crb_honest"][i, q], "nm")
            S.report("misalignment_mc_floor_%s_d%s_nm" % (tag, dtag), d["floor"][i, q], "nm")
    S.report("misalignment_delta0_identical", int(d["identical"][0]))
    S.report("fig8_n_patterns", P)
    S.report("fig8_n_rep", R)

    # ------------------------------------------------------------------ figure
    fig, ax = S.figure("double", aspect=0.34, nrows=1, ncols=3)
    mk = ["o", "s"]
    lab_pos = ["centre", r"$(L/4,0)$"]
    dd = np.linspace(0, dl.max() * 1.03, 50)
    a = ax[0]
    a.fill_between(dl, 0, d["floor"][:, 0], color="0.88", lw=0, label="MC noise floor")
    a.plot(dd, slopes[0] * dd, color=S.COLORS["naive"], lw=0.7, ls=":",
           label=r"$%.2f\,\delta$" % slopes[0])
    for q in range(2):
        for m in ("naive", "honest"):
            a.errorbar(dl, d[m + "_bias_abs_mean"][:, q], yerr=d[m + "_bias_abs_se"][:, q],
                       marker=mk[q], ms=3.3, capsize=0, color=S.COLORS[m],
                       mfc=S.COLORS[m] if q == 0 else "white", ls="-" if q == 0 else "--",
                       lw=0.9, label="%s, %s" % (m, lab_pos[q]))
    a.set_xlabel(r"zero displacement $\delta$ (nm)")
    a.set_ylabel(r"mean $|\mathrm{bias}|$ (nm)")
    a.set_ylim(bottom=0)
    a.legend(loc="upper left", fontsize=5.8, handlelength=1.8)

    for a, key, ylab in ((ax[1], "sigma_mean", r"$\sigma$ (nm)"), (ax[2], "rmse", r"rmse (nm)")):
        for q in range(2):
            for m in ("naive", "honest"):
                se = d[m + "_sigma_se"][:, q] if key == "sigma_mean" else None
                a.errorbar(dl, d["%s_%s" % (m, key)][:, q], yerr=se, marker=mk[q], ms=3.3,
                           capsize=0, color=S.COLORS[m], mfc=S.COLORS[m] if q == 0 else "white",
                           ls="-" if q == 0 else "--", lw=0.9, label="%s, %s" % (m, lab_pos[q]))
            a.plot(dl, d["crb_honest"][:, q], color=S.COLORS["crb"], lw=0.8,
                   ls="-" if q == 0 else "--", alpha=0.7,
                   label="honest CRB, %s" % lab_pos[q])
        a.set_xlabel(r"zero displacement $\delta$ (nm)")
        a.set_ylabel(ylab)
    ax[1].legend(loc="upper left", fontsize=5.8, handlelength=1.8)
    ax[2].text(0.04, 0.96, "$L=%g$ nm, $N=%d$, SBR $=%g$\n%d patterns $\\times$ %d reps"
               % (C.MIS["L"], C.MIS["N"], C.MIS["sbr"], P, R), transform=ax[2].transAxes,
               ha="left", va="top", fontsize=6.2)
    for axx, l in zip(ax, "abc"):
        S.panel_label(axx, l)
    S.savefig(fig, "fig8_misalignment")
    S.write_summary("fig8")
    print("fig8: total %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
