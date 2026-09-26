# -*- coding: utf-8 -*-
"""Figure 6 -- iterative MINFLUX versus an ideal camera.

Reads data/iterative_sweep.json (scripts/run_iterative_sweep.py: experiments.iterative_minflux
with _paperconfig.ITER -- 4 iterations, L = 150 -> 25 nm geometric, equal photon split,
re-centring on the previous estimate, MLE in a disk of radius 0.75 L_k, no background, seed 42,
ITER_N_REP = 10000 emitters uniform in a disk of radius L0/4 per N_total; SE by bootstrap).

(a) final sigma +- SE versus total photons N_total, with the ideal camera sigma_PSF/sqrt(N)
    (sigma_PSF = 100 nm), the fitted log-log slope (N >= 500) and the centre CRB if all photons
    were spent at L = 25 nm; the N_total = 1000 variants with SBR = 10 and with the adaptive
    L rule; the run without re-centring in grey, labelled as an artefact of the search disk
    (not physics).
(b) ratio sigma / sigma_camera versus N_total.
(c) sigma after each iteration versus L_k (N_total = 1000), with the centre CRB of each
    iteration.

If data/iterative_sweep.json is missing, a reduced version (1000 repetitions) is computed with
the same protocol and a warning is printed; it is cached in data/mc/fig6_fallback*.npz.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _paperconfig as C  # noqa: E402
import _paperstyle as S  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

from donutloc import experiments  # noqa: E402

SWEEP = os.path.join(C.DATA, "iterative_sweep.json")
_RUN_KEYS = ("sigma", "sigma_se", "rmse", "camera_sigma", "ratio_to_camera", "L", "N_k",
             "sigma_iter", "sigma_se_iter", "crb_center_iter", "crb_all_photons_Lmin")


def _run(N, **over):
    kw = dict(C.ITER)
    kw.update(over)
    out = experiments.iterative_minflux(int(N), n_rep=1000, seed=C.SEED,
                                        sigma_psf=C.SIGMA_PSF, **kw)
    return {k: out[k] for k in _RUN_KEYS if k in out}


def fallback(quick):
    """Same protocol as run_iterative_sweep.py but 1000 repetitions, no bootstrap CI."""
    Ns = np.array(C.ITER_N_LIST)
    runs = [_run(n) for n in Ns]
    out = {"N": Ns}
    for k in ("sigma", "sigma_se", "rmse", "camera_sigma", "ratio_to_camera",
              "crb_all_photons_Lmin"):
        out[k] = np.array([r[k] for r in runs])
    j = list(Ns).index(C.ITER_N_TOTAL)
    for k in ("L", "sigma_iter", "sigma_se_iter", "crb_center_iter"):
        out["ref_" + k] = np.asarray(runs[j][k])
    for name, over in (("sbr10", dict(sbr=10.0)), ("adaptive", dict(rule="adaptive")),
                       ("no_recenter", dict(recenter=False))):
        r = _run(C.ITER_N_TOTAL, **over)
        for k in ("sigma", "sigma_se", "L", "sigma_iter", "crb_center_iter"):
            out["%s_%s" % (name, k)] = np.asarray(r[k])
    return out


def load(args):
    if os.path.exists(SWEEP):
        with open(SWEEP, encoding="utf-8") as f:
            d = json.load(f)
        per = d["per_N"]
        Ns = np.asarray(d["N"])
        j = list(Ns).index(C.ITER_N_TOTAL)
        out = {"N": Ns, "source": "data/iterative_sweep.json", "quick": bool(d.get("quick")),
               "n_rep": d["params"]["n_rep"]}
        for k in ("sigma", "sigma_se", "rmse", "camera_sigma", "ratio_to_camera"):
            out[k] = np.asarray(d[k], float)
        out["ratio_se"] = np.array([p.get("ratio_se", np.nan) for p in per])
        out["crb_all_photons_Lmin"] = np.array([p["crb_all_photons_Lmin"] for p in per])
        for k in ("L", "sigma_iter", "sigma_se_iter", "crb_center_iter"):
            out["ref_" + k] = np.asarray(per[j][k], float)
        for name in ("sbr10", "adaptive", "no_recenter"):
            e = d["extras_N_ref"][name]
            for k in ("sigma", "sigma_se", "L", "sigma_iter", "crb_center_iter", "ratio_to_camera"):
                out["%s_%s" % (name, k)] = np.asarray(e[k], float)
        for s in ("slope_all", "slope_N_ge_500"):
            out[s] = d[s]["value"]
            out[s + "_ci"] = np.asarray(d[s]["ci95"], float)
            out[s + "_N"] = np.asarray(d[s]["N_used"])
        out["ratio_range"] = np.asarray(d["ratio_range"], float)
        return out
    print("WARNING: %s not found -- computing a reduced fallback (1000 repetitions) with the "
          "same protocol; rerun after scripts/run_iterative_sweep.py" % SWEEP)
    out = S.cached("fig6_fallback", lambda: fallback(args.quick), args.quick, args.no_cache)
    out = dict(out)
    out.update(source="fallback (1000 reps)", quick=True, n_rep=1000)
    out["ratio_se"] = out["sigma_se"] / out["camera_sigma"]
    for name in ("sbr10", "adaptive", "no_recenter"):
        out[name + "_ratio_to_camera"] = out[name + "_sigma"] / out["camera_sigma"][
            list(out["N"]).index(C.ITER_N_TOTAL)]
    lN = np.log(out["N"].astype(float))
    for s, m in (("slope_all", np.ones(lN.size, bool)), ("slope_N_ge_500", out["N"] >= 500)):
        out[s] = float(np.polyfit(lN[m], np.log(out["sigma"][m]), 1)[0])
        out[s + "_ci"] = np.array([np.nan, np.nan])
        out[s + "_N"] = out["N"][m]
    out["ratio_range"] = np.array([out["ratio_to_camera"].min(), out["ratio_to_camera"].max()])
    return out


def main():
    args = S.parse_args(__doc__.splitlines()[1])
    t0 = time.time()
    d = load(args)
    print("fig6: source = %s (n_rep = %s)" % (d["source"], d["n_rep"]))
    Ns = d["N"]
    j = list(Ns).index(C.ITER_N_TOTAL)

    # ------------------------------------------------------------------ numbers
    S.report("iterative_sigma_nm", d["sigma"][j], "nm", se=d["sigma_se"][j])
    S.report("camera_sigma_nm", d["camera_sigma"][j], "nm")
    S.report("iterative_ratio_to_camera_N1000", d["ratio_to_camera"][j], "", se=d["ratio_se"][j])
    for i, n in enumerate(Ns):
        S.report("fig6_iterative_sigma_N%d_nm" % n, d["sigma"][i], "nm", se=d["sigma_se"][i])
        S.report("fig6_ratio_to_camera_N%d" % n, d["ratio_to_camera"][i], "", se=d["ratio_se"][i])
    for s in ("slope_all", "slope_N_ge_500"):
        S.report("iterative_%s" % s, d[s], "")
        S.report("iterative_%s_ci95_lo" % s, float(d[s + "_ci"][0]), "")
        S.report("iterative_%s_ci95_hi" % s, float(d[s + "_ci"][1]), "")
    S.report("iterative_ratio_min", float(d["ratio_range"][0]), "")
    S.report("iterative_ratio_max", float(d["ratio_range"][1]), "")
    S.report("iterative_sbr10_sigma_nm", float(d["sbr10_sigma"]), "nm", se=float(d["sbr10_sigma_se"]))
    S.report("iterative_adaptive_sigma_nm", float(d["adaptive_sigma"]), "nm",
             se=float(d["adaptive_sigma_se"]))
    S.report("iterative_no_recenter_sigma_nm_ARTEFACT", float(d["no_recenter_sigma"]), "nm",
             se=float(d["no_recenter_sigma_se"]))
    S.report("fig6_crb_all_photons_L25_N1000_nm", float(d["crb_all_photons_Lmin"][j]), "nm")
    S.report("fig6_iterative_over_crb_all_photons_N1000",
             float(d["sigma"][j] / d["crb_all_photons_Lmin"][j]), "")
    S.report("fig6_n_rep", int(d["n_rep"]))

    # ------------------------------------------------------------------ figure
    fig, ax = S.figure("double", aspect=0.36, nrows=1, ncols=3,
                       gridspec_kw={"width_ratios": [1.25, 1.0, 1.0]})
    Nd = np.geomspace(Ns.min() * 0.8, Ns.max() * 1.25, 50)
    a = ax[0]
    a.plot(Nd, C.SIGMA_PSF / np.sqrt(Nd), color=S.COLORS["cam"], lw=1.0,
           label=r"ideal camera $\sigma_\mathrm{PSF}/\sqrt{N}$")
    a.errorbar(Ns, d["sigma"], yerr=d["sigma_se"], marker="o", ls="none", color=S.COLORS["lg"],
               label="iterative MINFLUX", zorder=4)
    m = Ns >= 500
    lN = np.log(Ns[m].astype(float))
    b = np.mean(np.log(d["sigma"][m])) - d["slope_N_ge_500"] * np.mean(lN)
    a.plot(Nd, np.exp(b) * Nd ** d["slope_N_ge_500"], color=S.COLORS["lg"], lw=0.8,
           label=r"fit $N\geq500$: slope $%.3f$" % d["slope_N_ge_500"])
    a.plot(Ns, d["crb_all_photons_Lmin"], color="0.45", ls="--", lw=0.8,
           label=r"CRB, all $N$ at $L=25$ nm")
    a.plot(C.ITER_N_TOTAL, d["sbr10_sigma"], marker="s", ls="none", color=S.COLORS["sbr"],
           mec="k", mew=0.3, label="SBR = 10", zorder=5)
    a.plot(C.ITER_N_TOTAL, d["adaptive_sigma"], marker="D", ls="none", color=S.COLORS["lms"],
           mec="k", mew=0.3, ms=3.2, label="adaptive $L_k$", zorder=5)
    a.plot(C.ITER_N_TOTAL, d["no_recenter_sigma"], marker="x", ls="none", color="0.55", ms=4)
    a.annotate("no re-centring\n(search-disk artefact)", (C.ITER_N_TOTAL, d["no_recenter_sigma"]),
               xytext=(1500, 8.1), fontsize=6, color="0.45", va="center", ha="left",
               arrowprops=dict(arrowstyle="-", lw=0.5, color="0.55"))
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_ylim(0.08, 1000)       # head-room: the legend sits above the annotated artefact point
    a.set_xlim(180, 1.5e4)       # keeps the "10000" tick label inside the panel
    a.set_xlabel(r"total photons $N_\mathrm{total}$")
    a.set_ylabel(r"localization precision $\sigma$ (nm)")
    a.legend(loc="upper right", fontsize=5.8, handlelength=1.5, borderaxespad=0.3)
    a.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
    a.yaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))

    a = ax[1]
    a.errorbar(Ns, d["ratio_to_camera"], yerr=d["ratio_se"], marker="o", color=S.COLORS["lg"],
               capsize=0)
    a.set_xscale("log")
    a.set_xlabel(r"total photons $N_\mathrm{total}$")
    a.set_ylabel(r"$\sigma\,/\,\sigma_\mathrm{camera}$")
    a.set_ylim(0.14, 0.20)
    a.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))

    a = ax[2]
    k = np.arange(1, d["ref_L"].size + 1)
    a.errorbar(d["ref_L"], d["ref_sigma_iter"], yerr=d["ref_sigma_se_iter"], marker="o",
               color=S.COLORS["lg"], capsize=0, label="fixed, no bkg")
    a.plot(d["ref_L"], d["ref_crb_center_iter"], color=S.COLORS["lg"], ls="--", lw=0.8,
           label="centre CRB")
    a.plot(d["sbr10_L"], d["sbr10_sigma_iter"], marker="s", color=S.COLORS["sbr"], mec="k",
           mew=0.3, label="fixed, SBR = 10")
    for i in range(k.size):
        a.annotate("%d" % k[i], (d["ref_L"][i], d["ref_sigma_iter"][i]), xytext=((5 if i == 0 else -3), -9),
                   textcoords="offset points", fontsize=6, color="0.3", ha=("left" if i == 0 else "right"))
    a.set_xscale("log")
    a.set_yscale("log")
    a.invert_xaxis()
    a.set_xlabel(r"TCP diameter $L_k$ (nm)")
    a.set_ylabel(r"$\sigma$ after iteration $k$ (nm)")
    a.set_title(r"$N_\mathrm{total}=%d$ ($4\times250$)" % C.ITER_N_TOTAL, fontsize=7)
    a.legend(loc="upper right", fontsize=6)
    a.set_xticks([150, 100, 50, 25])
    a.set_yticks([0.5, 1, 2, 5])
    a.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
    a.xaxis.set_minor_formatter(mticker.NullFormatter())
    a.yaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
    for axx, l in zip(ax, "abc"):
        S.panel_label(axx, l)
    S.savefig(fig, "fig6_iterative")
    S.write_summary("fig6")
    print("fig6: total %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
