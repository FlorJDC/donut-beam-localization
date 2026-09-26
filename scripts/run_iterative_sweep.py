# -*- coding: utf-8 -*-
"""Iterative MINFLUX vs. total photon budget (OBJECTIVE R4, ledger item W11).

Protocol (same as W10/W11 of the ledger): ``experiments.iterative_minflux`` with
``_paperconfig.ITER`` (4 iterations, L 150 -> 25 nm geometric, equal photon split, re-centring,
MLE), no background, seed 42 (the same seed for every N_total), ``ITER_N_LIST`` photon budgets and
``ITER_N_REP`` repetitions per budget.

Outputs ``data/iterative_sweep.json`` (``--quick``: 1000 repetitions, written to
``data/mc/iterative_sweep_quick.json`` so that the final file is never overwritten by a quick
run).  Per N: sigma, bootstrap SE (``N_BOOT`` resamples of the repetitions), Gaussian SE, rmse,
|bias|, per-iteration sigma and L_k, ideal camera sigma_PSF/sqrt(N), ratio.  Global: log-log
slope of sigma vs N over the full range and for N >= 500 with bootstrap 95 % CIs (independent
resampling of the repetitions of every N), ratio range, and at N = ITER_N_TOTAL the variants
SBR = 10, rule = "adaptive" and recenter = False (the latter flagged as an artefact of the
search disk, not physics).

Usage:  python scripts/run_iterative_sweep.py [--quick] [--out PATH]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _paperconfig as C  # noqa: E402

from donutloc import experiments, montecarlo  # noqa: E402

FINAL_PATH = os.path.join(C.DATA, "iterative_sweep.json")
QUICK_PATH = os.path.join(C.MC, "iterative_sweep_quick.json")


def _fit_slope(N, sig):
    return float(np.polyfit(np.log(np.asarray(N, float)), np.log(np.asarray(sig, float)), 1)[0])


def _run(N_total, n_rep, **over):
    kw = dict(C.ITER)
    kw.update(over)
    return experiments.iterative_minflux(int(N_total), n_rep=n_rep, seed=C.SEED,
                                         sigma_psf=C.SIGMA_PSF, **kw)


def _summ(res, n_boot, boot_seed):
    err = res["estimates"] - res["r_true"]
    se, samples = montecarlo.bootstrap_sigma_se(err, n_boot=n_boot, seed=boot_seed,
                                                return_samples=True)
    rad = np.hypot(err[:, 0], err[:, 1])
    d = {
        "sigma": float(res["sigma"]),
        "sigma_se": se,
        "sigma_se_gauss": float(res["sigma_se_gauss"]),
        "rmse": float(res["rmse"]),
        "bias_abs": float(res["bias_abs"]),
        "camera_sigma": float(res["camera_sigma"]),
        "ratio_to_camera": float(res["ratio_to_camera"]),
        "ratio_se": se / float(res["camera_sigma"]),
        "L": [float(v) for v in res["L"]],
        "N_k": [int(v) for v in res["N_k"]],
        "sigma_iter": [float(v) for v in res["sigma_iter"]],
        "sigma_se_iter": [float(v) for v in res["sigma_se_iter"]],
        "crb_center_iter": [float(v) for v in res["crb_center_iter"]],
        "crb_all_photons_Lmin": float(res["crb_all_photons_Lmin"]),
        "frac_err_gt_5nm": float(np.mean(rad > 5.0)),
        "max_err_nm": float(rad.max()),
    }
    return d, samples


def run_sweep(quick=False, n_rep=None, verbose=True):
    """Run the full sweep and return the JSON-able dict (does not write it)."""
    t0 = time.time()
    n_rep = (1000 if quick else C.ITER_N_REP) if n_rep is None else int(n_rep)
    n_boot = C.N_BOOT
    N_list = list(C.ITER_N_LIST)
    per_N, samples = [], []
    for i, N in enumerate(N_list):
        t = time.time()
        d, smp = _summ(_run(N, n_rep), n_boot, C.SEED + 1 + i)
        per_N.append(d)
        samples.append(smp)
        if verbose:
            print("N=%5d sigma=%.4f +- %.4f (gauss %.4f) ratio=%.4f  [%.1f s]"
                  % (N, d["sigma"], d["sigma_se"], d["sigma_se_gauss"], d["ratio_to_camera"],
                     time.time() - t), flush=True)
    N_arr = np.array(N_list, float)
    sig = np.array([d["sigma"] for d in per_N])
    S = np.stack(samples, axis=1)                          # (n_boot, n_N) bootstrap sigmas
    lN = np.log(N_arr)

    def slope_block(mask):
        val = _fit_slope(N_arr[mask], sig[mask])
        A = np.vstack([lN[mask], np.ones(mask.sum())]).T
        coef = np.linalg.lstsq(A, np.log(S[:, mask]).T, rcond=None)[0][0]
        lo, hi = np.percentile(coef, [2.5, 97.5])
        return {"value": val, "se": float(coef.std(ddof=1)), "ci95": [float(lo), float(hi)],
                "N_used": [int(n) for n in N_arr[mask]],
                "compatible_with_minus_half": bool(lo <= -0.5 <= hi)}

    slope_all = slope_block(np.ones(N_arr.size, bool))
    slope_500 = slope_block(N_arr >= 500)
    ratios = np.array([d["ratio_to_camera"] for d in per_N])

    extras = {}
    variants = [("sbr10", dict(sbr=10.0), False),
                ("adaptive", dict(rule="adaptive"), False),
                ("no_recenter", dict(recenter=False), True)]
    for j, (name, over, artefact) in enumerate(variants):
        t = time.time()
        d, _ = _summ(_run(C.ITER_N_TOTAL, n_rep, **over), n_boot, C.SEED + 101 + j)
        d["overrides"] = {k: v for k, v in over.items()}
        d["artefact"] = artefact
        if artefact:
            d["note"] = ("recenter=False: every TCP stays on the origin and the MLE search disk "
                         "(radius 0.75 L_k) is smaller than the emitter spread (L0/4 = 37.5 nm) "
                         "for the last L; the large final sigma is an artefact of the search-disk "
                         "truncation, not a physical effect of the zoom")
        extras[name] = d
        if verbose:
            print("N=%d %-12s sigma=%.4f +- %.4f  L=%s  [%.1f s]"
                  % (C.ITER_N_TOTAL, name, d["sigma"], d["sigma_se"],
                     np.round(d["L"], 2), time.time() - t), flush=True)

    out = {
        "script": "scripts/run_iterative_sweep.py",
        "quick": bool(quick),
        "params": {"iter": dict(C.ITER), "n_rep": n_rep, "seed": C.SEED, "n_boot": n_boot,
                   "sbr": None, "sigma_psf": C.SIGMA_PSF, "fwhm": C.FWHM,
                   "N_list": N_list, "N_ref": C.ITER_N_TOTAL,
                   "r0_spread": "L0/4", "search_radius": "0.75 L_k",
                   "boot_seeds": "SEED+1+i per N (i = index in N_list); SEED+101+j extras"},
        "N": N_list,
        "sigma": [d["sigma"] for d in per_N],
        "sigma_se": [d["sigma_se"] for d in per_N],
        "sigma_se_gauss": [d["sigma_se_gauss"] for d in per_N],
        "rmse": [d["rmse"] for d in per_N],
        "camera_sigma": [d["camera_sigma"] for d in per_N],
        "ratio_to_camera": [d["ratio_to_camera"] for d in per_N],
        "per_N": per_N,
        "slope_all": slope_all,
        "slope_N_ge_500": slope_500,
        "ratio_range": [float(ratios.min()), float(ratios.max())],
        "extras_N_ref": extras,
        "runtime_s": float(time.time() - t0),
    }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="1000 repetitions (separate file)")
    ap.add_argument("--out", default=None, help="output path (default: see module doc)")
    a = ap.parse_args(argv)
    out = run_sweep(quick=a.quick)
    path = a.out or (QUICK_PATH if a.quick else FINAL_PATH)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write("\n")
    s, s5 = out["slope_all"], out["slope_N_ge_500"]
    print("slope all   = %.4f  (95%% CI %.4f .. %.4f)" % (s["value"], *s["ci95"]))
    print("slope N>=500= %.4f  (95%% CI %.4f .. %.4f)" % (s5["value"], *s5["ci95"]))
    print("ratio range = %.4f .. %.4f" % tuple(out["ratio_range"]))
    print("wrote %s (%.0f s)" % (path, out["runtime_s"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
