# -*- coding: utf-8 -*-
"""Figura 9 -- efectos de temporización propios de p-MINFLUX (donutloc.pminflux).

Genera ``informe/figures/fig_pminflux_timing.{png,pdf}`` y ``informe/data/pminflux_timing.json``.

TCP L = 100 nm, dona LG fwhm = 300 nm, sin fondo, semilla 42.

(a) Flickering (telegrafía on/off, t_on = t_off = 100 µs, arranque estacionario), dwell total del
    patrón 400 µs, N medio = 100.  Secuencial: 4 r bloques consecutivos (r = 1, 2, 5, 10, 25).
    Entrelazado (p-MINFLUX): pulsos cada 50 ns por exposición, desfasados 12.5 ns; se cuentan
    exactamente los pulsos que caen en intervalos "on".  Tres estimaciones de sigma_fl:
      * ``pop``: error del MLE sobre conteos esperados (sin ruido de Poisson), RMS por eje -- el
        límite N -> inf, sigma_fl "puro";
      * ``fixedN``: conteos multinomiales con N = 100 exacto; sigma_fl^2 = STD^2 - STD_ctrl^2,
        con STD_ctrl la del mismo MLE sin flickering (mismas N); también frente a sigma_CRB(100);
      * ``poisson``: conteos Poisson (N variable; se descartan N = 0); control con los mismos N_i.
    STD "por eje" = sqrt((var_x + var_y)/2), mismo criterio que el CRB del proyecto.
(b) Cross-talk por tiempo de vida entre ventanas TCSPC (T = 12.5 ns): razón CRB(M)/CRB(ideal)
    en el centro y promedio en r <= L/2, para tau = 1..5 ns y dos órdenes de pulsos.
(c) Sesgo del MLE que ignora M (N = 500) a lo largo de x en [0, L/2]; MLE con M en el modelo.
"""
import json
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import _paperstyle as S  # noqa: E402  (pone src/ en sys.path y fija rcParams)

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from donutloc import beams, estimators, fisher, patterns, photons, pminflux as pm  # noqa: E402

ROOT = os.path.dirname(_HERE)
FIGDIR = os.path.join(ROOT, "informe", "figures")
DATADIR = os.path.join(ROOT, "informe", "data")

SEED = 42
L = 100.0
FWHM = 300.0
# flickering
T_TOTAL_US = 400.0
T_ON_US = 100.0
T_OFF_US = 100.0
N_MEAN = 100
REPS = (1, 2, 5, 10, 25)
N_LOC = 6000
POS_FL = {"centro": (0.0, 0.0), "x20": (20.0, 0.0)}
# cross-talk
TAUS = (1.0, 2.0, 3.0, 4.0, 5.0)
ORDERS = {"A": (0, 1, 2, 3), "B": (0, 2, 1, 3)}
N_CT = 500
R_CT = 2000
XS = np.linspace(0.0, L / 2, 11)
DISK_STEP = 2.5
SEARCH = 0.75 * L


def _axis_var(err):
    """Per-axis variance (var_x + var_y)/2 and its standard error."""
    e = err - err.mean(axis=0)
    d2 = 0.5 * (e[:, 0] ** 2 + e[:, 1] ** 2)
    return float(d2.mean()), float(d2.std(ddof=1) / np.sqrt(len(d2)))


def _sig_fl(v, v_ref, se, se_ref):
    d = v - v_ref
    s = float(np.sqrt(max(d, 0.0)))
    se_d = float(np.sqrt(se ** 2 + se_ref ** 2))
    return s, d, se_d


def flicker(rng, model, crb1):
    out = {}
    schemes = [("seq_r%d" % r, "sequential", r) for r in REPS] + [("interleaved", "interleaved", 1)]
    for pname, pos in POS_FL.items():
        pos = np.array(pos)
        p = model(pos)
        crbN = float(crb1[pname]) / np.sqrt(N_MEAN)
        res = {"sigma_crb_N100_nm": crbN}
        for key, scheme, reps in schemes:
            w = pm.exposure_weights(scheme, N_LOC, t_total=T_TOTAL_US, t_on=T_ON_US,
                                    t_off=T_OFF_US, reps=reps, rng=rng)
            on = (w * p[None, :]).sum(axis=1) > 0     # some signal expected
            d = {"n_loc": int(N_LOC), "n_no_signal": int((~on).sum())}
            # relative spread of the weights among the 4 exposures (localizations with >= 10 µs on)
            g = w.sum(axis=1) > 4 * 10.0 / (0.5 * T_TOTAL_US)
            rel = w[g] / w[g].mean(axis=1, keepdims=True)
            d["weight_rel_std"] = float(rel.std())
            d["weight_rel_maxdev"] = float(np.abs(rel - 1).max())
            # (i) population (noise-free) error
            lam = w[on] * p[None, :]
            est = estimators.mle(1e6 * lam / lam.sum(axis=1, keepdims=True), model, SEARCH)
            e = est - pos
            d["pop_sigma_fl_nm"] = float(np.sqrt(0.5 * np.mean(np.sum(e ** 2, axis=1))))
            d["pop_bias_nm"] = e.mean(axis=0).tolist()
            # (ii) fixed N = 100
            c = pm.simulate_flicker_counts(p, w[on], n_mean=N_MEAN, rng=rng, fixed_N=N_MEAN)
            c0 = photons.sample_counts(p, N_MEAN, size=int(on.sum()), rng=rng)
            e = estimators.mle(c, model, SEARCH) - pos
            e0 = estimators.mle(c0, model, SEARCH) - pos
            v, se = _axis_var(e)
            v0, se0 = _axis_var(e0)
            s, dd, sed = _sig_fl(v, v0, se, se0)
            s2, dd2, _ = _sig_fl(v, crbN ** 2, se, 0.0)
            d.update(fixedN_std_nm=float(np.sqrt(v)), fixedN_std_ctrl_nm=float(np.sqrt(v0)),
                     fixedN_sigma_fl_nm=s, fixedN_sigma_fl2_nm2=dd, fixedN_sigma_fl2_se_nm2=sed,
                     fixedN_sigma_fl_vs_crb_nm=s2, fixedN_bias_nm=e.mean(axis=0).tolist(),
                     fixedN_std_x_nm=float(e[:, 0].std(ddof=1)),
                     fixedN_std_y_nm=float(e[:, 1].std(ddof=1)))
            # (iii) Poisson counts, control with the same N_i
            c = pm.simulate_flicker_counts(p, w[on], n_mean=N_MEAN, rng=rng)
            Ni = c.sum(axis=1)
            k = Ni > 0
            c = c[k]
            Ni = Ni[k]
            c0 = np.array([rng.multinomial(int(n), p) for n in Ni])
            e = estimators.mle(c, model, SEARCH) - pos
            e0 = estimators.mle(c0, model, SEARCH) - pos
            v, se = _axis_var(e)
            v0, se0 = _axis_var(e0)
            s, dd, sed = _sig_fl(v, v0, se, se0)
            crb_eff = float(crb1[pname]) * np.sqrt(np.mean(1.0 / Ni))
            d.update(poisson_n_used=int(k.sum()), poisson_N_mean=float(Ni.mean()),
                     poisson_std_nm=float(np.sqrt(v)), poisson_std_ctrl_nm=float(np.sqrt(v0)),
                     poisson_sigma_fl_nm=s, poisson_sigma_fl2_nm2=dd,
                     poisson_sigma_fl2_se_nm2=sed, poisson_sigma_crb_eff_nm=crb_eff,
                     poisson_std_x_nm=float(e[:, 0].std(ddof=1)),
                     poisson_std_y_nm=float(e[:, 1].std(ddof=1)))
            res[key] = d
            print("flicker %-7s %-12s pop=%.3f fixedN: std=%.3f ctrl=%.3f crb=%.3f sfl=%.3f | "
                  "poisson: std=%.3f ctrl=%.3f sfl=%.3f" % (
                      pname, key, d["pop_sigma_fl_nm"], d["fixedN_std_nm"],
                      d["fixedN_std_ctrl_nm"], crbN, d["fixedN_sigma_fl_nm"],
                      d["poisson_std_nm"], d["poisson_std_ctrl_nm"], d["poisson_sigma_fl_nm"]))
        out[pname] = res
    return out


def crosstalk(rng, c, beam, ideal):
    g = np.arange(-L / 2, L / 2 + 1e-9, DISK_STEP)
    X, Y = np.meshgrid(g, g)
    pts = np.stack([X.ravel(), Y.ravel()], axis=-1)
    pts = pts[np.sum(pts ** 2, axis=1) <= (L / 2) ** 2 + 1e-9]
    crb_id = fisher.crb(ideal, pts, 1.0, zero_policy="limit")
    crb_id0 = float(fisher.crb(ideal, np.zeros(2), 1.0, zero_policy="limit"))
    crb_id0p = float(fisher.crb(ideal, np.zeros(2), 1.0, zero_policy="point"))
    out = {"disk_n_points": int(len(pts)), "crb1_ideal_centre_limit_nm": crb_id0,
           "crb1_ideal_centre_point_nm": crb_id0p,
           "crb1_ideal_disk_mean_nm": float(crb_id.mean())}
    pos = np.stack([XS, np.zeros_like(XS)], axis=1)
    # reference: no cross-talk (tau = 0), ideal counts and ideal MLE -> intrinsic finite-N bias
    rb, rse = [], []
    for x0 in pos:
        cnt = photons.sample_counts(ideal(x0), N_CT, size=R_CT, rng=rng)
        e = estimators.mle(cnt, ideal, SEARCH) - x0
        rb.append(e.mean(axis=0))
        rse.append(e.std(axis=0, ddof=1) / np.sqrt(R_CT))
    rb, rse = np.array(rb), np.array(rse)
    out["ref_tau0"] = {"x_nm": XS.tolist(), "mc_bias_nm": rb.tolist(), "mc_bias_se_nm": rse.tolist(),
                       "max_abs_bias_nm": float(np.abs(rb).max()),
                       "max_abs_z": float(np.abs(rb / rse).max()),
                       "chi2": float(np.sum((rb / rse) ** 2)), "chi2_dof": int(rb.size)}
    print("xtalk ref tau=0: max|b|=%.3f max|z|=%.2f chi2=%.1f/%d" % (
        out["ref_tau0"]["max_abs_bias_nm"], out["ref_tau0"]["max_abs_z"], out["ref_tau0"]["chi2"],
        rb.size))
    for oname, order in ORDERS.items():
        for tau in TAUS:
            f = pm.make_model_crosstalk(c, beam, tau, order=order)
            key = "order%s_tau%g" % (oname, tau)
            d = {"order": list(order), "tau_ns": tau, "M": f.M.tolist(),
                 "q": float(np.exp(-pm.WINDOW_NS / tau))}
            cm = fisher.crb(f, pts, 1.0, zero_policy="limit")
            cc = float(fisher.crb(f, np.zeros(2), 1.0, zero_policy="limit"))
            d["crb_ratio_centre"] = cc / crb_id0             # vs ideal limit r -> 0
            d["crb_ratio_centre_vs_point"] = cc / crb_id0p   # vs ideal point value (Eq. S27)
            # noise-free check: the MLE with M recovers the true position
            popm = estimators.mle(1e6 * f(pos), f, SEARCH, tol=1e-5) - pos
            d["withM_pop_max_abs_bias_nm"] = float(np.abs(popm).max())
            d["crb_ratio_disk_of_means"] = float(cm.mean() / crb_id.mean())
            d["crb_ratio_disk_mean_of_ratios"] = float(np.mean(cm / crb_id))
            # population bias of the naive MLE (noise-free expected counts)
            pop = estimators.mle(1e6 * f(pos), ideal, SEARCH, tol=1e-5) - pos
            d["x_nm"] = XS.tolist()
            d["naive_pop_bias_x_nm"] = pop[:, 0].tolist()
            d["naive_pop_bias_y_nm"] = pop[:, 1].tolist()
            nb, nse, mb, mse, crbN, nrm, mrm = [], [], [], [], [], [], []
            for x0 in pos:
                cnt = photons.sample_counts(f(x0), N_CT, size=R_CT, rng=rng)
                en = estimators.mle(cnt, ideal, SEARCH) - x0
                em = estimators.mle(cnt, f, SEARCH) - x0
                nb.append(en.mean(axis=0))
                nrm.append(float(np.sqrt(0.5 * np.mean(np.sum(en ** 2, axis=1)))))
                mrm.append(float(np.sqrt(0.5 * np.mean(np.sum(em ** 2, axis=1)))))
                nse.append(en.std(axis=0, ddof=1) / np.sqrt(R_CT))
                mb.append(em.mean(axis=0))
                mse.append(em.std(axis=0, ddof=1) / np.sqrt(R_CT))
                crbN.append(float(fisher.crb(f, x0, N_CT, zero_policy="limit")))
            nb, nse, mb, mse = map(np.array, (nb, nse, mb, mse))
            d.update(naive_mc_bias_nm=nb.tolist(), naive_mc_bias_se_nm=nse.tolist(),
                     withM_mc_bias_nm=mb.tolist(), withM_mc_bias_se_nm=mse.tolist(),
                     crb_withM_N500_nm=crbN, naive_mc_rmse_axis_nm=nrm,
                     withM_mc_rmse_axis_nm=mrm)
            dz = (mb - rb) / np.sqrt(mse ** 2 + rse ** 2)
            d["withM_minus_ref_max_abs_z"] = float(np.abs(dz).max())
            d["withM_minus_ref_chi2"] = float(np.sum(dz ** 2))
            z = mb / mse
            d["withM_max_abs_z"] = float(np.abs(z).max())
            d["withM_chi2"] = float(np.sum(z ** 2))
            d["withM_chi2_dof"] = int(z.size)
            d["withM_max_abs_bias_nm"] = float(np.abs(mb).max())
            d["naive_pop_max_abs_bias_nm"] = float(np.max(np.hypot(pop[:, 0], pop[:, 1])))
            d["naive_pop_bias_centre_nm"] = pop[0].tolist()
            out[key] = d
            print("xtalk %s tau=%g crb ratio centre=%.4f disk=%.4f | naive pop |b|max=%.3f "
                  "centre=(%.3f,%.3f) | withM max|z|=%.2f chi2=%.1f/%d; vs ref chi2=%.1f pop=%.1e" % (
                      oname, tau, d["crb_ratio_centre"], d["crb_ratio_disk_of_means"],
                      d["naive_pop_max_abs_bias_nm"], pop[0, 0], pop[0, 1],
                      d["withM_max_abs_z"], d["withM_chi2"], d["withM_chi2_dof"],
                      d["withM_minus_ref_chi2"], d["withM_pop_max_abs_bias_nm"]))
    return out


def make_figure(fl, ct):
    cm = 1 / 2.54
    fig, axs = plt.subplots(1, 3, figsize=(17 * cm, 5.6 * cm), constrained_layout=True)
    # (a) flickering
    ax = axs[0]
    n_int = int(round(T_TOTAL_US / (pm.PERIOD_NS * 1e-3)))
    xr = list(REPS) + [n_int]
    cols = {"centro": S.CYCLE[0], "x20": S.CYCLE[1]}
    labs = {"centro": "centro", "x20": "x = 20 nm"}
    for pname in POS_FL:
        keys = ["seq_r%d" % r for r in REPS] + ["interleaved"]
        pop = [fl[pname][k]["pop_sigma_fl_nm"] for k in keys]
        fN = [fl[pname][k]["fixedN_sigma_fl_nm"] for k in keys]
        ax.plot(xr[:-1], pop[:-1], "-", color=cols[pname], label="%s (sin ruido)" % labs[pname])
        ax.plot(xr[:-1], fN[:-1], "o", mfc="none", color=cols[pname],
                label="%s (MC, N = 100)" % labs[pname])
        ax.plot([xr[-1]], [max(pop[-1], 1e-3)], "*", ms=7, color=cols[pname])
        ax.axhline(fl[pname]["sigma_crb_N100_nm"], color=cols[pname], ls=":", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 30)
    ax.set_xlabel("repeticiones del patrón r")
    ax.set_ylabel(r"$\sigma_{fl}$ (nm)")
    ax.annotate("p-MINFLUX\n(entrelazado)", (n_int, 2e-3), ha="center", va="bottom", fontsize=6)
    ax.text(1.0, 0.72, r"$\sigma_{CRB}$(N = 100) (punteado)", fontsize=6, va="top", ha="right", transform=ax.transAxes)
    ax.legend(fontsize=5.5, loc="lower left")
    ax.set_title("(a) Flickering: secuencial vs entrelazado", fontsize=7.5, loc="left")
    # (b) CRB ratio
    ax = axs[1]
    t = np.array(TAUS)
    for oname, ls in (("A", "-"), ("B", "--")):
        rc = [ct["order%s_tau%g" % (oname, x)]["crb_ratio_centre"] for x in TAUS]
        rd = [ct["order%s_tau%g" % (oname, x)]["crb_ratio_disk_of_means"] for x in TAUS]
        rp = [ct["order%s_tau%g" % (oname, x)]["crb_ratio_centre_vs_point"] for x in TAUS]
        if oname == "A":
            ax.plot(t, rp, "^:", color=S.CYCLE[1], label="centro vs S27 (valor puntual)")
        ax.plot(t, rd, "s" + ls, color=S.CYCLE[2], mfc="none" if oname == "B" else None,
                label=r"media $r\leq L/2$, orden %s" % oname)
    ax.set_xlabel(r"tiempo de vida $\tau$ (ns), T = 12.5 ns")
    ax.set_ylabel("CRB(con M) / CRB(ideal)")
    ax.legend(fontsize=5.5)
    ax.set_title("(b) Cross-talk: pérdida de precisión", fontsize=7.5, loc="left")
    # (c) naive bias
    ax = axs[2]
    for i, tau in enumerate(TAUS):
        col = S.CYCLE[i]
        for oname, ls in (("A", "-"), ("B", "--")):
            d = ct["order%s_tau%g" % (oname, tau)]
            b = np.hypot(d["naive_pop_bias_x_nm"], d["naive_pop_bias_y_nm"])
            ax.plot(XS, b, ls, color=col, lw=1.0,
                    label=r"$\tau$ = %g ns" % tau if oname == "A" else None)
            if oname == "A":
                mb = np.array(d["naive_mc_bias_nm"])
                ax.plot(XS, np.hypot(mb[:, 0], mb[:, 1]), "o", ms=2.5, mfc="none", color=col)
    d = ct["orderA_tau5"]
    mb = np.array(d["withM_mc_bias_nm"])
    ax.plot(XS, np.hypot(mb[:, 0], mb[:, 1]), "x", ms=3, color="k",
            label=r"MLE con M ($\tau$=5, A)")
    ax.set_xlabel("posición x (nm), y = 0")
    ax.text(0.40, 0.97, "línea: sin ruido (orden A)\n--: orden B;  o: MC", transform=ax.transAxes, fontsize=5.5, va="top")
    ax.set_ylabel("|sesgo| del MLE (nm)")
    ax.legend(fontsize=5.5, ncol=1, loc="upper left")
    ax.set_title("(c) Sesgo del MLE que ignora M (N = 500)", fontsize=7.5, loc="left")
    os.makedirs(FIGDIR, exist_ok=True)
    png = os.path.join(FIGDIR, "fig_pminflux_timing.png")
    pdf = os.path.join(FIGDIR, "fig_pminflux_timing.pdf")
    fig.savefig(png, dpi=200)
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    print("wrote", png)
    print("wrote", pdf)


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    c = patterns.tcp_centers(L)
    beam = beams.make_beam("donut", fwhm=FWHM)
    ideal = photons.make_model(c, beam)
    crb1 = {k: float(fisher.crb(ideal, np.array(v), 1.0, zero_policy="limit"))
            for k, v in POS_FL.items()}
    fl = flicker(rng, ideal, crb1)
    ct = crosstalk(rng, c, beam, ideal)
    data = {
        "_script": "scripts/fig_9_pminflux_timing.py",
        "_notes": ("CRB con zero_policy='limit'; STD y sigma por eje = sqrt((var_x+var_y)/2); "
                   "sigma_fl^2 = STD^2 - STD_ctrl^2 (ctrl: mismo MLE, mismos N, sin flickering); "
                   "pop = MLE sobre conteos esperados (sin ruido de Poisson). Orden de pulsos: "
                   "order[s] = exposición disparada en la ranura s (índice 3 = centro)."),
        "params": {"seed": SEED, "L_nm": L, "fwhm_nm": FWHM, "background": "none",
                   "t_total_us": T_TOTAL_US, "t_on_us": T_ON_US, "t_off_us": T_OFF_US,
                   "telegraph_start": "stationary", "n_mean": N_MEAN, "reps": list(REPS),
                   "n_loc": N_LOC, "positions_flicker_nm": {k: list(v) for k, v in POS_FL.items()},
                   "period_ns": pm.PERIOD_NS, "window_ns": pm.WINDOW_NS, "taus_ns": list(TAUS),
                   "orders": {k: list(v) for k, v in ORDERS.items()}, "N_crosstalk": N_CT,
                   "R_crosstalk": R_CT, "disk_step_nm": DISK_STEP, "mle_search_radius_nm": SEARCH,
                   "estimator": "estimators.mle (multinomial, grid + pattern search)"},
        "flicker": fl,
        "crosstalk": ct,
    }
    os.makedirs(DATADIR, exist_ok=True)
    path = os.path.join(DATADIR, "pminflux_timing.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, sort_keys=True, ensure_ascii=False)
    print("wrote", path)
    make_figure(fl, ct)
    print("elapsed %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
