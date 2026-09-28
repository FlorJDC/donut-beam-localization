# -*- coding: utf-8 -*-
"""Figura 10 -- sesgo por fondo en p-MINFLUX con TCP fijo (analogo TCP de SimuFLUX Fig. 2d).

Caso: TCP fijo (3 donas LG en un circulo de diametro L + centro, 4 exposiciones de igual
duracion, no iterativo), L = 100 nm, fwhm = 300 nm, N = 500 fotones detectados en total
(senal + fondo; multinomial condicionada a N).

Convencion de fondo: fondo constante b por exposicion (Balzarotti2017 Eq. S28, modelo
``photons.probabilities(..., bg_per_exposure=b)``), igual en las 4 ventanas (TCSPC).  b se fija
para que la SBR en el centro del patron (Eq. S29, SBR = sum_j lambda_j / (K b)) sea 5 o 20
(``background.bg_from_sbr``); lejos del centro la SBR cambia y se reporta (``sbr_at_x``).

Posiciones verdaderas: (x, 0), x = 0, 5, 10, 20, 30, 40, 50 nm.  En cada punto se simulan
N_REP vectores de cuentas (mismas cuentas para todos los estimadores) y se estiman con MLE:

  sin_fondo   modelo sin fondo (p_i proporcional a lambda_i);
  fondo_exacto modelo con el b verdadero;
  fondo_libre b como parametro libre (MLE en x, y, b; ``background.mle_free_bg``);
  fondo_m30 / fondo_p30  modelo con b mal estimado en -30 % / +30 %.

Se reporta sesgo (media - verdad) por eje con su error estandar, STD por eje, sigma =
sqrt((var_x+var_y)/2), RMSE por eje (``montecarlo.run_mc``: sqrt(mean|e|^2/2)), el sesgo
"poblacional" (MLE sobre las cuentas esperadas N p, sin ruido) y los CRB: con b conocido
(``fisher.crb`` del modelo correcto) y con b libre (bloque xy de la inversa de la Fisher 3x3,
``background.crb_free_bg``).

Salidas: informe/figures/fig_background_bias.{png,pdf}, informe/data/background_bias.json.
Semilla: np.random.default_rng([42, i_sbr, i_x]) por punto.  Uso: python scripts/fig_10_background_bias.py [--quick]
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _paperconfig as C  # noqa: E402  (pone src/ en sys.path)
import _paperstyle as S  # noqa: E402  (rcParams)

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from donutloc import background, beams, estimators, fisher, patterns, photons  # noqa: E402

SEED = C.SEED
L = 100.0
FWHM = 300.0
N = 500
SBR_CENTER = (5.0, 20.0)
XS = (0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0)
N_REP = 4000
SEARCH_RADIUS = 1.5 * L          # disco de busqueda de todos los MLE, centrado en el TCP
B_MAX = 0.3                      # extremo de la grilla inicial de b (SBR en el centro ~ 0.5)
MIS = 0.3                        # error relativo del fondo mal estimado
EST = ("sin_fondo", "fondo_exacto", "fondo_libre", "fondo_m30", "fondo_p30")
LABEL = {"sin_fondo": "MLE sin fondo", "fondo_exacto": "MLE fondo exacto",
         "fondo_libre": "MLE fondo libre (x, y, b)", "fondo_m30": u"MLE fondo −30 %",
         "fondo_p30": "MLE fondo +30 %"}
COLOR = {"sin_fondo": "#D55E00", "fondo_exacto": "#0072B2", "fondo_libre": "#009E73",
         "fondo_m30": "#CC79A7", "fondo_p30": "#E69F00"}
ROOT = C.ROOT
OUT_FIG = os.path.join(ROOT, "informe", "figures")
OUT_DATA = os.path.join(ROOT, "informe", "data")


def estimator_fns(centers, beam, b):
    def mk(bg):
        return photons.make_model(centers, beam, bg_per_exposure=bg)

    models = {"sin_fondo": photons.make_model(centers, beam), "fondo_exacto": mk(b),
              "fondo_m30": mk((1 - MIS) * b), "fondo_p30": mk((1 + MIS) * b)}
    fns = {k: (lambda cts, m=m: estimators.mle(cts, m, SEARCH_RADIUS)) for k, m in models.items()}
    fns["fondo_libre"] = lambda cts: background.mle_free_bg(cts, centers, beam, SEARCH_RADIUS, B_MAX)
    return fns


def stats(est, r_true):
    e = est[:, :2] - r_true
    M = e.shape[0]
    var = e.var(axis=0, ddof=1)
    sigma = float(np.sqrt(0.5 * var.sum()))
    return {"bias_x": float(e[:, 0].mean()), "bias_y": float(e[:, 1].mean()),
            "bias_x_se": float(np.sqrt(var[0] / M)), "bias_y_se": float(np.sqrt(var[1] / M)),
            "std_x": float(np.sqrt(var[0])), "std_y": float(np.sqrt(var[1])),
            "sigma": sigma, "sigma_se": sigma / (2.0 * np.sqrt(M)),
            "rmse": float(np.sqrt(0.5 * np.mean(np.sum(e ** 2, axis=1))))}


def compute(n_rep):
    centers = patterns.tcp_centers(L)
    beam = beams.make_beam("donut", fwhm=FWHM)
    out = {}
    for i_s, sbr in enumerate(SBR_CENTER):
        b = background.bg_from_sbr(centers, beam, sbr)
        fns = estimator_fns(centers, beam, b)
        model = photons.make_model(centers, beam, bg_per_exposure=b)
        res = {"bg_per_exposure": b, "sbr_at_x": [], "bg_fraction_at_x": [],
               "crb_known_b_nm": [], "crb_free_b_nm": [], "crb_free_b_x_nm": [],
               "crb_free_b_y_nm": [], "crb_free_b_sigma_b": [],
               "fondo_libre_b_hat_mean": [], "fondo_libre_b_hat_std": [],
               "fondo_libre_frac_b_zero": []}
        for k in EST:
            res[k] = {}
        for i_x, x in enumerate(XS):
            r_true = np.array([x, 0.0])
            p = model(r_true)
            sb = float(photons.sbr_at(r_true, centers, beam, b))
            res["sbr_at_x"].append(sb)
            res["bg_fraction_at_x"].append(1.0 / (1.0 + sb))
            res["crb_known_b_nm"].append(float(fisher.crb(model, r_true, N)))
            cf = background.crb_free_bg(centers, beam, r_true, b, N)
            res["crb_free_b_nm"].append(cf["sigma"])
            res["crb_free_b_x_nm"].append(cf["sigma_x"])
            res["crb_free_b_y_nm"].append(cf["sigma_y"])
            res["crb_free_b_sigma_b"].append(cf["sigma_b"])
            rng = np.random.default_rng([SEED, i_s, i_x])
            counts = rng.multinomial(N, p, size=n_rep)
            for k in EST:
                est = np.atleast_2d(fns[k](counts))
                st = stats(est, r_true)
                pop = np.atleast_1d(fns[k](N * p))
                st["bias_pop_x"] = float(pop[0] - x)
                st["bias_pop_y"] = float(pop[1])
                for kk, v in st.items():
                    res[k].setdefault(kk, []).append(v)
                if k == "fondo_libre":
                    res["fondo_libre_b_hat_mean"].append(float(est[:, 2].mean()))
                    res["fondo_libre_b_hat_std"].append(float(est[:, 2].std(ddof=1)))
                    res["fondo_libre_frac_b_zero"].append(float(np.mean(est[:, 2] < 1e-6 * b)))
            print("SBR0=%g x=%g listo" % (sbr, x), flush=True)
        out["sbr%g" % sbr] = res
    return out


def plot(data, path_base):
    fig, axes = plt.subplots(1, 2, figsize=(17 / 2.54, 6.2 / 2.54), constrained_layout=True)
    xs = np.asarray(XS)
    mk = {5.0: ("o", "-"), 20.0: ("s", "--")}
    for sbr in SBR_CENTER:
        d = data["sbr%g" % sbr]
        m, ls = mk[sbr]
        for k in EST:
            r = d[k]
            fill = COLOR[k] if sbr == 5.0 else "white"
            axes[0].errorbar(xs, r["bias_x"], yerr=r["bias_x_se"], color=COLOR[k], ls=ls, marker=m,
                             mfc=fill, lw=0.9, ms=3.5, capsize=0)
            axes[1].plot(xs, r["rmse"], color=COLOR[k], ls=ls, marker=m, mfc=fill, lw=0.9, ms=3.5)
            if k == "sin_fondo":      # unico caso con RMSE >> STD: se muestra tambien la STD
                axes[1].plot(xs, r["sigma"], color=COLOR[k], ls=":", lw=0.9)
        axes[1].plot(xs, d["crb_known_b_nm"], color="k", ls=ls, lw=1.4)
        axes[1].plot(xs, d["crb_free_b_nm"], color="0.55", ls=ls, lw=1.4)
    axes[0].axhline(0.0, color="0.6", lw=0.5)
    axes[0].set_xlabel(u"posici\u00f3n verdadera $x$ (nm)")
    axes[0].set_ylabel(r"sesgo en $x$ (nm)")
    axes[1].set_xlabel(u"posición verdadera $x$ (nm)")
    axes[1].set_ylabel(u"RMSE por eje (nm)")
    axes[1].set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    axes[1].yaxis.set_major_locator(FixedLocator([1.5, 2, 3, 4, 6, 8]))
    axes[1].yaxis.set_major_formatter(ScalarFormatter())
    axes[1].yaxis.set_minor_formatter(NullFormatter())
    axes[0].set_ylim(-1.5, 8.5)
    S.panel_label(axes[0], "a")
    S.panel_label(axes[1], "b")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=COLOR[k], lw=1.2, label=LABEL[k]) for k in EST]
    axes[0].legend(handles=h, loc="upper center", fontsize=6, ncol=2)
    h2 = [Line2D([], [], color="k", marker="o", ls="-", label=u"SBR$_0$ = 5"),
          Line2D([], [], color="k", marker="s", mfc="white", ls="--", label=u"SBR$_0$ = 20"),
          Line2D([], [], color="k", lw=1.4, label=u"$\\sigma_{CRB}$, b conocido"),
          Line2D([], [], color="0.55", lw=1.4, label=u"$\\sigma_{CRB}$, b libre"),
          Line2D([], [], color=COLOR["sin_fondo"], ls=":", label=u"STD, MLE sin fondo")]
    axes[1].legend(handles=h2, loc="upper left", fontsize=6, ncol=2)
    axes[1].set_ylim(1.4, 13)
    os.makedirs(os.path.dirname(path_base), exist_ok=True)
    fig.savefig(path_base + ".png", dpi=200)
    fig.savefig(path_base + ".pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def main():
    quick = "--quick" in sys.argv
    n_rep = 300 if quick else N_REP
    t0 = time.time()
    data = compute(n_rep)
    meta = {"L_nm": L, "fwhm_nm": FWHM, "N_total": N, "sbr_center": list(SBR_CENTER),
            "x_nm": list(XS), "n_rep": n_rep, "seed": SEED,
            "seed_scheme": "np.random.default_rng([42, i_sbr, i_x]); mismas cuentas para todos los estimadores",
            "search_radius_nm": SEARCH_RADIUS, "b_max_free_fit": B_MAX, "misestimate_rel": MIS,
            "estimators": list(EST), "photon_model": "multinomial, N total fijo, fondo b igual por exposicion (Balzarotti Eq. S28)",
            "sbr_definition": "SBR(r) = sum_j lambda_j(r) / (K b), K = 4 (Balzarotti Eq. S29); SBR_0 en r = 0",
            "rmse_definition": "sqrt(mean(|r_hat - r|^2)/2)", "sigma_definition": "sqrt((var_x+var_y)/2)",
            "quick": quick}
    data["meta"] = meta
    sub = "quick" if quick else ""
    dpath = os.path.join(OUT_DATA, sub, "background_bias.json")
    os.makedirs(os.path.dirname(dpath), exist_ok=True)
    with open(dpath, "w") as f:
        json.dump(data, f, indent=1, sort_keys=True)
    plot(data, os.path.join(OUT_FIG, sub, "fig_background_bias"))
    print("escrito %s (%.1f s)" % (dpath, time.time() - t0))


if __name__ == "__main__":
    main()
