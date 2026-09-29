# -*- coding: utf-8 -*-
"""Explorador guiado de las simulaciones de donut-beam-localization.

Uso:
    python explore/donut_explorer.py                 # menú interactivo
    python explore/donut_explorer.py --mode 3        # un modo, con preguntas
    python explore/donut_explorer.py --mode 3 --yes  # un modo, con los valores por defecto
    python explore/donut_explorer.py --all --yes --no-show --save out_explorer
                                                     # todos los modos, sin ventanas, guarda PNG

Es standalone: agrega ``src/`` al ``sys.path`` y solo necesita numpy, scipy y matplotlib.
La guía de cada modo está en ``explore/guia_explorador.html``. Unidades: nm.
"""

from __future__ import print_function

import argparse
import os
import sys
import textwrap
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

import numpy as np  # noqa: E402

from donutloc import (background, beams, estimators, experiments, fisher,  # noqa: E402
                      patterns, photons, pminflux)

SEED = 42
OPTS = {"yes": False, "show": True, "save": None}


# --------------------------------------------------------------------------- utilidades
def say(text):
    print(textwrap.fill(" ".join(ln.strip() for ln in text.strip().splitlines()), width=92))


def header(title):
    print("\n" + "=" * 92 + "\n" + title + "\n" + "=" * 92)


def ask(prompt, default, cast=float):
    """Pregunta un valor; Enter (o --yes) usa el valor por defecto."""
    if OPTS["yes"] or not sys.stdin.isatty():
        print("  %s [%s]" % (prompt, default))
        return cast(default)
    while True:
        raw = input("  %s [%s]: " % (prompt, default)).strip()
        if raw == "":
            return cast(default)
        try:
            return cast(raw)
        except ValueError:
            print("    valor no válido, probá de nuevo")


def ask_list(prompt, default):
    return [float(v) for v in str(ask(prompt, default, str)).replace(",", " ").split()]


def lesson(lines):
    print("\nQué aprender de esto:")
    for ln in lines:
        print(textwrap.fill(ln, width=92, initial_indent="  - ", subsequent_indent="    "))


def finish(fig, name):
    import matplotlib.pyplot as plt
    fig.tight_layout()
    if OPTS["save"]:
        os.makedirs(OPTS["save"], exist_ok=True)
        path = os.path.join(OPTS["save"], name + ".png")
        fig.savefig(path, dpi=130)
        print("\nFigura guardada en", path)
    if OPTS["show"]:
        plt.show()
    plt.close(fig)


def plt_module():
    import matplotlib
    if not OPTS["show"]:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def tcp_model(L, fwhm, sbr=None, eps=0.0):
    beam = beams.make_beam("donut", fwhm, eps=eps)
    centers = patterns.tcp_centers(L)
    return centers, beam, photons.make_model(centers, beam, sbr=sbr)


# --------------------------------------------------------------------------- modos
def mode_beam():
    header("1. La dona y el patrón TCP")
    say("""Una dona Laguerre-Gauss tiene intensidad cero en el centro y un anillo de máximo. En MINFLUX
    se ubica la dona en K posiciones (el TCP: 3 en un círculo de diámetro L y 1 en el centro) y se
    cuenta cuántos fotones emite la molécula en cada exposición.""")
    fwhm = ask("fwhm de la dona (nm)", 300)
    L = ask("diámetro L del TCP (nm)", 100)
    eps = ask("cero residual eps = I(0)/Imax (0 = cero perfecto)", 0.0)
    plt = plt_module()
    beam = beams.make_beam("donut", fwhm, eps=eps)
    centers = patterns.tcp_centers(L)
    rr = np.linspace(0, 2.5 * fwhm, 400)
    x = np.linspace(-1.5 * L, 1.5 * L, 201)
    X, Y = np.meshgrid(x, x)
    I = photons.intensities(np.stack([X, Y], -1), centers, beam)
    p = I / I.sum(-1, keepdims=True)
    print("\n  radio del anillo: %.1f nm  (0.6006·fwhm)" % beams.ring_radius(fwhm))
    for name, r in (("centro", (0, 0)), ("x = L/4", (L / 4, 0))):
        pr = photons.probabilities(np.array(r, float), centers, beam)
        print("  p_i en %-8s: %s" % (name, np.array2string(pr, precision=3)))
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].plot(rr, beam(rr, 0 * rr))
    ax[0].set(xlabel="r (nm)", ylabel="I / pico", title="perfil radial")
    im = ax[1].imshow(p[..., 3], extent=[x[0], x[-1], x[0], x[-1]], origin="lower")
    ax[1].plot(centers[:, 0], centers[:, 1], "wx")
    ax[1].set(title="p de la exposición central", xlabel="x (nm)", ylabel="y (nm)")
    fig.colorbar(im, ax=ax[1])
    im = ax[2].imshow(p[..., 0], extent=[x[0], x[-1], x[0], x[-1]], origin="lower")
    ax[2].plot(centers[:, 0], centers[:, 1], "wx")
    ax[2].set(title="p de la exposición 0 (arriba)", xlabel="x (nm)")
    fig.colorbar(im, ax=ax[2])
    lesson(["En el centro la exposición central da p = 0 (si eps = 0): es la que más informa, porque "
            "crece como r² apenas la molécula se aleja.",
            "La información está en cómo cambian las fracciones p_i con la posición, no en el brillo "
            "absoluto: MINFLUX mide proporciones.",
            "Probá eps = 0.01: el mínimo ya no es cero y la exposición central pierde contraste."])
    finish(fig, "modo1_dona")


def mode_crb():
    header("2. Límite de Cramér-Rao (CRB): mapa y escalamiento")
    say("""El CRB es la mejor precisión posible (desviación estándar por eje) para N fotones. Se calcula
    con la información de Fisher del modelo multinomial.""")
    L = ask("L (nm)", 50)
    N = ask("fotones N", 100, int)
    sbr = ask("SBR (0 = sin fondo)", 0)
    fwhm = ask("fwhm (nm)", 300)
    plt = plt_module()
    _, _, p = tcp_model(L, fwhm, sbr=(sbr if sbr > 0 else None))
    x = np.linspace(-L, L, 61)
    X, Y = np.meshgrid(x, x)
    t0 = time.time()
    m = fisher.crb_map(p, x, x, N)
    c_lim = fisher.crb(p, np.zeros(2), N, zero_policy="limit")
    c_pt = fisher.crb(p, np.zeros(2), N, zero_policy="point")
    print("\n  CRB en el centro: límite r->0 = %.3f nm, valor puntual (S27) = %.3f nm" % (c_lim, c_pt))
    Ls = np.geomspace(10, 200, 12)
    cL = [fisher.crb(tcp_model(l, fwhm, sbr=(sbr if sbr > 0 else None))[2], np.zeros(2), N,
                     zero_policy="limit") for l in Ls]
    print("  pendiente log(CRB)/log(L) para L <= 50: %.2f   (%.1f s)" %
          (np.polyfit(np.log(Ls[Ls <= 50]), np.log(np.array(cL)[Ls <= 50]), 1)[0], time.time() - t0))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    im = ax[0].imshow(np.log10(m), extent=[x[0], x[-1], x[0], x[-1]], origin="lower")
    fig.colorbar(im, ax=ax[0], label="log10 CRB (nm)")
    ax[0].set(title="CRB, L=%g, N=%d" % (L, N), xlabel="x (nm)", ylabel="y (nm)")
    ax[1].loglog(Ls, cL, "o-")
    ax[1].set(xlabel="L (nm)", ylabel="CRB en el centro (nm)", title="escalamiento con L")
    lesson(["La precisión es mejor dentro del TCP y empeora rápido fuera: el campo útil es ~L.",
            "Sin fondo el CRB escala ~ L/sqrt(N): achicar L mejora linealmente. Con fondo (SBR > 0) ese "
            "beneficio se satura a L chico.",
            "Sin fondo, el valor puntual en el centro (S27) es mayor que el límite r->0 (límite/S27 ~ 2/sqrt(5) "
            "= 0.89): el CRB es discontinuo allí."])
    finish(fig, "modo2_crb")


def mode_estimator():
    header("3. Estimador de máxima verosimilitud (MLE) frente al CRB")
    say("""Simulamos muchas localizaciones de una molécula fija y las estimamos con el MLE. Si el
    estimador es bueno, su dispersión coincide con el CRB y su sesgo es ~0.""")
    L = ask("L (nm)", 100)
    N = ask("fotones N", 200, int)
    sbr = ask("SBR (0 = sin fondo)", 10)
    x0 = ask("posición x de la molécula (nm)", 20)
    n_rep = ask("repeticiones", 1000, int)
    plt = plt_module()
    _, _, p = tcp_model(L, 300, sbr=(sbr if sbr > 0 else None))
    r0 = np.array([x0, 0.0])
    rng = np.random.default_rng(SEED)
    counts = photons.sample_counts(p(r0), N, size=n_rep, rng=rng)
    est = estimators.mle(counts, p, search_radius=1.5 * L)
    err = est - r0
    sd = np.sqrt((err[:, 0].var() + err[:, 1].var()) / 2)
    crb = fisher.crb(p, r0, N, zero_policy="limit")
    print("\n  sesgo = (%.2f, %.2f) nm   STD por eje = %.2f nm   CRB = %.2f nm   STD/CRB = %.2f"
          % (err[:, 0].mean(), err[:, 1].mean(), sd, crb, sd / crb))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(est[:, 0], est[:, 1], ".", ms=2, alpha=0.4)
    ax[0].plot(*r0, "r+", ms=14)
    ax[0].set(aspect="equal", xlabel="x (nm)", ylabel="y (nm)", title="estimaciones")
    ax[1].hist(err[:, 0], 40, density=True, alpha=0.7, label="error en x")
    g = np.linspace(err[:, 0].min(), err[:, 0].max(), 200)
    ax[1].plot(g, np.exp(-g ** 2 / (2 * crb ** 2)) / np.sqrt(2 * np.pi) / crb, "k", label="Gauss(CRB)")
    ax[1].legend()
    ax[1].set(xlabel="error (nm)")
    lesson(["Con fondo y la molécula dentro del TCP, el MLE alcanza el CRB (STD/CRB ~ 1).",
            "Sin fondo y justo en el centro (x = 0, SBR = 0) el modelo no es regular: el MLE puede dar "
            "STD < CRB, pero sesgado hacia el centro. Probalo.",
            "Fuera del TCP (x > L/2) aparecen sesgo y colas: el estimador pierde información."])
    finish(fig, "modo3_estimador")


def mode_zero():
    header("4. Cero imperfecto: ¿conviene achicar L?")
    say("""Con un cero residual eps, la exposición central nunca llega a 0 fotones. A partir de cierto
    punto achicar L no mejora: hay un L óptimo.""")
    eps_list = ask_list("valores de eps (separados por espacio)", "0.002 0.01 0.05")
    N = ask("fotones N", 100, int)
    plt = plt_module()
    Ls = np.geomspace(3, 300, 40)
    fig, ax = plt.subplots(figsize=(6, 4))
    print()
    for eps in eps_list:
        c = [experiments.crb_center(l, N, eps=eps) for l in Ls]
        o = experiments.optimal_L(eps, N=N)
        ax.loglog(Ls, c, label="eps = %g" % eps)
        ax.plot(o["L_opt"], o["crb_opt"], "k*")
        print("  eps = %-6g L_opt = %6.1f nm  CRB_opt = %.2f nm   0.78·fwhm·sqrt(eps) = %.1f nm"
              % (eps, o["L_opt"], o["crb_opt"], 0.78 * 300 * np.sqrt(eps)))
    ax.set(xlabel="L (nm)", ylabel="CRB en el centro (nm)", title="N = %d" % N)
    ax.legend()
    lesson(["El L óptimo sigue ~0.78·fwhm·sqrt(eps) mientras eps <~ 0.01; para ceros peores cae por "
            "debajo de esa regla.",
            "Para un experimento con L fijo (p-MINFLUX): medir eps primero y elegir L cerca del óptimo.",
            "A la derecha del óptimo el CRB crece ~L; a la izquierda lo domina el pedestal del cero."])
    finish(fig, "modo4_cero")


def mode_misalignment():
    header("5. Patrón desalineado: sesgo si se usa el patrón nominal")
    say("""Las donas reales no están exactamente donde se cree. Movemos cada dona una distancia delta
    en una dirección al azar y estimamos con el patrón nominal (ingenuo) y con el verdadero.""")
    delta = ask("desplazamiento delta de cada dona (nm)", 5)
    N = ask("fotones N", 500, int)
    n_pat = ask("patrones de desalineación al azar", 30, int)
    plt = plt_module()
    L, fwhm, sbr = 100.0, 300.0, 10
    beam = beams.make_beam("donut", fwhm)
    nominal = patterns.tcp_centers(L)
    p_nom = photons.make_model(nominal, beam, sbr=sbr)
    rng = np.random.default_rng(SEED)
    b_naive, b_true = [], []
    for _ in range(n_pat):
        true = patterns.perturb_centers(nominal, delta, rng=rng)
        p_true = photons.make_model(true, beam, sbr=sbr)
        expct = N * p_true(np.zeros(2))[None, :]
        b_naive.append(estimators.mle(expct, p_nom, 1.5 * L)[0])
        b_true.append(estimators.mle(expct, p_true, 1.5 * L)[0])
    bn = np.hypot(*np.array(b_naive).T)
    bt = np.hypot(*np.array(b_true).T)
    print("\n  |sesgo| medio en el centro: ingenuo = %.2f nm (%.2f·delta), patrón verdadero = %.3f nm"
          % (bn.mean(), bn.mean() / delta, bt.mean()))
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(*np.array(b_naive).T, "o", label="patrón nominal")
    ax.plot(*np.array(b_true).T, "x", label="patrón verdadero")
    ax.set(aspect="equal", xlabel="sesgo x (nm)", ylabel="sesgo y (nm)", title="delta = %g nm" % delta)
    ax.legend()
    lesson(["El sesgo del estimador ingenuo es ~0.8·delta en el centro (0.78 en promedio sobre muchos patrones; "
            "con pocos patrones varía): es sistemático, no se promedia.",
            "Con el patrón verdadero el sesgo desaparece: la solución es calibrar la posición de las donas.",
            "Regla práctica: para sesgo < 1 nm, las donas tienen que conocerse con error <~ 1.3 nm."])
    finish(fig, "modo5_desalineacion")


def mode_background():
    header("6. Fondo: qué pasa si el estimador lo ignora")
    say("""Agregamos un fondo igual en las 4 exposiciones y estimamos con tres modelos: sin fondo,
    con el fondo exacto y con el fondo como parámetro libre. Usamos las cuentas esperadas (sin
    ruido), así el sesgo que aparece es del modelo, no del azar.""")
    sbr0 = ask("SBR en el centro", 5)
    N = ask("fotones totales N", 500, int)
    plt = plt_module()
    L, fwhm = 100.0, 300.0
    beam = beams.make_beam("donut", fwhm)
    centers = patterns.tcp_centers(L)
    b = background.bg_from_sbr(centers, beam, sbr0)
    p_true = photons.make_model(centers, beam, bg_per_exposure=b)
    p_nobg = photons.make_model(centers, beam)
    xs = np.array([5.0, 10, 20, 30, 40, 50])
    rs = np.stack([xs, 0 * xs], -1)
    expct = N * p_true(rs)
    e_nobg = estimators.mle(expct, p_nobg, 1.5 * L)
    e_true = estimators.mle(expct, p_true, 1.5 * L)
    e_free = background.mle_free_bg(expct, centers, beam, 1.5 * L, b_max=20 * b)[:, :2]
    print("\n  x (nm) | sesgo x sin fondo | con fondo exacto | fondo libre | CRB b conocido | CRB b libre")
    crb_k, crb_f = [], []
    for i, x in enumerate(xs):
        c = background.crb_free_bg(centers, beam, rs[i], b, N)
        crb_k.append(c["sigma_known_b"])
        crb_f.append(c["sigma"])
        print("  %6.0f | %17.2f | %16.3f | %11.3f | %14.2f | %11.2f" % (
            x, e_nobg[i, 0] - x, e_true[i, 0] - x, e_free[i, 0] - x, crb_k[-1], crb_f[-1]))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(xs, e_nobg[:, 0] - xs, "o-", label="sin fondo")
    ax[0].plot(xs, e_true[:, 0] - xs, "s-", label="fondo exacto")
    ax[0].plot(xs, e_free[:, 0] - xs, "^-", label="fondo libre")
    ax[0].set(xlabel="x (nm)", ylabel="sesgo en x (nm)", title="SBR0 = %g" % sbr0)
    ax[0].legend()
    ax[1].plot(xs, crb_k, "o-", label="b conocido")
    ax[1].plot(xs, crb_f, "s-", label="b libre")
    ax[1].set(xlabel="x (nm)", ylabel="CRB (nm)", title="costo de ajustar el fondo")
    ax[1].legend()
    lesson(["Ignorar el fondo sesga hacia afuera: la señal extra en la exposición central se interpreta "
            "como un emisor más lejos del centro (entre +1.6 y +6.8 nm con SBR = 5, según la posición; "
            "también hay sesgo en y).",
            "Ajustar el fondo como parámetro libre elimina el sesgo y cuesta poco cerca del centro.",
            "En x = 0 el modelo sin fondo tiene 3 soluciones simétricas: por eso no se muestra ese punto."])
    finish(fig, "modo6_fondo")


def mode_flicker():
    header("7. Parpadeo del fluoróforo: MINFLUX secuencial frente a p-MINFLUX entrelazado")
    say("""Si la molécula se apaga y se prende (t_on, t_off) mientras se mide, cada exposición ve un
    tiempo 'on' distinto en el modo secuencial. En p-MINFLUX las exposiciones se alternan cada
    12.5 ns, así que todas ven casi lo mismo.""")
    t_on = ask("t_on (µs)", 100)
    t_off = ask("t_off (µs)", 100)
    reps = ask("repeticiones del patrón en modo secuencial", 1, int)
    n = ask("localizaciones", 1500, int)
    plt = plt_module()
    L = 100.0
    _, _, p = tcp_model(L, 300)
    r0 = np.array([20.0, 0.0])
    p0 = p(r0)
    rng = np.random.default_rng(SEED)
    res = {}
    for scheme in ("none", "sequential", "interleaved"):
        w = pminflux.exposure_weights(scheme, n, t_on=t_on, t_off=t_off, reps=reps, rng=rng)
        cnt = pminflux.simulate_flicker_counts(p0, w, 100.0, rng=rng, fixed_N=100)
        e = estimators.mle(cnt, p, 1.5 * L) - r0
        res[scheme] = e
        print("  %-12s STD por eje = %6.2f nm" % (scheme, np.sqrt((e[:, 0].var() + e[:, 1].var()) / 2)))
    s0 = np.sqrt((res["none"][:, 0].var() + res["none"][:, 1].var()) / 2)
    for scheme in ("sequential", "interleaved"):
        s = np.sqrt((res[scheme][:, 0].var() + res[scheme][:, 1].var()) / 2)
        print("  sigma_fl (%s) = sqrt(STD² - STD_sin_parpadeo²) = %.2f nm"
              % (scheme, np.sqrt(max(s ** 2 - s0 ** 2, 0))))
    fig, ax = plt.subplots(figsize=(6, 4))
    for scheme, lab in (("none", "sin parpadeo"), ("sequential", "secuencial, r = %d" % reps),
                        ("interleaved", "entrelazado (p-MINFLUX)")):
        ax.hist(res[scheme][:, 0], 60, histtype="step", density=True, label=lab)
    ax.set(xlabel="error en x (nm)", title="N = 100, x0 = 20 nm")
    ax.legend()
    lesson(["En modo secuencial el parpadeo agrega un error que no baja con más fotones (decenas de nm "
            "con 1 pasada del patrón); repetir el patrón (r) lo reduce.",
            "Con pulsos entrelazados (p-MINFLUX) el error extra es prácticamente cero.",
            "El entrelazado solo protege si el parpadeo es mucho más lento que el período de 50 ns."])
    finish(fig, "modo7_parpadeo")


def mode_crosstalk():
    header("8. Cruce entre ventanas por el tiempo de vida (p-MINFLUX a 20 MHz)")
    say("""Cada pulso excita la molécula, que emite con un retardo exponencial de tiempo de vida tau. Si
    el retardo supera la ventana de 12.5 ns, el fotón se cuenta en la ventana del pulso siguiente.
    Eso mezcla las probabilidades: p_obs = M p.""")
    tau = ask("tiempo de vida tau (ns)", 4)
    N = ask("fotones N (para el CRB)", 500, int)
    plt = plt_module()
    L, fwhm = 100.0, 300.0
    beam = beams.make_beam("donut", fwhm)
    centers = patterns.tcp_centers(L)
    M = pminflux.crosstalk_matrix(tau)
    print("\n  matriz M (columna j = pulso, fila i = ventana donde se cuenta):")
    print(np.array2string(M, precision=4, suppress_small=True, prefix="  "))
    p_ideal = photons.make_model(centers, beam, sbr=20)
    p_ct = pminflux.make_model_crosstalk(centers, beam, tau, sbr=20)
    xs = np.linspace(0, 50, 11)
    rs = np.stack([xs, 0 * xs], -1)
    expct = N * p_ct(rs)
    e_naive = estimators.mle(expct, p_ideal, 1.5 * L)
    e_ok = estimators.mle(expct, p_ct, 1.5 * L)
    bias_naive = np.hypot(*(e_naive - rs).T)
    bias_ok = np.hypot(*(e_ok - rs).T)
    ratio = [fisher.crb(p_ct, r, N) / fisher.crb(p_ideal, r, N) for r in rs]
    print("\n  |sesgo| máximo ignorando M = %.2f nm; con M = %.4f nm; CRB(con M)/CRB(ideal) medio = %.3f"
          % (bias_naive.max(), bias_ok.max(), np.mean(ratio)))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(xs, bias_naive, "o-", label="MLE que ignora M")
    ax[0].plot(xs, bias_ok, "s-", label="MLE con M")
    ax[0].set(xlabel="x (nm)", ylabel="|sesgo| (nm)", title="tau = %g ns, SBR = 20" % tau)
    ax[0].legend()
    ax[1].plot(xs, ratio, "o-")
    ax[1].set(xlabel="x (nm)", ylabel="CRB(con M) / CRB(ideal)", title="pérdida de precisión")
    lesson(["La fracción que cae en la ventana siguiente es ~exp(-12.5/tau): crece rápido con tau.",
            "Ignorar M produce un sesgo sistemático de varios nm para tau >~ 3 ns.",
            "Incluir M en el modelo elimina el sesgo; el costo en precisión es moderado (aquí, sobre x = 0–50 "
            "con SBR = 20: +3 %, +7 % y +14 % para tau = 3, 4 y 5 ns). M se mide del histograma TCSPC."])
    finish(fig, "modo8_cruce")


MODES = [
    ("La dona y el patrón TCP", mode_beam),
    ("Límite de Cramér-Rao: mapa y escalamiento", mode_crb),
    ("MLE frente al CRB (Monte Carlo)", mode_estimator),
    ("Cero imperfecto y L óptimo", mode_zero),
    ("Patrón desalineado", mode_misalignment),
    ("Fondo no modelado", mode_background),
    ("Parpadeo: secuencial vs entrelazado", mode_flicker),
    ("Cruce entre ventanas (tiempo de vida)", mode_crosstalk),
]


def menu():
    while True:
        header("Explorador de simulaciones MINFLUX / p-MINFLUX (donutloc)")
        for i, (name, _) in enumerate(MODES, 1):
            print("  %d. %s" % (i, name))
        print("  0. salir")
        if OPTS["yes"] or not sys.stdin.isatty():
            print("  (sin terminal interactiva: usá --mode N o --all para correr modos)")
            return
        choice = ask("elegí un modo", 0, int)
        if choice == 0:
            return
        if 1 <= choice <= len(MODES):
            MODES[choice - 1][1]()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Explorador guiado de donutloc (ver guia_explorador.html)")
    ap.add_argument("--mode", type=int, help="número de modo (1-%d)" % len(MODES))
    ap.add_argument("--all", action="store_true", help="correr todos los modos en orden")
    ap.add_argument("--yes", action="store_true", help="usar los valores por defecto sin preguntar")
    ap.add_argument("--no-show", action="store_true", help="no abrir ventanas de figuras")
    ap.add_argument("--save", metavar="DIR", help="guardar cada figura como PNG en DIR")
    a = ap.parse_args(argv)
    OPTS.update(yes=a.yes, show=not a.no_show, save=a.save)
    if a.all:
        for _, fn in MODES:
            t0 = time.time()
            fn()
            print("  (%.1f s)" % (time.time() - t0))
    elif a.mode:
        if not 1 <= a.mode <= len(MODES):
            ap.error("--mode debe estar entre 1 y %d" % len(MODES))
        MODES[a.mode - 1][1]()
    else:
        menu()


if __name__ == "__main__":
    main()
