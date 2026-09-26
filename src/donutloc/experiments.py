# -*- coding: utf-8 -*-
"""Simulation drivers (OBJECTIVE R4 and R5): iterative MINFLUX vs. camera, finite zero depth
(eps x L), and TCP misalignment with a naive estimator.

All drivers are pure functions that return ``dict``s of numpy arrays / floats plus a
``"params"`` entry recording every parameter and the seed (default 42).  They only use the
public API of ``patterns``, ``beams``, ``photons``, ``fisher``, ``estimators`` and
``montecarlo``.  Any beam enters as an argument (default: LG donut of Balzarotti2017 Eq. S17
with fwhm = 300 nm), so that e.g. a vectorial donut can be passed later.

Statistics convention (same as the CRB, Balzarotti2017 Eq. S13): for a set of errors
``e = r_hat - r_true`` (R x 2),

* ``bias  = mean(e)`` (2,),
* ``sigma = sqrt((var_x + var_y)/2)`` (ddof = 1),
* ``rmse  = sqrt(mean |e|^2 / 2)`` (per axis, includes bias),
* standard error of sigma ``sigma_se``: nonparametric bootstrap
  (``montecarlo.bootstrap_sigma_se``); the Gaussian formula ``sigma / (2 sqrt(R))`` is kept as
  ``sigma_se_gauss`` (it underestimates the SE for heavy-tailed errors, e.g. iterative MINFLUX).

Camera reference: the ideal camera ``sigma_PSF / sqrt(N)`` with ``sigma_PSF = 100 nm``
(Balzarotti2017 p. 1 and p. 32), computed inline.
"""

import numpy as np

from . import beams, estimators, fisher, montecarlo, patterns, photons

__all__ = ["iterative_minflux", "iterative_vs_photons", "l_schedule", "eps_L_sweep",
           "optimal_L", "crb_center", "misalignment_study", "error_stats"]

DEFAULT_SEED = 42
DEFAULT_FWHM = 300.0
DEFAULT_SIGMA_PSF = 100.0       # nm, Balzarotti2017 p. 1 / p. 32
DEFAULT_KAPPA = 6.0             # adaptive rule L_{k+1} = max(L_min, kappa sigma_k), see l_schedule
DEFAULT_L_MIN = 25.0            # nm, smallest L of Balzarotti2017 Fig. S1 (p. 53)


# --------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------

def _default_beam(beam, fwhm=DEFAULT_FWHM):
    return beams.make_beam("donut", fwhm=fwhm) if beam is None else beam


def error_stats(est, r_true, n_boot=1000, seed=DEFAULT_SEED):
    """Bias / sigma / rmse of estimates ``est`` (R, 2) w.r.t. ``r_true`` ((2,) or (R, 2)).

    Returns a dict with ``bias`` (2,), ``bias_abs`` (|bias|), ``sigma``, ``rmse``,
    ``sigma_se`` (bootstrap standard error of sigma, :func:`montecarlo.bootstrap_sigma_se` with
    ``n_boot`` resamples and ``seed``; robust to heavy tails), ``sigma_se_gauss``
    (= sigma / (2 sqrt R), valid only for Gaussian errors), ``bias_se`` (2,) (= std / sqrt R)
    and ``n``.  ``n_boot=0`` skips the bootstrap and sets ``sigma_se = sigma_se_gauss``.
    """
    e = np.asarray(est, float) - np.asarray(r_true, float)
    R = e.shape[0]
    if R < 2:
        raise ValueError("need at least 2 estimates")
    bias = e.mean(axis=0)
    var = e.var(axis=0, ddof=1)
    sigma = float(np.sqrt(0.5 * var.sum()))
    se_g = sigma / (2.0 * np.sqrt(R))
    se_b = (montecarlo.bootstrap_sigma_se(e, n_boot=n_boot, seed=seed) if int(n_boot) >= 2
            else se_g)
    return {
        "bias": bias,
        "bias_abs": float(np.hypot(bias[0], bias[1])),
        "bias_se": np.sqrt(var / R),
        "sigma": sigma,
        "sigma_se": se_b,
        "sigma_se_gauss": se_g,
        "rmse": float(np.sqrt(0.5 * np.mean(np.sum(e * e, axis=1)))),
        "n": int(R),
    }


def _split_photons(N_total, n_iter, photon_split):
    N_total = int(N_total)
    if N_total < n_iter:
        raise ValueError("N_total (%d) must be >= number of iterations (%d)" % (N_total, n_iter))
    if isinstance(photon_split, str):
        if photon_split != "equal":
            raise ValueError("photon_split must be 'equal' or a sequence of weights")
        w = np.ones(n_iter)
    else:
        w = np.asarray(photon_split, float)
        if w.shape != (n_iter,) or np.any(w <= 0) or not np.all(np.isfinite(w)):
            raise ValueError("photon_split weights must be %d positive numbers" % n_iter)
    Nk = np.floor(N_total * w / w.sum()).astype(int)
    Nk[-1] += N_total - Nk.sum()           # remainder to the last (finest) iteration
    if np.any(Nk < 1):
        raise ValueError("every iteration needs >= 1 photon, got %r" % (Nk,))
    return Nk


def _center_crb(p_fn, N):
    """CRB at the pattern centre: the r -> 0 limit (fisher.crb_limit, well defined also at a
    perfect zero without background; equal to the point value wherever that is continuous)."""
    return fisher.crb_limit(p_fn, N)


# --------------------------------------------------------------------------------------------
# (a) iterative MINFLUX
# --------------------------------------------------------------------------------------------

def l_schedule(n_iter=4, L0=150.0, L_min=None, rule="fixed", N_k=None, beam=None, sbr=None,
               bg_per_exposure=None, kappa=DEFAULT_KAPPA):
    """Sequence of TCP diameters L_k (nm).

    * ``rule="fixed"``: geometric from ``L0`` to ``L_min``, ``L_k = L0 (L_min/L0)^(k/(n-1))``.
    * ``rule="adaptive"``: ``L_{k+1} = min(L_k, max(L_min, kappa * sigma_k))`` where
      ``sigma_k`` is the centre CRB (r -> 0 limit) of iteration k (TCP of diameter ``L_k``,
      ``N_k[k]`` photons).  The default ``kappa = 6`` sets the next pattern radius to
      ``L_{k+1}/2 = 3 sigma_k``.  Caveat: ``sigma_k`` is the *centre* CRB, which
      **underestimates** the real error of iteration k, because the emitter is not at the
      pattern centre (up to ``L0/4`` off-centre in iteration 0, where the CRB is larger) and the
      MLE is not exactly efficient.  Measured with the defaults (``L0 = 150``, 4 x 250 photons,
      no background, seed 42, 10^4 emitters; ``tests/test_experiments.py``): the centre CRB of
      iteration 0 is 3.47 nm while the real per-axis error is 4.26 nm; ``3 sigma_k`` alone would
      contain only ~94 % of the emitters, and with the ``L_min = 25 nm`` floor (``L_1/2 = 12.5``
      nm) **~97 %** of the emitters (2.65 % outside) lie inside the next pattern radius -- not
      the 99 % a 2D Gaussian of per-axis sigma_k would suggest.  This rule is our own design
      choice: no iterative recipe is given in the sources
      (docs/literature/A_minflux_theory.md section 5).  The schedule is deterministic (same for
      every repetition).

    ``L_min`` defaults to 25 nm (smallest L of Balzarotti2017 Fig. S1).
    """
    n_iter = int(n_iter)
    if n_iter < 1:
        raise ValueError("n_iter must be >= 1")
    L0 = float(L0)
    L_min = DEFAULT_L_MIN if L_min is None else float(L_min)
    if not (L0 > 0 and L_min > 0):
        raise ValueError("L0 and L_min must be > 0")
    if rule == "fixed":
        if n_iter == 1:
            return np.array([L0])
        return L0 * (L_min / L0) ** (np.arange(n_iter) / (n_iter - 1.0))
    if rule == "adaptive":
        if N_k is None:
            raise ValueError("rule='adaptive' needs N_k (photons per iteration)")
        beam = _default_beam(beam)
        Ls = [L0]
        for k in range(n_iter - 1):
            p_fn = photons.make_model(patterns.tcp_centers(Ls[-1]), beam, sbr=sbr,
                                      bg_per_exposure=bg_per_exposure)
            s_k = _center_crb(p_fn, N_k[k])
            Ls.append(min(Ls[-1], max(L_min, float(kappa) * s_k)))
        return np.array(Ls)
    raise ValueError("rule must be 'fixed' or 'adaptive', got %r" % (rule,))


def iterative_minflux(N_total, L_schedule=None, n_iter=4, L0=150.0, L_min=None,
                      photon_split="equal", beam=None, sbr=None, r0_spread=None, n_rep=500,
                      seed=DEFAULT_SEED, recenter=True, rule="fixed", kappa=DEFAULT_KAPPA,
                      bg_per_exposure=None, sigma_psf=DEFAULT_SIGMA_PSF, search_margin=0.25,
                      mle_kw=None):
    """Iterative MINFLUX with re-centring, Monte Carlo over ``n_rep`` emitters.

    Protocol (our own design; the sources only describe the "zoom-in" idea qualitatively,
    Balzarotti2017 pp. 6-7):

    1. The true position of each repetition is uniform in a disk of radius ``r0_spread``
       (default ``L0/4``) around the origin.  The first TCP is centred on the origin.
    2. Iteration k: TCP of diameter ``L_k`` centred on ``c_k`` (``c_0 = 0``); ``N_k`` photons are
       drawn (multinomial, Balzarotti2017 Eq. S3) and the position is estimated by MLE
       (``estimators.mle``) in the disk of radius ``L_k/2 + search_margin*L_k`` around ``c_k``.
    3. ``recenter=True``: ``c_{k+1} = r_hat_k``; ``False``: all TCPs stay on the origin.
    4. The final estimate is that of the last iteration only (no fusion of earlier estimates).

    Background: ``sbr`` = fixed SBR in every iteration (Balzarotti2017 Eq. S30), or
    ``bg_per_exposure`` = fixed background per exposure (Eq. S28), in which case the SBR drops as
    L shrinks (Eq. S32); the centre SBR of each iteration is returned in ``sbr_center``.

    Implementation: the beam is shift invariant (``photons.intensities`` uses
    ``beam(x - cx, y - cy)``), so the model of a TCP centred on ``c`` evaluated at ``r`` equals
    the model of the TCP centred on the origin evaluated at ``r - c``.  All repetitions are
    therefore simulated and estimated as one batch in pattern-relative coordinates.

    Returns a dict with ``sigma``, ``bias``, ``rmse``, ``sigma_se`` (bootstrap, 1000 resamples),
    ``sigma_se_gauss`` (final, per axis, w.r.t. the truth); per-iteration arrays ``L``, ``N_k``,
    ``sigma_iter``, ``sigma_se_iter``, ``rmse_iter``, ``bias_iter``,
    ``crb_center_iter``, ``sbr_center``; ``camera_sigma = sigma_psf/sqrt(N_total)``,
    ``ratio_to_camera = sigma/camera_sigma``, ``crb_all_photons_Lmin`` (centre CRB if all
    N_total photons were spent at the last L), ``r_true``, ``estimates`` and ``params``.
    """
    beam = _default_beam(beam)
    if L_schedule is not None:
        Ls_fixed = np.asarray(L_schedule, float).ravel()
        if Ls_fixed.size < 1 or np.any(Ls_fixed <= 0):
            raise ValueError("L_schedule must be positive")
        n_iter = Ls_fixed.size
    n_iter = int(n_iter)
    Nk = _split_photons(N_total, n_iter, photon_split)
    if L_schedule is not None:
        Ls = Ls_fixed
    else:
        Ls = l_schedule(n_iter, L0, L_min, rule=rule, N_k=Nk, beam=beam, sbr=sbr,
                        bg_per_exposure=bg_per_exposure, kappa=kappa)
    n_rep = int(n_rep)
    if n_rep < 2:
        raise ValueError("n_rep must be >= 2")
    r0_spread = float(Ls[0]) / 4.0 if r0_spread is None else float(r0_spread)
    if r0_spread < 0:
        raise ValueError("r0_spread must be >= 0")
    mle_kw = {} if mle_kw is None else dict(mle_kw)

    rng = np.random.default_rng(seed)
    rad = r0_spread * np.sqrt(rng.uniform(size=n_rep))
    phi = rng.uniform(0.0, 2 * np.pi, size=n_rep)
    r_true = np.stack([rad * np.cos(phi), rad * np.sin(phi)], axis=1)

    c = np.zeros((n_rep, 2))
    est_iter = np.empty((n_iter, n_rep, 2))
    crb_iter = np.empty(n_iter)
    sbr_c = np.empty(n_iter)
    for k in range(n_iter):
        centers = patterns.tcp_centers(Ls[k])
        p_fn = photons.make_model(centers, beam, sbr=sbr, bg_per_exposure=bg_per_exposure)
        rel = r_true - c
        P = np.asarray(p_fn(rel), float)
        P = P / P.sum(axis=1, keepdims=True)
        counts = rng.multinomial(int(Nk[k]), P)
        radius = Ls[k] * (0.5 + float(search_margin))
        est_rel = estimators.mle(counts, p_fn, search_radius=radius, **mle_kw)
        est = c + est_rel
        est_iter[k] = est
        crb_iter[k] = _center_crb(p_fn, Nk[k])
        if bg_per_exposure is not None:
            sbr_c[k] = float(photons.sbr_at(np.zeros(2), centers, beam, bg_per_exposure))
        else:
            sbr_c[k] = np.inf if sbr is None else float(sbr)
        if recenter:
            c = est

    stats = [error_stats(est_iter[k], r_true) for k in range(n_iter)]
    fin = stats[-1]
    cam = float(sigma_psf) / np.sqrt(float(N_total))
    p_last = photons.make_model(patterns.tcp_centers(Ls[-1]), beam, sbr=sbr,
                                bg_per_exposure=bg_per_exposure)
    return {
        "sigma": fin["sigma"],
        "sigma_se": fin["sigma_se"],
        "sigma_se_gauss": fin["sigma_se_gauss"],
        "bias": fin["bias"],
        "bias_abs": fin["bias_abs"],
        "rmse": fin["rmse"],
        "L": np.asarray(Ls, float),
        "N_k": Nk,
        "sigma_iter": np.array([s["sigma"] for s in stats]),
        "sigma_se_iter": np.array([s["sigma_se"] for s in stats]),
        "rmse_iter": np.array([s["rmse"] for s in stats]),
        "bias_iter": np.array([s["bias"] for s in stats]),
        "crb_center_iter": crb_iter,
        "sbr_center": sbr_c,
        "camera_sigma": cam,
        "ratio_to_camera": fin["sigma"] / cam,
        "crb_all_photons_Lmin": _center_crb(p_last, int(N_total)),
        "r_true": r_true,
        "estimates": est_iter[-1],
        "params": {
            "N_total": int(N_total), "n_iter": n_iter, "L0": float(Ls[0]),
            "L_min": float(Ls[-1]), "rule": rule if L_schedule is None else "explicit",
            "kappa": float(kappa), "photon_split": photon_split, "sbr": sbr,
            "bg_per_exposure": bg_per_exposure, "r0_spread": r0_spread, "n_rep": n_rep,
            "seed": seed, "recenter": bool(recenter), "sigma_psf": float(sigma_psf),
            "search_margin": float(search_margin), "beam": _beam_desc(beam),
        },
    }


def iterative_vs_photons(N_list, sigma_psf=DEFAULT_SIGMA_PSF, **kw):
    """Run :func:`iterative_minflux` for each ``N_total`` in ``N_list`` (same seed for each).

    Returns ``N``, ``sigma``, ``sigma_se``, ``rmse``, ``camera_sigma``, ``ratio_to_camera``,
    ``slope`` (least-squares log-log slope of sigma vs N; camera = -0.5), ``slope_rmse``,
    ``runs`` (the individual dicts without the big arrays) and ``params``.
    """
    N_list = np.asarray(N_list)
    runs = []
    for N in N_list:
        out = iterative_minflux(int(N), sigma_psf=sigma_psf, **kw)
        for key in ("r_true", "estimates"):
            out.pop(key)
        runs.append(out)
    sig = np.array([r["sigma"] for r in runs])
    rm = np.array([r["rmse"] for r in runs])
    cam = np.array([r["camera_sigma"] for r in runs])
    lN = np.log(N_list.astype(float))
    slope = float(np.polyfit(lN, np.log(sig), 1)[0]) if len(N_list) > 1 else np.nan
    slope_rm = float(np.polyfit(lN, np.log(rm), 1)[0]) if len(N_list) > 1 else np.nan
    return {
        "N": N_list.astype(int),
        "sigma": sig,
        "sigma_se": np.array([r["sigma_se"] for r in runs]),
        "rmse": rm,
        "camera_sigma": cam,
        "ratio_to_camera": sig / cam,
        "slope": slope,
        "slope_rmse": slope_rm,
        "runs": runs,
        "params": dict(runs[0]["params"], N_list=[int(n) for n in N_list]),
    }


def _beam_desc(beam):
    return {k: getattr(beam, k) for k in ("kind", "fwhm", "eps", "zero_model", "power")
            if hasattr(beam, k)}


# --------------------------------------------------------------------------------------------
# (b) finite zero depth: eps x L
# --------------------------------------------------------------------------------------------

def crb_center(L, N=100, fwhm=DEFAULT_FWHM, eps=0.0, sbr=None, zero_model="gaussian",
               bg_per_exposure=None, return_method=False):
    """Centre CRB (nm) of the TCP of diameter ``L`` with an LG donut with residual zero ``eps``.

    With ``eps > 0`` or any background the CRB is continuous at the centre and the point value
    ``fisher.crb(p_fn, 0, N)`` is used.  With ``eps = 0`` and no background the centre is a
    perfect zero (0/0 term, see fisher module docstring) and the r -> 0 limit
    ``fisher.crb_limit`` is used.
    """
    beam = beams.make_beam("donut", fwhm=fwhm, eps=eps, zero_model=zero_model)
    p_fn = photons.make_model(patterns.tcp_centers(L), beam, sbr=sbr,
                              bg_per_exposure=bg_per_exposure)
    perfect = float(eps) == 0.0 and (sbr is None or np.isinf(sbr)) and not bg_per_exposure
    if perfect:
        val, how = fisher.crb_limit(p_fn, N), "limit"
    else:
        val, how = float(fisher.crb(p_fn, np.zeros(2), N)), "point"
    return (val, how) if return_method else val


def _fov_mean_crb(L, N, fwhm, eps, sbr, zero_model, bg_per_exposure, radius, n_r=8, n_a=12):
    """Area-weighted mean CRB over the disk of ``radius`` (polar midpoint rule, origin
    excluded, so the 0/0 point never appears)."""
    beam = beams.make_beam("donut", fwhm=fwhm, eps=eps, zero_model=zero_model)
    p_fn = photons.make_model(patterns.tcp_centers(L), beam, sbr=sbr,
                              bg_per_exposure=bg_per_exposure)
    rr = radius * (np.arange(n_r) + 0.5) / n_r
    aa = 2 * np.pi * np.arange(n_a) / n_a
    R_, A_ = np.meshgrid(rr, aa, indexing="ij")
    pts = np.stack([R_ * np.cos(A_), R_ * np.sin(A_)], axis=-1)
    vals = fisher.crb(p_fn, pts, N)
    w = np.broadcast_to(rr[:, None], vals.shape)
    return float(np.sum(vals * w) / np.sum(w))


def eps_L_sweep(eps_list, L_list, N=100, fwhm=DEFAULT_FWHM, sbr=None, zero_model="gaussian",
                fov_radius=None, bg_per_exposure=None, fov=False):
    """Centre CRB for every (eps, L) pair.

    Returns ``crb_center`` (n_eps, n_L), ``method`` (n_eps, n_L) of "point"/"limit" (see
    :func:`crb_center`), ``L_best`` (n_eps,) = grid argmin over ``L_list``, and if ``fov=True``
    ``crb_fov_mean`` (n_eps, n_L): area-weighted mean CRB over a disk of radius ``fov_radius``
    (default ``L/4``).  ``params`` records the inputs.
    """
    eps_list = np.atleast_1d(np.asarray(eps_list, float))
    L_list = np.atleast_1d(np.asarray(L_list, float))
    out = np.empty((eps_list.size, L_list.size))
    meth = np.empty(out.shape, dtype=object)
    fovm = np.empty(out.shape) if fov else None
    for i, e in enumerate(eps_list):
        for j, L in enumerate(L_list):
            out[i, j], meth[i, j] = crb_center(L, N, fwhm, e, sbr, zero_model,
                                               bg_per_exposure, return_method=True)
            if fov:
                rad = L / 4.0 if fov_radius is None else float(fov_radius)
                fovm[i, j] = _fov_mean_crb(L, N, fwhm, e, sbr, zero_model, bg_per_exposure,
                                           rad)
    res = {
        "eps": eps_list, "L": L_list, "crb_center": out, "method": meth,
        "L_best": L_list[np.argmin(out, axis=1)],
        "params": {"N": N, "fwhm": float(fwhm), "sbr": sbr, "zero_model": zero_model,
                   "bg_per_exposure": bg_per_exposure, "fov_radius": fov_radius},
    }
    if fov:
        res["crb_fov_mean"] = fovm
    return res


def optimal_L(eps, N=100, fwhm=DEFAULT_FWHM, sbr=None, zero_model="gaussian",
              bg_per_exposure=None, L_bounds=(1.0, None), n_grid=60):
    """L minimizing the centre CRB (:func:`crb_center`) for residual zero ``eps``.

    Coarse log-spaced scan of ``n_grid`` points in ``L_bounds`` (default ``(1 nm, 1.5 fwhm)``)
    followed by a bounded scalar minimization in log L between the neighbours of the best grid
    point (``scipy.optimize.minimize_scalar``).  Returns a dict with ``L_opt``, ``crb_opt``,
    ``interior`` (False if the minimum sits at a bound, e.g. eps = 0 without a background that
    grows as L shrinks, where the CRB decreases monotonically as L -> 0) and ``params``.
    """
    from scipy.optimize import minimize_scalar
    lo = float(L_bounds[0])
    hi = 1.5 * float(fwhm) if L_bounds[1] is None else float(L_bounds[1])
    if not 0 < lo < hi:
        raise ValueError("need 0 < L_min < L_max")
    grid = np.geomspace(lo, hi, int(n_grid))

    def f(L):
        return crb_center(L, N, fwhm, eps, sbr, zero_model, bg_per_exposure)

    vals = np.array([f(L) for L in grid])
    j = int(np.argmin(vals))
    interior = 0 < j < grid.size - 1
    if interior:
        res = minimize_scalar(lambda t: f(np.exp(t)), bounds=(np.log(grid[j - 1]),
                              np.log(grid[j + 1])), method="bounded",
                              options={"xatol": 1e-6})
        L_opt, c_opt = float(np.exp(res.x)), float(res.fun)
    else:
        L_opt, c_opt = float(grid[j]), float(vals[j])
    return {"L_opt": L_opt, "crb_opt": c_opt, "interior": bool(interior),
            "params": {"eps": float(eps), "N": N, "fwhm": float(fwhm), "sbr": sbr,
                       "zero_model": zero_model, "bg_per_exposure": bg_per_exposure,
                       "L_bounds": (lo, hi)}}


# --------------------------------------------------------------------------------------------
# (c) misalignment with naive vs honest estimator
# --------------------------------------------------------------------------------------------

def misalignment_study(displacement_list, L=100.0, N=500, fwhm=DEFAULT_FWHM, sbr=10,
                       n_patterns=20, n_rep=200, seed=DEFAULT_SEED, estimator="mle",
                       positions=None, beam=None, search_margin=0.25):
    """Bias and precision loss of a naive estimator under TCP misalignment.

    For each displacement delta (nm) in ``displacement_list``:

    * ``n_patterns`` perturbed TCPs are drawn with ``patterns.perturb_centers`` (every one of
      the K = 4 zeros moved by delta in an independent uniform random direction; Masullo-type
      model).  The directions come from ``default_rng(seed)`` re-created for every delta, so
      the same directions are used for all deltas (common random numbers; delta = 0 gives the
      ideal TCP exactly).
    * At each true position in ``positions`` (default ``[(0, 0), (L/4, 0)]``) counts are simulated
      with the **true** (perturbed) model by ``montecarlo.run_mc`` (``n_rep`` repetitions).
    * The same counts (identical seed) are estimated with the **honest** model (knows the
      perturbed centres) and the **naive** model (assumes the ideal TCP).

    ``estimator``: ``"mle"`` (``estimators.mle`` in the disk of radius ``L/2 + search_margin*L``
    around the origin) or ``"lms"`` (``estimators.lms`` linearized at the origin).

    Returns a dict with ``delta`` (n_d,), ``positions`` (n_pos, 2) and, for ``"honest"`` and
    ``"naive"``, arrays of shape (n_d, n_pos):

    * ``bias_abs_mean``: mean over patterns of |bias| (nm),
    * ``bias_abs_se``: std over patterns / sqrt(n_patterns),
    * ``bias_chi2``: mean over patterns of |bias|^2 / (var_x/R + var_y/R); its expectation is
      1 for an unbiased estimator (sd ~ 1/sqrt(n_patterns)), >> 1 for a biased one,
    * ``sigma_mean``: mean over patterns of sigma,
    * ``sigma_se``: standard error of ``sigma_mean`` = std over patterns of the per-pattern
      sigmas / sqrt(n_patterns) (includes the pattern-to-pattern variance, which dominates;
      NaN if n_patterns = 1),
    * ``sigma_se_within``: the within-pattern Gaussian formula sigma_mean / (2 sqrt(R P)); it
      ignores the between-pattern variance and underestimates the SE (kept for reference),
    * ``rmse``: sqrt(mean over patterns of rmse^2);

    plus ``crb_honest`` (n_d, n_pos): mean over patterns of the honest CRB (r -> r_true limit,
    ``fisher.crb_limit``), ``bias_noise_floor`` (n_d, n_pos): expected mean |bias| of an unbiased
    estimator from MC noise alone, sigma sqrt(pi/(2R)) (Rayleigh mean, per-axis sd sigma/sqrt R),
    ``identical`` (n_d,): True where honest and naive estimates coincide bit for bit, and
    ``params``.
    """
    beam = _default_beam(beam, fwhm)
    deltas = np.atleast_1d(np.asarray(displacement_list, float))
    if np.any(deltas < 0):
        raise ValueError("displacements must be >= 0")
    if positions is None:
        positions = [(0.0, 0.0), (L / 4.0, 0.0)]
    pos = np.atleast_2d(np.asarray(positions, float))
    P, R = int(n_patterns), int(n_rep)
    if P < 1 or R < 2:
        raise ValueError("n_patterns >= 1 and n_rep >= 2 required")
    ideal = patterns.tcp_centers(L)
    p_naive = photons.make_model(ideal, beam, sbr=sbr)
    radius = L * (0.5 + float(search_margin))

    def make_est(p_fn):
        if estimator == "mle":
            return lambda cnt: estimators.mle(cnt, p_fn, search_radius=radius)
        if estimator == "lms":
            return lambda cnt: estimators.lms(cnt, p_fn)
        raise ValueError("estimator must be 'mle' or 'lms', got %r" % (estimator,))

    est_naive = make_est(p_naive)
    shape = (deltas.size, pos.shape[0])
    keys = ("bias_abs_mean", "bias_abs_se", "bias_chi2", "sigma_mean", "sigma_se",
            "sigma_se_within", "rmse")
    res = {m: {k: np.empty(shape) for k in keys} for m in ("honest", "naive")}
    crb_h = np.empty(shape)
    floor = np.empty(shape)
    identical = np.ones(deltas.size, dtype=bool)
    ss = np.random.SeedSequence(seed)
    mc_seeds = ss.generate_state(P * pos.shape[0]).reshape(P, pos.shape[0])

    for i, d in enumerate(deltas):
        rng = np.random.default_rng(seed)
        acc = {m: {"b": [], "chi": [], "s": [], "r2": []} for m in res}
        acc_crb = np.zeros(pos.shape[0])
        for j in range(P):
            true_c = patterns.perturb_centers(ideal, d, rng=rng)
            p_true = photons.make_model(true_c, beam, sbr=sbr)
            est_honest = make_est(p_true)
            for q in range(pos.shape[0]):
                s_ = int(mc_seeds[j, q])
                outs = {"honest": montecarlo.run_mc(est_honest, p_true, pos[q], N, R, seed=s_),
                        "naive": montecarlo.run_mc(est_naive, p_true, pos[q], N, R, seed=s_)}
                if not np.array_equal(outs["honest"]["estimates"], outs["naive"]["estimates"]):
                    identical[i] = False
                for m, o in outs.items():
                    st = error_stats(o["estimates"], pos[q], n_boot=0)
                    se2 = np.sum(st["bias_se"] ** 2)
                    acc[m]["b"].append(st["bias_abs"])
                    acc[m]["chi"].append(st["bias_abs"] ** 2 / se2 if se2 > 0 else np.nan)
                    acc[m]["s"].append(st["sigma"])
                    acc[m]["r2"].append(st["rmse"] ** 2)
                acc_crb[q] += fisher.crb_limit(p_true, N, r_center=pos[q])
        for m in res:
            for name in ("b", "chi", "s", "r2"):
                acc[m][name] = np.asarray(acc[m][name]).reshape(P, pos.shape[0])
            res[m]["bias_abs_mean"][i] = acc[m]["b"].mean(axis=0)
            res[m]["bias_abs_se"][i] = (acc[m]["b"].std(axis=0, ddof=1) / np.sqrt(P)
                                        if P > 1 else np.nan)
            res[m]["bias_chi2"][i] = acc[m]["chi"].mean(axis=0)
            res[m]["sigma_mean"][i] = acc[m]["s"].mean(axis=0)
            res[m]["sigma_se"][i] = (acc[m]["s"].std(axis=0, ddof=1) / np.sqrt(P)
                                     if P > 1 else np.nan)
            res[m]["sigma_se_within"][i] = acc[m]["s"].mean(axis=0) / (2.0 * np.sqrt(R * P))
            res[m]["rmse"][i] = np.sqrt(acc[m]["r2"].mean(axis=0))
        crb_h[i] = acc_crb / P
        floor[i] = res["honest"]["sigma_mean"][i] * np.sqrt(np.pi / (2.0 * R))

    out = {"delta": deltas, "positions": pos, "honest": res["honest"], "naive": res["naive"],
           "crb_honest": crb_h, "bias_noise_floor": floor, "identical": identical,
           "params": {"L": float(L), "N": int(N), "fwhm": float(fwhm), "sbr": sbr,
                      "n_patterns": P, "n_rep": R, "seed": seed, "estimator": estimator,
                      "search_margin": float(search_margin), "beam": _beam_desc(beam)}}
    return out
