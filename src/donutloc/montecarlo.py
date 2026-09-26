# -*- coding: utf-8 -*-
"""Monte Carlo evaluation of position estimators.

:func:`run_mc` draws ``n_rep`` photon-count vectors at a fixed true position, runs a *batched*
estimator on all of them at once and returns bias / spread / RMSE statistics
(Balzarotti2017 Eq. S41-S42 for bias and standard deviation).
"""

import numpy as np

__all__ = ["run_mc", "sample_counts_from_p"]


def sample_counts_from_p(p, N, n_rep, rng, mode="multinomial"):
    """Draw ``(n_rep, K)`` integer counts from probabilities ``p`` (K,).

    ``mode="multinomial"``: exactly N photons per repetition.  ``mode="poisson"``: independent
    Poisson counts with means ``N * p_i`` (total photons fluctuate, mean N).
    """
    p = np.asarray(p, dtype=float)
    p = p / p.sum()
    if mode == "multinomial":
        return rng.multinomial(int(N), p, size=int(n_rep))
    if mode == "poisson":
        return rng.poisson(float(N) * p, size=(int(n_rep), p.size))
    raise ValueError("mode must be 'multinomial' or 'poisson'")


def run_mc(estimator, p_fn, r_true, N, n_rep, seed=42, mode="multinomial"):
    """Monte Carlo of ``estimator`` at the true position ``r_true``.

    Parameters
    ----------
    estimator : callable ``counts (M, K) -> estimates (M, 2)`` (batched).
    p_fn : model callable ``(..., 2) -> (..., K)``.
    r_true : (2,) true position (nm).
    N : photons per repetition (mean N for ``mode="poisson"``).
    n_rep : number of repetitions.
    seed : seed for ``np.random.default_rng``.
    mode : "multinomial" or "poisson".

    Returns
    -------
    dict with
      ``bias`` (2,)   mean(estimate) - r_true            (Eq. S41),
      ``std`` (2,)    per-axis sample std, ddof=1         (Eq. S42),
      ``sigma``       sqrt((var_x + var_y) / 2)  -- same convention as the CRB (Eq. S13),
      ``rmse``        sqrt(mean(|r_hat - r_true|^2) / 2), per-axis RMS error w.r.t. the true
                      position (includes bias; per-axis normalisation so that for an unbiased
                      estimator rmse ~ sigma, comparable with the CRB),
      ``sigma_err``   approximate standard error of ``sigma``, sigma / sqrt(2 n_valid)
                      (conservative: treats the two axes as one sample of size n_valid),
      ``n_rep``, ``n_valid`` (finite estimates used in the statistics), ``seed``,
      ``estimates`` (n_rep, 2), ``counts`` (n_rep, K).
    """
    rng = np.random.default_rng(seed)
    r_true = np.asarray(r_true, dtype=float)
    p = np.asarray(p_fn(r_true), dtype=float)
    counts = sample_counts_from_p(p, N, n_rep, rng, mode=mode)
    est = np.asarray(estimator(counts), dtype=float).reshape(int(n_rep), 2)
    ok = np.all(np.isfinite(est), axis=1)
    e = est[ok]
    n_valid = int(ok.sum())
    if n_valid < 2:
        raise RuntimeError("fewer than 2 finite estimates")
    bias = e.mean(axis=0) - r_true
    var = e.var(axis=0, ddof=1)
    sigma = float(np.sqrt(0.5 * var.sum()))
    rmse = float(np.sqrt(0.5 * np.mean(np.sum((e - r_true) ** 2, axis=1))))
    return {
        "bias": bias,
        "std": np.sqrt(var),
        "sigma": sigma,
        "rmse": rmse,
        "sigma_err": sigma / np.sqrt(2.0 * n_valid),
        "n_rep": int(n_rep),
        "n_valid": n_valid,
        "seed": seed,
        "estimates": est,
        "counts": counts,
    }
