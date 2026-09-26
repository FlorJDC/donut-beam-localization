# -*- coding: utf-8 -*-
"""Monte Carlo evaluation of position estimators.

:func:`run_mc` draws ``n_rep`` photon-count vectors at a fixed true position, runs a *batched*
estimator on all of them at once and returns bias / spread / RMSE statistics
(Balzarotti2017 Eq. S41-S42 for bias and standard deviation).
"""

import numpy as np

__all__ = ["run_mc", "sample_counts_from_p", "bootstrap_sigma_se", "sigma_of_errors"]


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


def sigma_of_errors(errors):
    """Per-axis sigma of an error sample ``errors`` (R, 2): sqrt((var_x + var_y)/2), ddof=1."""
    e = np.asarray(errors, dtype=float)
    return float(np.sqrt(0.5 * e.var(axis=0, ddof=1).sum()))


def bootstrap_sigma_se(errors, n_boot=1000, seed=42, return_samples=False, max_block=2 ** 22):
    """Nonparametric bootstrap standard error of ``sigma = sqrt((var_x + var_y)/2)`` (ddof=1).

    ``errors``: (R, 2) array of estimation errors (or estimates: sigma is shift invariant).
    The R rows are resampled with replacement ``n_boot`` times (``np.random.default_rng(seed)``)
    and the sample std (ddof=1) of the ``n_boot`` bootstrap sigmas is returned.  Unlike the
    Gaussian formula ``sigma / (2 sqrt R)`` it is valid for heavy-tailed / heteroscedastic error
    distributions (e.g. the iterative MINFLUX errors, whose kurtosis exceeds the Gaussian one).
    Non-finite rows are dropped first.  With ``return_samples=True`` returns
    ``(se, sigma_boot)``.  Work is done in blocks of at most ``max_block`` resampled rows.
    """
    e = np.asarray(errors, dtype=float)
    if e.ndim != 2 or e.shape[1] != 2:
        raise ValueError("errors must have shape (R, 2), got %r" % (e.shape,))
    e = e[np.all(np.isfinite(e), axis=1)]
    R = e.shape[0]
    n_boot = int(n_boot)
    if R < 2:
        raise ValueError("need at least 2 finite rows")
    if n_boot < 2:
        raise ValueError("n_boot must be >= 2")
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    chunk = max(1, int(max_block) // R)
    for i0 in range(0, n_boot, chunk):
        m = min(chunk, n_boot - i0)
        idx = rng.integers(0, R, size=(m, R))
        s = e[idx]                                     # (m, R, 2)
        out[i0:i0 + m] = np.sqrt(0.5 * s.var(axis=1, ddof=1).sum(axis=1))
    se = float(out.std(ddof=1))
    return (se, out) if return_samples else se


def run_mc(estimator, p_fn, r_true, N, n_rep, seed=42, mode="multinomial", n_boot=1000):
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
                      position (includes bias).  Normalisation: divided by d = 2, the same
                      per-axis convention as the project CRB sqrt(tr Sigma / 2)
                      (Balzarotti2017 Eq. S13), so that for an unbiased estimator rmse ~ sigma
                      and rmse / CRB is an efficiency; rmse^2 = sigma^2 (n-1)/n + |bias|^2/2.
                      (Masullo's Eq. 4.2 is not in docs/literature and is not cited.)
      ``sigma_err``   approximate standard error of ``sigma``: sigma / (2 sqrt(n_valid)).  For a
                      Gaussian sample the std of a sample std is sigma/sqrt(2n); pooling the two
                      axes (2 n_valid values) gives sigma/sqrt(4 n_valid).  (R1 used the
                      conservative sigma/sqrt(2 n_valid), sqrt(2) too large.)
      ``sigma_se_boot`` bootstrap standard error of ``sigma`` (:func:`bootstrap_sigma_se`,
                      ``n_boot`` resamples, default 1000, seed = ``seed``); robust to
                      non-Gaussian errors.  ``n_boot=0`` skips it (NaN).
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
        "sigma_err": sigma / (2.0 * np.sqrt(n_valid)),
        "sigma_se_boot": (bootstrap_sigma_se(e, n_boot=n_boot, seed=seed) if int(n_boot) >= 2
                          else float("nan")),
        "n_rep": int(n_rep),
        "n_valid": n_valid,
        "seed": seed,
        "estimates": est,
        "counts": counts,
    }
