# -*- coding: utf-8 -*-
"""Position estimators for MINFLUX-type measurements.

All estimators take photon counts ``counts`` of shape ``(K,)`` (one measurement) or ``(M, K)``
(a batch of M measurements) and return ``(2,)`` or ``(M, 2)`` positions in nm.

* :func:`mle` -- maximum likelihood (Balzarotti2017 Eq. S33/S37) for any model ``p_fn``;
  vectorized over the batch.
* :func:`lms` -- linearized least-mean-squares, general form (Balzarotti2017 Eq. S48) with the
  full numerical Jacobian of ``p_fn``.
* :func:`lms_tcp` -- closed form for the TCP with an LG donut (Balzarotti2017 Eq. S49-S50).
* :func:`mlms_tcp` -- modified LMS (Balzarotti2017 Eq. S51, k = 1).

TCP convention of this project: the three peripheral zeros sit on a circle of **diameter** L at
angles ``rotation + 2 pi k / 3`` (k = 0, 1, 2; default ``rotation = pi/2``) and the central
exposure is the **last** one (index 3).  Balzarotti numbers the centre 0 and puts the peripheral
zeros at ``2 pi i / 3``; the closed forms below are written with the explicit vectors r_b_i,
so the two conventions give the same estimate for the same physical geometry.
"""

import numpy as np

__all__ = ["neg_loglike", "mle", "lms", "lms_tcp", "mlms_tcp"]

_P_FLOOR = 1e-300


def neg_loglike(r, counts, p_fn):
    """Negative multinomial log-likelihood up to an r-independent constant,
    ``-sum_i n_i ln p_i(r)`` (Balzarotti2017 Eq. S37), with p clipped to 1e-300.

    ``r`` is ``(..., 2)``; ``counts`` is ``(K,)`` or broadcastable against ``p_fn(r)``
    (shape ``(..., K)``).  Returns shape ``(...)``.
    """
    p = np.asarray(p_fn(np.asarray(r, dtype=float)), dtype=float)
    return -np.sum(np.asarray(counts, dtype=float) * np.log(np.maximum(p, _P_FLOOR)), axis=-1)


def _disk_grid(center, radius, step):
    """Square-lattice points inside the closed disk, sorted by distance to the centre."""
    n = int(np.floor(radius / step + 1e-9))
    g = step * np.arange(-n, n + 1)
    X, Y = np.meshgrid(g, g, indexing="xy")
    pts = np.stack([X.ravel(), Y.ravel()], axis=-1)
    d2 = np.sum(pts ** 2, axis=-1)
    inside = d2 <= radius ** 2 * (1 + 1e-12)
    pts, d2 = pts[inside], d2[inside]
    order = np.argsort(d2, kind="stable")
    return pts[order] + np.asarray(center, dtype=float)


_MEM_BUDGET = 64 * 2 ** 20          # bytes for the (chunk, G) log-likelihood block of mle


def mle(counts, p_fn, search_radius, center=(0.0, 0.0), grid_step=None, refine=True,
        tol=1e-4, chunk=4096, mem_budget=_MEM_BUDGET):
    """Maximum-likelihood position estimate (Balzarotti2017 Eq. S33, S37) inside a disk.

    Algorithm (fully vectorized over the batch):

    1. **Global grid**: the log-likelihood ``sum_i n_i ln p_i`` is evaluated on a square
       lattice of step ``grid_step`` (default ``search_radius / 25``) restricted to the disk
       ``|r - center| <= search_radius``, as one matrix product ``counts @ log p(grid)^T``.
    2. **Refinement** (``refine=True``): vectorized coarse-to-fine pattern search, in the spirit of
       the successive grids of Balzarotti2017 (p. 22, p. 42): around the current best point a
       5x5 local grid of half-width 2*step is evaluated, the best point is kept, and the step is
       halved, until ``step < tol`` (nm).  Candidates outside the disk are rejected, so the result
       always lies in the disk.  This is a continuous (not lattice-limited) optimizer with
       resolution ``tol``; it replaces a per-sample ``scipy.optimize.minimize`` for speed.
       ``refine="scipy"`` instead polishes each sample with ``scipy.optimize.minimize``
       (Nelder-Mead on the negative log-likelihood, starting from the grid optimum, result
       clipped to the disk).  ``refine=False`` returns the grid optimum.

    Branch choice with multiple maxima: the estimate is the refinement of the **global maximum
    over the grid**; exact ties (e.g. all counts zero, or symmetric likelihoods) are broken in
    favour of the grid point **closest to** ``center`` (the lattice is sorted by radius).
    Refinement is local, so a secondary maximum narrower than the grid step can be missed.

    Robustness: ``n_i = 0`` terms contribute 0; ``p_i = 0`` is clipped to 1e-300 so the
    log-likelihood stays finite; an all-zero count vector returns ``center``.

    Validation (``ValueError``): ``search_radius``, ``grid_step`` (if given) and ``tol`` must be
    finite and > 0 (``tol <= 0`` would never stop the pattern search); ``chunk`` >= 1;
    ``mem_budget`` > 0; counts must be finite and >= 0; the number of exposures K of ``counts``
    must equal the length of ``p_fn(center)``.

    Memory: the global grid is processed in batch chunks of
    ``min(chunk, max(1, floor(mem_budget / (8 G))))`` rows, G = number of grid points, so the
    ``(rows, G)`` float64 log-likelihood block stays below ``mem_budget`` bytes (default 64 MiB)
    whatever ``grid_step``.

    Parameters
    ----------
    counts : array ``(K,)`` or ``(M, K)``.
    p_fn : model callable ``(..., 2) -> (..., K)``.
    search_radius : float, radius (nm) of the search disk.
    center : (2,), centre of the search disk.
    grid_step : float or None.
    refine : bool or "scipy".
    tol : float, final step of the pattern search (nm).
    chunk : int, maximum batch chunk size for the global grid.
    mem_budget : float, bytes allowed for one ``(rows, G)`` block of the global grid
        (default 64 MiB); the effective chunk is ``min(chunk, max(1, floor(mem_budget/(8 G))))``.

    Returns
    -------
    ``(2,)`` or ``(M, 2)`` estimates in nm.
    """
    counts = np.asarray(counts, dtype=float)
    if counts.ndim not in (1, 2):
        raise ValueError("counts must have shape (K,) or (M, K), got %r" % (counts.shape,))
    single = counts.ndim == 1
    C = np.atleast_2d(counts)
    M = C.shape[0]
    center = np.asarray(center, dtype=float)
    if center.shape != (2,):
        raise ValueError("center must have shape (2,), got %r" % (center.shape,))
    _check_pos("search_radius", search_radius)
    _check_pos("tol", tol)
    if grid_step is not None:
        _check_pos("grid_step", grid_step)
    _check_pos("mem_budget", mem_budget)
    if int(chunk) < 1:
        raise ValueError("chunk must be >= 1, got %r" % (chunk,))
    if not np.all(np.isfinite(C)):
        raise ValueError("counts must be finite")
    if np.any(C < 0):
        raise ValueError("counts must be >= 0")
    K = np.asarray(p_fn(center), float).shape[-1]
    if C.shape[1] != K:
        raise ValueError("counts have K = %d exposures but p_fn returns K = %d probabilities"
                         % (C.shape[1], K))
    step = float(search_radius) / 25.0 if grid_step is None else float(grid_step)

    grid = _disk_grid(center, search_radius, step)                   # (G, 2)
    logp = np.log(np.maximum(np.asarray(p_fn(grid), float), _P_FLOOR))  # (G, K)
    G = grid.shape[0]
    chunk = int(min(int(chunk), max(1, int(np.floor(float(mem_budget) / (8.0 * G))))))
    best = np.empty((M, 2))
    for s in range(0, M, chunk):
        ll = C[s:s + chunk] @ logp.T                                  # (m, G)
        best[s:s + chunk] = grid[np.argmax(ll, axis=1)]

    if refine is True:
        best = _pattern_search(C, p_fn, best, center, float(search_radius), step, tol)
    elif refine == "scipy":
        best = _scipy_refine(C, p_fn, best, center, float(search_radius))
    elif refine is not False:
        raise ValueError("refine must be True, False or 'scipy'")
    return best[0] if single else best


def _check_pos(name, v):
    try:
        ok = np.isfinite(float(v)) and float(v) > 0
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise ValueError("%s must be a finite number > 0, got %r" % (name, v))


def _pattern_search(C, p_fn, start, center, radius, step, tol):
    off = np.arange(-2, 3, dtype=float)
    OX, OY = np.meshgrid(off, off, indexing="xy")
    offsets = np.stack([OX.ravel(), OY.ravel()], axis=-1)            # (25, 2); index 12 = (0,0)
    # put the zero offset first so that ties keep the current point
    order = np.argsort(np.sum(offsets ** 2, axis=-1), kind="stable")
    offsets = offsets[order]
    cur = start.copy()
    s = step / 2.0
    r2max = radius ** 2 * (1 + 1e-12)
    while True:
        cand = cur[:, None, :] + s * offsets[None, :, :]              # (M, 25, 2)
        p = np.asarray(p_fn(cand), float)                             # (M, 25, K)
        ll = np.sum(C[:, None, :] * np.log(np.maximum(p, _P_FLOOR)), axis=-1)
        out = np.sum((cand - center) ** 2, axis=-1) > r2max
        ll[out] = -np.inf
        cur = cand[np.arange(cand.shape[0]), np.argmax(ll, axis=1)]
        if s < tol:
            break
        s /= 2.0
    return cur


def _scipy_refine(C, p_fn, start, center, radius):
    from scipy.optimize import minimize
    out = start.copy()
    for m in range(C.shape[0]):
        c = C[m]

        def f(x, c=c):
            if np.sum((x - center) ** 2) > radius ** 2:
                return np.inf
            return float(neg_loglike(x, c, p_fn))

        res = minimize(f, start[m], method="Nelder-Mead",
                       options={"xatol": 1e-5, "fatol": 1e-10, "maxiter": 2000})
        if np.isfinite(res.fun) and res.fun <= f(start[m]):
            out[m] = res.x
    return out


def _phat(counts):
    counts = np.asarray(counts, dtype=float)
    single = counts.ndim == 1
    C = np.atleast_2d(counts)
    Ntot = C.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        P = C / Ntot
    return P, single, Ntot[:, 0]


def lms(counts, p_fn, r_lin=(0.0, 0.0), h=1e-3):
    """Linearized LMS estimator, general form (Balzarotti2017 Eq. S48):

        r_hat = r_lin + (J^T J)^-1 J^T (p_hat - p(r_lin)),   p_hat = n / N  (Eq. S45),

    with the full K x 2 Jacobian ``J = d p / d r`` at ``r_lin`` from centred differences of
    step ``h``.  If ``p_fn`` contains background (Eq. S30) the Jacobian is scaled by
    s = SBR/(SBR+1) and the 1/s correction (note A section 3.2) appears automatically.
    Measurements with zero total counts return ``r_lin``.
    """
    P, single, Ntot = _phat(counts)
    r_lin = np.asarray(r_lin, dtype=float)
    pts = np.stack([r_lin, r_lin + [h, 0.0], r_lin - [h, 0.0], r_lin + [0.0, h], r_lin - [0.0, h]])
    Q = np.asarray(p_fn(pts), float)
    p0 = Q[0]
    J = np.stack([(Q[1] - Q[2]) / (2 * h), (Q[3] - Q[4]) / (2 * h)], axis=-1)   # (K, 2)
    A = np.linalg.solve(J.T @ J, J.T)                                           # (2, K)
    est = r_lin + (P - p0) @ A.T
    est[Ntot <= 0] = r_lin
    return est[0] if single else est


def _tcp_ring(L, rotation):
    ang = rotation + 2.0 * np.pi * np.arange(3) / 3.0
    return 0.5 * L * np.stack([np.cos(ang), np.sin(ang)], axis=-1)              # (3, 2)


def _tcp_gain(L, fwhm, sbr):
    q = L ** 2 * np.log(2.0) / fwhm ** 2 if np.isfinite(fwhm) else 0.0
    if q >= 1.0:
        raise ValueError("LMS undefined for L^2 ln2 / fwhm^2 >= 1 (L=%g, fwhm=%g)" % (L, fwhm))
    g = 1.0 / (1.0 - q)
    if sbr is not None and np.isfinite(sbr):
        g *= (sbr + 1.0) / sbr                      # 1/s, s = SBR/(SBR+1)  (note A section 3.2)
    return g


def lms_tcp(counts, L, fwhm, sbr=None, rotation=np.pi / 2):
    """Closed-form LMS for the TCP with LG donut (Balzarotti2017 Eq. S49-S50), centre last:

        r_hat = -(1 - L^2 ln2/fwhm^2)^-1 * (1/s) * sum_{i=peripheral} p_hat_i r_b_i ,

    where ``r_b_i = (L/2)(cos a_i, sin a_i)``, ``a_i = rotation + 2 pi i/3``, and
    ``s = SBR/(SBR+1)`` if ``sbr`` is given (``None``/``inf``: no background; ``s = 1``).
    ``fwhm = np.inf`` gives the quadratic limit.  The central count enters only through N.
    Zero total counts return (0, 0).
    """
    P, single, Ntot = _phat(counts)
    if P.shape[1] != 4:
        raise ValueError("lms_tcp expects K = 4 counts (3 peripheral + centre last)")
    rb = _tcp_ring(float(L), rotation)
    est = -_tcp_gain(float(L), float(fwhm), sbr) * (P[:, :3] @ rb)
    est[Ntot <= 0] = 0.0
    return est[0] if single else est


def mlms_tcp(counts, L, fwhm, beta=(1.27, 3.8), sbr=None, rotation=np.pi / 2):
    """Modified LMS (Balzarotti2017 Eq. S51):

        r_hat = (sum_j beta_j p_hat_0^j) * r_hat_LMS ,

    with ``p_hat_0`` the fraction of photons in the **central** exposure (last index here) and
    ``r_hat_LMS`` from :func:`lms_tcp` (including the 1/s factor if ``sbr`` is given; the
    polynomial uses the raw ``p_hat_0``).  ``beta`` may have any length (order k = len-1);
    default k = 1, beta = (1.27, 3.8) as used for live tracking (Balzarotti2017 p. 35, 44).
    The estimator is biased by design (Eq. S53).
    """
    P, single, Ntot = _phat(counts)
    base = np.atleast_2d(lms_tcp(counts, L, fwhm, sbr=sbr, rotation=rotation))
    p0 = np.nan_to_num(P[:, 3])
    poly = np.zeros_like(p0)
    for j, b in enumerate(beta):
        poly = poly + b * p0 ** j
    est = poly[:, None] * base
    return est[0] if single else est
