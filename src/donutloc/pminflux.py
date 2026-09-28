# -*- coding: utf-8 -*-
"""Timing effects specific to pulsed-interleaved MINFLUX (p-MINFLUX, Masullo et al. 2021).

p-MINFLUX: a pulsed laser (repetition period ``P``, e.g. 50 ns at 20 MHz) is split into ``K``
interleaved pulse trains (here the 3 TCP donuts + the centre, K = 4), delayed by ``T = P/K``
(12.5 ns).  Photons are assigned to exposures by time-correlated single-photon counting (TCSPC)
in ``K`` time windows of equal width ``T``.  All exposures last the whole dwell (not
iterative, L fixed).  Times in this module are in **microseconds** for the fluorophore kinetics
and in **nanoseconds** for the lifetime / window arithmetic (each function says which).

Two effects are modelled.

1. **Cross-talk between time windows (fluorescence lifetime).**  Slot ``s`` of the period
   (``s = 0 .. K-1``) fires the exposure ``order[s]``; its pulse arrives at ``s T`` (mod P) and
   its time window is ``[s T, (s+1) T)`` (windows aligned with the pulses, no IRF, no jitter).
   A photon excited by that pulse is emitted after a delay ``t ~ Exp(tau)``, so it lands in
   window ``(s + floor(t/T)) mod K`` (the train is periodic, so delays longer than one period
   wrap around).  With ``q = exp(-T/tau)``,

       P(floor(t/T) mod K = m) = (1 - q) q^m / (1 - q^K),   m = 0 .. K-1,

   and the mixing matrix in **exposure** indexing is

       M[i, j] = P(photon of exposure j is counted as exposure i)
               = (1 - q) q^((slot(i) - slot(j)) mod K) / (1 - q^K).

   ``M`` is circulant in slot order, hence doubly stochastic (rows and columns sum to 1): a
   time-uniform background stays uniform.  ``tau -> 0`` gives the identity.  The observed
   probabilities are ``p_obs = M p`` (:func:`make_model_crosstalk`).  Assumptions: the
   excitation probability per pulse is small (no saturation / no re-excitation of an excited
   molecule), windows exactly aligned with the pulses, mono-exponential decay.

   Note on the order: only the *cyclic* order of the pulses enters ``M``.  For the TCP every
   cyclic order of {0, 1, 2, centre} is mapped onto any other by a symmetry of the pattern
   (3-fold rotation / mirror), so the order changes the *orientation* of the bias map, not the
   CRB at the centre or averaged over a centred disk.

2. **Fluorophore flickering: sequential versus interleaved exposures.**  The emitter follows a
   two-state telegraph process (exponential on/off dwell times ``t_on``, ``t_off``, stationary
   start).  In a *sequential* scheme the dwell ``t_total`` is split into ``4 r`` consecutive
   blocks (``r`` pattern repetitions), so each exposure sees a different "on" time and the
   photon fractions are biased (SimuFLUX, Marin & Ries 2026, Fig. 2e / SI Fig. 8).  In the
   *interleaved* scheme each exposure receives its pulses every ``P`` during the whole dwell;
   the number of "on" pulses per exposure is computed exactly (:func:`interleaved_on_pulses`)
   and differs between exposures by at most one pulse per on-interval.
"""

import numpy as np

from . import photons

__all__ = ["PERIOD_NS", "WINDOW_NS", "crosstalk_matrix", "crosstalk_mc", "make_model_crosstalk",
           "telegraph_on_intervals", "sequential_on_times", "interleaved_on_pulses",
           "exposure_weights", "simulate_flicker_counts"]

PERIOD_NS = 50.0     # 20 MHz repetition period (Masullo et al. 2021)
WINDOW_NS = 12.5     # TCSPC window width = pulse separation T = P / 4


def _order(order, K=None):
    o = np.asarray(order, dtype=int)
    if o.ndim != 1 or o.size < 1:
        raise ValueError("order must be a 1-D sequence of exposure indices")
    if K is not None and o.size != K:
        raise ValueError("order has %d entries, expected K = %d" % (o.size, K))
    if sorted(o.tolist()) != list(range(o.size)):
        raise ValueError("order must be a permutation of 0..K-1, got %r" % (order,))
    return o


def crosstalk_matrix(tau, T=WINDOW_NS, order=(0, 1, 2, 3)):
    """Closed-form lifetime cross-talk matrix ``M`` (K x K), exposure indexing.

    Parameters
    ----------
    tau : float
        Fluorescence lifetime (ns), ``>= 0``.  ``tau = 0`` returns the identity.
    T : float
        Window width = pulse separation (ns), ``> 0``; the period is ``K T``.
    order : sequence of int
        ``order[s]`` is the exposure fired in slot ``s`` (pulse at ``s T`` within the period).
        Default ``(0, 1, 2, 3)``: the three TCP donuts then the centre (index 3).

    Returns
    -------
    M : ndarray (K, K), ``M[i, j]`` = probability that a photon excited by exposure ``j`` is
        counted in the window of exposure ``i``.  ``p_obs = M @ p``.
    """
    o = _order(order)
    K = o.size
    tau = float(tau)
    T = float(T)
    if not T > 0:
        raise ValueError("T must be > 0, got %r" % (T,))
    if not tau >= 0:
        raise ValueError("tau must be >= 0, got %r" % (tau,))
    m = np.arange(K)
    if tau == 0.0:
        c = (m == 0).astype(float)
    else:
        q = np.exp(-T / tau)
        c = (1.0 - q) * q ** m / (1.0 - q ** K)
    slot = np.empty(K, dtype=int)
    slot[o] = np.arange(K)                       # slot of each exposure
    shift = (slot[:, None] - slot[None, :]) % K  # (i, j) -> m
    return c[shift]


def crosstalk_mc(tau, T=WINDOW_NS, order=(0, 1, 2, 3), n=10 ** 6, rng=None):
    """Monte Carlo estimate of :func:`crosstalk_matrix` (independent route, for checks).

    For each exposure ``j`` simulates ``n`` photons: arrival time ``slot(j) T + t`` with
    ``t ~ Exp(tau)`` (ns), window ``floor(arrival / T) mod K``; returns the empirical
    ``M_hat[i, j]``.
    """
    o = _order(order)
    K = o.size
    if not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    Mh = np.zeros((K, K))
    for s, j in enumerate(o):
        t = s * T + rng.exponential(tau, size=int(n)) if tau > 0 else np.full(int(n), s * T)
        w = np.floor(t / T).astype(np.int64) % K
        cnt = np.bincount(w, minlength=K) / float(n)
        Mh[o, j] = cnt            # window w belongs to exposure o[w]
    return Mh


def make_model_crosstalk(centers, beam, tau, T=WINDOW_NS, order=(0, 1, 2, 3), sbr=None,
                         bg_per_exposure=None):
    """Model callable ``p_fn(r) -> (..., K)`` with lifetime cross-talk: ``p_obs = M p``.

    ``p`` is :func:`donutloc.photons.probabilities` (background included as there; since ``M``
    is doubly stochastic a uniform background term is unchanged).  The callable carries the
    attributes of :func:`photons.make_model` plus ``M, tau, T, order``.
    """
    base = photons.make_model(centers, beam, sbr=sbr, bg_per_exposure=bg_per_exposure)
    M = crosstalk_matrix(tau, T=T, order=order)
    if M.shape[0] != base.K:
        raise ValueError("order has %d entries but the pattern has K = %d exposures"
                         % (M.shape[0], base.K))

    def p_fn(r):
        return np.asarray(base(r), dtype=float) @ M.T

    for k in ("centers", "beam", "sbr", "bg_per_exposure", "K"):
        setattr(p_fn, k, getattr(base, k))
    p_fn.M = M
    p_fn.tau = float(tau)
    p_fn.T = float(T)
    p_fn.order = tuple(int(v) for v in order)
    return p_fn


# ---------------------------------------------------------------------------- flickering
def telegraph_on_intervals(n, t_total, t_on, t_off, rng=None, start="stationary"):
    """Sample ``n`` independent telegraph (on/off) trajectories on ``[0, t_total]`` (µs).

    Dwell times are exponential with means ``t_on`` (on) and ``t_off`` (off).  ``start``:
    ``"stationary"`` (on with probability ``t_on/(t_on+t_off)``), ``"on"`` or ``"off"``.

    Returns ``(starts, ends)``, arrays ``(n, S)`` of the on-intervals clipped to
    ``[0, t_total]``; padding / off intervals have ``start == end`` (zero length).
    """
    if not (t_total > 0 and t_on > 0 and t_off > 0):
        raise ValueError("t_total, t_on and t_off must be > 0")
    if not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    n = int(n)
    if start == "stationary":
        s0 = rng.random(n) < t_on / (t_on + t_off)
    elif start == "on":
        s0 = np.ones(n, dtype=bool)
    elif start == "off":
        s0 = np.zeros(n, dtype=bool)
    else:
        raise ValueError("start must be 'stationary', 'on' or 'off'")
    S = int(8 + 4 * t_total / min(t_on, t_off))
    while True:
        k = np.arange(S)
        on_k = s0[:, None] ^ (k[None, :] % 2 == 1)            # state of interval k
        mean = np.where(on_k, t_on, t_off)
        d = rng.exponential(1.0, size=(n, S)) * mean
        c = np.concatenate([np.zeros((n, 1)), np.cumsum(d, axis=1)], axis=1)
        if np.all(c[:, -1] >= t_total):
            break
        S *= 2
    a = np.minimum(c[:, :-1], t_total)
    b = np.minimum(c[:, 1:], t_total)
    b = np.where(on_k, b, a)
    return a, b


def sequential_on_times(starts, ends, t_total, reps, order=(0, 1, 2, 3)):
    """On-time (µs) of each exposure for a sequential scheme: ``[0, t_total]`` is split into
    ``K reps`` equal consecutive blocks, block ``k`` belonging to exposure ``order[k % K]``.

    Returns ``(n, K)``.  Without flickering every entry would be ``t_total / K``.
    """
    o = _order(order)
    K = o.size
    reps = int(reps)
    if reps < 1:
        raise ValueError("reps must be >= 1")
    edges = np.linspace(0.0, float(t_total), K * reps + 1)
    a = np.asarray(starts, float)[..., None]
    b = np.asarray(ends, float)[..., None]
    ov = np.clip(np.minimum(b, edges[1:]) - np.maximum(a, edges[:-1]), 0.0, None)  # (n, S, B)
    per_block = ov.sum(axis=1)                                                    # (n, B)
    out = np.zeros((per_block.shape[0], K))
    for blk in range(K * reps):
        out[:, o[blk % K]] += per_block[:, blk]
    return out


def interleaved_on_pulses(starts, ends, t_total, period=PERIOD_NS * 1e-3,
                          T=WINDOW_NS * 1e-3, order=(0, 1, 2, 3)):
    """Exact number of pulses of each exposure that hit the emitter while it is "on".

    Exposure ``order[s]`` is pulsed at times ``s T + k period`` (µs), ``k = 0 .. n_per - 1``,
    ``n_per = round(t_total / period)``.  The number of such pulses in ``[a, b)`` is
    ``G(b) - G(a)`` with ``G(t) = clip(ceil((t - s T)/period), 0, n_per)``.

    Returns ``(n, K)`` integer counts (without flickering each would be ``n_per``).
    """
    o = _order(order)
    K = o.size
    n_per = int(round(float(t_total) / float(period)))
    a = np.asarray(starts, float)
    b = np.asarray(ends, float)
    out = np.zeros((a.shape[0], K), dtype=np.int64)
    for s, j in enumerate(o):
        def G(t):
            return np.clip(np.ceil((t - s * T) / period), 0, n_per).astype(np.int64)
        out[:, j] = np.sum(G(b) - G(a), axis=1)
    return out


def exposure_weights(scheme, n, t_total=400.0, t_on=100.0, t_off=100.0, reps=1,
                     order=(0, 1, 2, 3), rng=None, period=PERIOD_NS * 1e-3,
                     T=WINDOW_NS * 1e-3):
    """Relative "on" exposure of each of the K exposures for ``n`` localizations.

    ``scheme``: ``"sequential"`` (``reps`` pattern repetitions), ``"interleaved"`` (p-MINFLUX
    pulses) or ``"none"`` (no flickering, constant brightness).  Returns ``w`` of shape
    ``(n, K)`` normalized so that ``E[w] = 1`` for every exposure: sequential on-times are
    divided by ``duty * t_total / K`` and interleaved on-pulse counts by ``duty * n_per``
    (``duty = t_on / (t_on + t_off)``, ``n_per = round(t_total / period)``); ``"none"``: 1.
    """
    K = _order(order).size
    if scheme == "none":
        return np.ones((int(n), K))
    a, b = telegraph_on_intervals(n, t_total, t_on, t_off, rng=rng)
    duty = t_on / float(t_on + t_off)
    if scheme == "sequential":
        w = sequential_on_times(a, b, t_total, reps, order=order)
        return w / (duty * float(t_total) / K)
    if scheme == "interleaved":
        n_per = int(round(float(t_total) / float(period)))
        w = interleaved_on_pulses(a, b, t_total, period=period, T=T, order=order)
        return w / (duty * n_per)
    raise ValueError("scheme must be 'sequential', 'interleaved' or 'none'")


def simulate_flicker_counts(p, w, n_mean=100.0, rng=None, fixed_N=None):
    """Photon counts with flickering weights ``w`` (from :func:`exposure_weights`).

    ``p`` : ``(K,)`` ideal probabilities at the true position.  The expected counts are
    ``n_mean * w[:, i] * p[i]`` (so ``E[N] = n_mean``).  Default: independent Poisson counts.
    ``fixed_N``: instead draw a multinomial with exactly ``fixed_N`` photons and probabilities
    ``w p / sum(w p)`` (rows with ``sum(w p) = 0`` get zero counts).  Returns ``(n, K)`` int64.
    """
    if not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    p = np.asarray(p, float)
    w = np.asarray(w, float)
    lam = n_mean * w * p[None, :]
    if fixed_N is None:
        return rng.poisson(lam).astype(np.int64)
    tot = lam.sum(axis=1, keepdims=True)
    q = np.where(tot > 0, lam / np.where(tot > 0, tot, 1.0), 0.0)
    out = np.zeros(lam.shape, dtype=np.int64)
    good = tot[:, 0] > 0
    q = q[good]
    q = q / q.sum(axis=1, keepdims=True)
    out[good] = rng.multinomial(int(fixed_N), q)
    return out
