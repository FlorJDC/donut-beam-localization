# -*- coding: utf-8 -*-
"""Exposure-pattern geometry (positions of the beam zeros), in nm.

Convention (project-wide): the TCP (Balzarotti2017 Eq. S24) is 3 donuts on a circle of
**diameter L** at angles ``rotation + 2 pi k/3`` (default ``rotation = pi/2``, vertex up) plus
the central donut **last** (index 3). Balzarotti numbers the center 0; our index ``k`` for
k = 0, 1, 2 corresponds to Balzarotti's k + 1, and our index 3 to Balzarotti's 0.
"""

import numpy as np


def polygon_centers(L, M, rotation=np.pi / 2, center=True):
    """Centers of ``M`` exposures on a circle of diameter ``L`` (nm), angles
    ``rotation + 2 pi k/M``, plus the origin appended last if ``center``.

    Returns an array of shape ``(M + 1, 2)`` (or ``(M, 2)`` if ``center=False``).
    """
    M = int(M)
    if M < 1:
        raise ValueError("M must be >= 1, got %r" % (M,))
    L = float(L)
    if L < 0:
        raise ValueError("L must be >= 0, got %r" % (L,))
    ang = rotation + 2 * np.pi * np.arange(M) / M
    ring = np.stack([0.5 * L * np.cos(ang), 0.5 * L * np.sin(ang)], axis=1)
    if center:
        return np.vstack([ring, [[0.0, 0.0]]])
    return ring


def tcp_centers(L, rotation=np.pi / 2, center=True):
    """Triangular coordinate pattern (Balzarotti2017 Eq. S24) with L = diameter.

    Returns ``(4, 2)``: three peripheral zeros at ``rotation + 2 pi k/3`` then the center
    (or ``(3, 2)`` if ``center=False``). With the defaults this is bit-identical to ``_tcp(L)``
    of ``tests/test_acceptance.py``.
    """
    return polygon_centers(L, 3, rotation=rotation, center=center)


def perturb_centers(centers, displacement, rng=None, mode="random_direction"):
    """Misalign a pattern (Masullo-type model): return a displaced copy of ``centers``.

    Parameters
    ----------
    centers : array_like, shape (K, 2)
    displacement : float, array (K,), or array (K, 2)
        - scalar or shape ``(K,)``: magnitude (nm) of the shift of each center; each center is
          moved by that distance in an independent, uniformly random direction;
        - shape ``(K, 2)``: explicit displacement vectors (``rng`` and ``mode`` unused).
    rng : numpy.random.Generator, int, or None
        Generator or seed; ``None`` uses a fresh unseeded ``default_rng()``.
    mode : {"random_direction"}
    """
    c = np.array(centers, dtype=float)
    if c.ndim != 2 or c.shape[1] != 2:
        raise ValueError("centers must have shape (K, 2), got %r" % (c.shape,))
    K = c.shape[0]
    d = np.asarray(displacement, dtype=float)
    if d.shape == (K, 2):
        return c + d
    if mode != "random_direction":
        raise ValueError("mode must be 'random_direction', got %r" % (mode,))
    if d.ndim == 0:
        d = np.full(K, float(d))
    if d.shape != (K,):
        raise ValueError("displacement must be scalar, (K,) or (K, 2); got %r" % (d.shape,))
    if not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    phi = rng.uniform(0.0, 2 * np.pi, size=K)
    return c + d[:, None] * np.stack([np.cos(phi), np.sin(phi)], axis=1)
