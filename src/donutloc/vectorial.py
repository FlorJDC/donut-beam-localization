# -*- coding: utf-8 -*-
"""Vectorial (Richards-Wolf / Debye) focal field of a vortex-phase donut, in nm.

Own implementation (no third-party code) of the aplanatic focusing integral of
Caprile2022 Eq. 1 (docs/literature/B_donut_optics.md section 2.1)::

    E(rho, phi, z) = (1/2pi) int_0^{2pi} int_0^alpha E_0(theta, phi')
                     exp{i k [z cos(theta) + rho sin(theta) cos(phi' - phi)]}
                     sqrt(cos theta) sin(theta) dtheta dphi'

with k = 2 pi n / lambda, sin(alpha) = NA/n, and the constant prefactor -i k f exp(-i k f)
dropped (fields are returned **unnormalized**; only ratios matter). The refracted amplitude is
``E_0 = (E_i . phi') phi' + (E_i . rho') theta`` with
``rho' = (cos phi', sin phi', 0)``, ``phi' = (-sin phi', cos phi', 0)``,
``theta = (cos th cos phi', cos th sin phi', -sin th)``.

Pupil field (sine condition rho' = f sin theta, filling F = w0/h, so f sin th / w0 =
sin th / (F sin alpha))::

    E_i = a(theta) exp(i l phi') e,     a(theta) = exp(-(sin th / (F sin alpha))^2)

with topological charge ``l = charge`` and polarization vector ``e``:

- ``"circular"``: ``e = (x + i s y)/sqrt 2``, ``s = handedness = +-1``. The "same hand"
  (perfect zero) case is ``l s = +1`` (note B section 2.2 rule); with the default ``charge=1``,
  ``handedness=+1`` is the correct hand and ``handedness=-1`` the opposite hand.
- ``"linear"``: ``e = (cos g, sin g, 0)`` with ``g = pol_angle`` (rad).

Bessel reduction (derived here; generalizes note B section 2.2)
---------------------------------------------------------------
For ``e = (x + i s y)/sqrt 2``: ``e.rho' = exp(i s phi')/sqrt 2`` and
``e.phi' = i s exp(i s phi')/sqrt 2``, so ``E_0 = a exp(i m phi') (i s phi' + theta)/sqrt 2``
with ``m = l + s``. Writing cos/sin phi' as exponentials, the Cartesian components of
``(i s phi' + theta)`` are::

    x: [(cos th - s) e^{+i phi'} + (cos th + s) e^{-i phi'}] / 2
    y: (i/2) [(s - cos th) e^{+i phi'} + (s + cos th) e^{-i phi'}]
    z: -sin th

and the identity int_0^{2pi} e^{i p phi'} e^{i u cos(phi'-phi)} dphi' = 2 pi i^p e^{i p phi}
J_p(u) (Caprile2022 Eq. 12) gives, with g(th) = a sqrt(cos th) sin th e^{i k z cos th},
u = k rho sin th::

    P = int g (cos th - s) J_{m+1}(u) dth,  Q = int g (cos th + s) J_{m-1}(u) dth,
    R = int g sin th J_m(u) dth

    Ex = [ i^{m+1} e^{i(m+1)phi} P + i^{m-1} e^{i(m-1)phi} Q ] / (2 sqrt 2)
    Ey = (i/(2 sqrt 2)) [ -i^{m+1} e^{i(m+1)phi} P + i^{m-1} e^{i(m-1)phi} Q ]
    Ez = -i^m e^{i m phi} R / sqrt 2

    I = |Ex|^2 + |Ey|^2 + |Ez|^2 = (|P|^2 + |Q|^2)/4 + |R|^2/2     (no phi dependence)

* ``l = 1, s = +1`` (m = 2, correct hand): P = -I_B, Q = I_A, R = I_C of note B section 2.2,
  i.e. ``Ex = (i/(2 sqrt 2))(I_A e^{i phi} + I_B e^{3 i phi})``,
  ``Ey = -(1/(2 sqrt 2))(I_A e^{i phi} - I_B e^{3 i phi})``, ``Ez = I_C e^{2 i phi}/sqrt 2``
  (identical to note B; **not** the E_y of the printed Caprile2022 Eq. 15). Orders J1, J3
  (transverse) and J2 (axial): all vanish at rho = 0, so I(0, z) = 0 for every z.
* ``l = 1, s = -1`` (m = 0, opposite hand): with J_{-1} = -J_1,
  ``Ex = (i/(2 sqrt 2)) (I_A e^{i phi} - I_B' e^{-i phi})``,
  ``Ey = -(1/(2 sqrt 2)) (I_A e^{i phi} + I_B' e^{-i phi})``,
  ``Ez = -I_C'/sqrt 2`` with ``I_A = int g (1 + cos th) J1``, ``I_B' = int g (1 - cos th) J1``,
  ``I_C' = int g sin th J0``. Transverse components are pure J1 (zero on axis); the axial one is
  J0, maximal on axis: ``I(0) = |I_C'(0)|^2 / 2``.

Linear polarization is not rotationally symmetric. It is computed exactly as the superposition
``(cos g, sin g) = e^{-i g}(x + i y)/2 + e^{i g}(x - i y)/2``, i.e.
``E_lin = [e^{-i g} E(s=+1) + e^{i g} E(s=-1)] / sqrt 2`` (fields add, not intensities), which
is again a sum of 1-D Bessel integrals with explicit phi harmonics. The direct 2-D pupil
integral in (theta, phi') is available as ``method="2d"`` for every polarization and is the
independent check (tests/test_vectorial.py).

Quadrature: trapezoid in theta on ``n_theta`` points (default 801, note B section 6.2); for
``method="2d"`` a periodic trapezoid (spectrally accurate) with ``n_phi`` points in phi'.
"""

import numpy as np
from scipy.special import jv
from scipy.optimize import minimize_scalar, minimize
from scipy.interpolate import RegularGridInterpolator

_LN2 = np.log(2.0)
_POLS = ("circular", "linear")
_OPT_KEYS = ("z", "wavelength", "NA", "n", "filling", "polarization", "handedness", "charge",
             "n_theta", "pol_angle", "method", "n_phi")


# ----------------------------------------------------------------------------- parameters
def _params(wavelength=640.0, NA=1.4, n=1.518, filling=5.0 / 3.0, n_theta=801):
    wavelength = float(wavelength)
    NA = float(NA)
    n = float(n)
    filling = float(filling)
    n_theta = int(n_theta)
    if not wavelength > 0:
        raise ValueError("wavelength must be > 0, got %r" % (wavelength,))
    if not (NA > 0 and NA < n):
        raise ValueError("need 0 < NA < n, got NA=%r, n=%r" % (NA, n))
    if not filling > 0:
        raise ValueError("filling must be > 0, got %r" % (filling,))
    if n_theta < 3:
        raise ValueError("n_theta must be >= 3, got %r" % (n_theta,))
    k = 2.0 * np.pi * n / wavelength
    alpha = np.arcsin(NA / n)
    th = np.linspace(0.0, alpha, n_theta)
    w = np.full(n_theta, th[1] - th[0])
    w[0] *= 0.5
    w[-1] *= 0.5
    a = np.exp(-(np.sin(th) / (filling * NA / n)) ** 2)
    return k, th, w, a


def _check_pol(polarization, handedness, charge):
    if polarization not in _POLS:
        raise ValueError("polarization must be one of %s, got %r" % (_POLS, polarization))
    if int(handedness) not in (1, -1) or float(handedness) != int(handedness):
        raise ValueError("handedness must be +1 or -1, got %r" % (handedness,))
    if float(charge) != int(charge):
        raise ValueError("charge must be an integer, got %r" % (charge,))
    return int(handedness), int(charge)


def _split_opts(opt):
    bad = set(opt) - set(_OPT_KEYS)
    if bad:
        raise TypeError("unexpected option(s): %s" % sorted(bad))
    return dict(opt)


# ----------------------------------------------------------------------------- Bessel route
def _pqr(rho, z, k, th, w, a, s, l, chunk=1024):
    """Radial integrals P, Q, R (module docstring) for 1-D ``rho`` (nm). Complex arrays."""
    rho = np.asarray(rho, dtype=float).ravel()
    m = l + s
    st, ct = np.sin(th), np.cos(th)
    g = a * np.sqrt(ct) * st * np.exp(1j * k * float(z) * ct) * w  # quadrature weights folded
    gp = g * (ct - s)
    gq = g * (ct + s)
    gr = g * st
    P = np.empty(rho.size, complex)
    Q = np.empty(rho.size, complex)
    R = np.empty(rho.size, complex)
    for i0 in range(0, rho.size, chunk):
        u = k * rho[i0:i0 + chunk, None] * st[None, :]
        P[i0:i0 + chunk] = jv(m + 1, u) @ gp
        Q[i0:i0 + chunk] = jv(m - 1, u) @ gq
        R[i0:i0 + chunk] = jv(m, u) @ gr
    return P, Q, R


def _assemble(P, Q, R, phi, s, l):
    m = l + s
    ep = (1j ** (m + 1)) * np.exp(1j * (m + 1) * phi) * P
    eq = (1j ** (m - 1)) * np.exp(1j * (m - 1) * phi) * Q
    Ex = (ep + eq) / (2.0 * np.sqrt(2.0))
    Ey = 1j * (-ep + eq) / (2.0 * np.sqrt(2.0))
    Ez = -(1j ** m) * np.exp(1j * m * phi) * R / np.sqrt(2.0)
    return Ex, Ey, Ez


def _field_bessel(rho, phi, z, k, th, w, a, pol, s, l, gamma):
    shape = np.broadcast(rho, phi).shape
    rho = np.broadcast_to(np.asarray(rho, float), shape).ravel()
    phi = np.broadcast_to(np.asarray(phi, float), shape).ravel()
    ur, inv = np.unique(rho, return_inverse=True)
    if pol == "circular":
        P, Q, R = _pqr(ur, z, k, th, w, a, s, l)
        out = _assemble(P[inv], Q[inv], R[inv], phi, s, l)
    else:
        Pp, Qp, Rp = _pqr(ur, z, k, th, w, a, +1, l)
        Pm, Qm, Rm = _pqr(ur, z, k, th, w, a, -1, l)
        Fp = _assemble(Pp[inv], Qp[inv], Rp[inv], phi, +1, l)
        Fm = _assemble(Pm[inv], Qm[inv], Rm[inv], phi, -1, l)
        cp = np.exp(-1j * gamma) / np.sqrt(2.0)
        cm = np.exp(1j * gamma) / np.sqrt(2.0)
        out = tuple(cp * fp + cm * fm for fp, fm in zip(Fp, Fm))
    return tuple(c.reshape(shape) for c in out)


# ----------------------------------------------------------------------------- direct 2-D route
def _field_2d(rho, phi, z, k, th, w, a, pol, s, l, gamma, n_phi):
    """Direct 2-D quadrature of Caprile2022 Eq. 1 (independent check; slow)."""
    shape = np.broadcast(rho, phi).shape
    rho = np.broadcast_to(np.asarray(rho, float), shape).ravel()
    phi = np.broadcast_to(np.asarray(phi, float), shape).ravel()
    n_phi = int(n_phi)
    pp = 2.0 * np.pi * np.arange(n_phi) / n_phi
    dpp = 2.0 * np.pi / n_phi
    T, PP = np.meshgrid(th, pp, indexing="ij")          # (nt, nphi)
    st, ct = np.sin(T), np.cos(T)
    cp, sp = np.cos(PP), np.sin(PP)
    if pol == "circular":
        ex, ey = 1.0 / np.sqrt(2.0), 1j * s / np.sqrt(2.0)
    else:
        ex, ey = np.cos(gamma), np.sin(gamma)
    amp = a[:, None] * np.exp(1j * l * PP)
    e_rho = amp * (ex * cp + ey * sp)                   # E_i . rho'
    e_phi = amp * (-ex * sp + ey * cp)                  # E_i . phi'
    E0x = e_phi * (-sp) + e_rho * ct * cp
    E0y = e_phi * cp + e_rho * ct * sp
    E0z = -e_rho * st
    wt = (w[:, None] * dpp) * np.sqrt(ct) * st * np.exp(1j * k * float(z) * ct) / (2.0 * np.pi)
    out = [np.empty(rho.size, complex) for _ in range(3)]
    for i in range(rho.size):
        ph = np.exp(1j * k * rho[i] * st * np.cos(PP - phi[i])) * wt
        out[0][i] = np.sum(E0x * ph)
        out[1][i] = np.sum(E0y * ph)
        out[2][i] = np.sum(E0z * ph)
    return tuple(c.reshape(shape) for c in out)


# ----------------------------------------------------------------------------- public API
def focal_field(rho, phi, z=0.0, wavelength=640.0, NA=1.4, n=1.518, filling=5.0 / 3.0,
                polarization="circular", handedness=+1, charge=1, n_theta=801, pol_angle=0.0,
                method="bessel", n_phi=200):
    """Focal field ``(Ex, Ey, Ez)`` (complex, unnormalized) at polar points ``(rho, phi)`` (nm, rad)
    in the plane ``z`` (nm; z = 0 is the focus).

    Parameters
    ----------
    wavelength : vacuum wavelength (nm). NA, n : numerical aperture and medium index.
    filling : F = w0/h (Gaussian 1/e amplitude radius over aperture radius); large F -> uniform.
    polarization : {"circular", "linear"}.
    handedness : +1 or -1, sign s of ``(x + i s y)/sqrt 2`` (circular only). Correct hand: s = charge.
    charge : integer topological charge l of the vortex ``exp(i l phi')`` (0: no vortex).
    pol_angle : angle of the linear polarization with x (rad), linear only.
    method : "bessel" (1-D integrals, default) or "2d" (direct pupil quadrature, ``n_phi``
        points in phi'; slow, for checks).
    """
    k, th, w, a = _params(wavelength, NA, n, filling, n_theta)
    s, l = _check_pol(polarization, handedness, charge)
    if method == "bessel":
        return _field_bessel(rho, phi, z, k, th, w, a, polarization, s, l, float(pol_angle))
    if method == "2d":
        return _field_2d(rho, phi, z, k, th, w, a, polarization, s, l, float(pol_angle), n_phi)
    raise ValueError("method must be 'bessel' or '2d', got %r" % (method,))


def intensity(x, y, z=0.0, **opt):
    """Total intensity |Ex|^2 + |Ey|^2 + |Ez|^2 (unnormalized) at Cartesian ``(x, y)`` (nm),
    vectorized with broadcasting. ``opt`` as in :func:`focal_field`."""
    _split_opts(opt)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    Ex, Ey, Ez = focal_field(np.hypot(x, y), np.arctan2(y, x), z=z, **opt)
    return (Ex * Ex.conjugate()).real + (Ey * Ey.conjugate()).real + (Ez * Ez.conjugate()).real


def _radial_I(rho, z=0.0, **opt):
    """Rotationally symmetric intensity (P, Q, R route) on 1-D rho; circular only."""
    pol = opt.get("polarization", "circular")
    if pol != "circular":
        raise ValueError("radial_profile is only defined for rotationally symmetric cases "
                         "(polarization='circular'); got %r. Use intensity(x, y) instead." % (pol,))
    if opt.get("method", "bessel") != "bessel":
        raise ValueError("radial profile uses method='bessel'")
    k, th, w, a = _params(opt.get("wavelength", 640.0), opt.get("NA", 1.4),
                          opt.get("n", 1.518), opt.get("filling", 5.0 / 3.0),
                          opt.get("n_theta", 801))
    s, l = _check_pol(pol, opt.get("handedness", 1), opt.get("charge", 1))
    rho = np.asarray(rho, dtype=float)
    P, Q, R = _pqr(rho.ravel(), z, k, th, w, a, s, l)
    I = 0.25 * (np.abs(P) ** 2 + np.abs(Q) ** 2) + 0.5 * np.abs(R) ** 2
    return I.reshape(rho.shape)


def radial_profile(rho, **opt):
    """Radial intensity I(rho) (unnormalized) for the rotationally symmetric cases
    (``polarization="circular"``, either handedness). Raises ``ValueError`` for "linear"."""
    _split_opts(opt)
    z = opt.pop("z", 0.0)
    return _radial_I(rho, z=z, **opt)


def _is_circ(opt):
    return opt.get("polarization", "circular") == "circular"


_MAX_CACHE = {}


def _max_search(opt, rho_max=1000.0):
    """Return (I_max, rho_at_max, phi_at_max) of the total intensity in the plane z (memoized
    on the option values)."""
    try:
        key = (tuple(sorted((k, float(val) if not isinstance(val, str) else val)
                            for k, val in opt.items() if k != "method")), float(rho_max))
    except (TypeError, ValueError):
        key = None
    if key is not None and key in _MAX_CACHE:
        return _MAX_CACHE[key]
    res = _max_search_raw(opt, rho_max)
    if key is not None:
        if len(_MAX_CACHE) > 256:
            _MAX_CACHE.clear()
        _MAX_CACHE[key] = res
    return res


def _max_search_raw(opt, rho_max):
    o = dict(opt)
    o.pop("method", None)
    if _is_circ(o):
        z = o.pop("z", 0.0)
        r = np.arange(0.0, rho_max + 1e-9, 4.0)
        I = _radial_I(r, z=z, **o)
        i = int(np.argmax(I))
        if i == 0:
            return float(I[0]), 0.0, 0.0
        lo, hi = r[max(i - 1, 0)], r[min(i + 1, r.size - 1)]
        res = minimize_scalar(lambda t: -float(_radial_I(np.array([t]), z=z, **o)[0]),
                              bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
        return max(-float(res.fun), float(I[i])), float(res.x), 0.0
    r = np.arange(0.0, rho_max + 1e-9, 4.0)
    ph = np.linspace(0.0, 2.0 * np.pi, 145)[:-1]
    R, PH = np.meshgrid(r, ph, indexing="ij")
    I = intensity(R * np.cos(PH), R * np.sin(PH), **o)
    i, j = np.unravel_index(int(np.argmax(I)), I.shape)
    x0 = np.array([R[i, j] * np.cos(PH[i, j]), R[i, j] * np.sin(PH[i, j])])
    res = minimize(lambda p: -float(intensity(p[0], p[1], **o)), x0, method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-14 * float(I[i, j]), "maxiter": 2000})
    xb = res.x if -res.fun >= I[i, j] else x0
    return max(-float(res.fun), float(I[i, j])), float(np.hypot(*xb)), float(np.arctan2(xb[1], xb[0]))


def _I0(opt):
    o = dict(opt)
    o.pop("method", None)
    return float(intensity(0.0, 0.0, **o))


def zero_depth(**opt):
    """``I(0)/I_max`` in the plane ``z`` (default focus). I_max: radial search (circular) or 2-D
    grid search + Nelder-Mead refinement (linear)."""
    _split_opts(opt)
    Imax = _max_search(opt)[0]
    return _I0(opt) / Imax


def peak_to_peak_diameter(**opt):
    """Peak-to-peak ring diameter (nm), ``2 rho_peak``.

    Circular: ``rho_peak`` = radius of the maximum of I(rho). Linear: radius of the maximum of the
    azimuthally averaged intensity (the averaged profile of a non-symmetric focus).
    """
    _split_opts(opt)
    if _is_circ(opt):
        return 2.0 * _max_search(opt)[1]
    o = dict(opt)
    o.pop("method", None)
    ph = np.linspace(0.0, 2.0 * np.pi, 73)[:-1]

    def avg(r):
        return float(np.mean(intensity(r * np.cos(ph), r * np.sin(ph), **o)))

    r = np.arange(0.0, 1000.0, 2.0)
    vals = np.array([avg(t) for t in r])
    i = int(np.argmax(vals))
    lo, hi = r[max(i - 1, 0)], r[min(i + 1, r.size - 1)]
    res = minimize_scalar(lambda t: -avg(t), bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-5})
    return 2.0 * float(res.x)


def zero_curvature(**opt):
    """Coefficient c (nm^-2) of ``I(rho)/I_max ~= I(0)/I_max + c rho^2`` near the center.

    Least-squares fit of ``I/I_max = c0 + c rho^2 + d rho^4`` on rho = 1..10 nm (azimuthal average
    over 36 directions for "linear"). For the correct-hand circular case this equals the
    analytic small-argument limit ``(k^2/16) |int g (1 + cos th) sin th dth|^2 / I_max``.
    """
    _split_opts(opt)
    Imax = _max_search(opt)[0]
    o = dict(opt)
    o.pop("method", None)
    r = np.linspace(1.0, 10.0, 10)
    if _is_circ(o):
        z = o.pop("z", 0.0)
        I = _radial_I(r, z=z, **o)
    else:
        ph = np.linspace(0.0, 2.0 * np.pi, 37)[:-1]
        R, PH = np.meshgrid(r, ph, indexing="ij")
        I = intensity(R * np.cos(PH), R * np.sin(PH), **o).mean(axis=1)
    A = np.stack([np.ones_like(r), r ** 2, r ** 4], axis=1)
    coef = np.linalg.lstsq(A, I / Imax, rcond=None)[0]
    return float(coef[1])


def lg_equivalent_fwhm(beam_opts=None, match="curvature"):
    """``fwhm`` (nm) of the project LG donut (beams.lg_donut, Balzarotti Eq. S17) equivalent to the
    vectorial donut defined by ``beam_opts`` (dict of focal-field options).

    - ``match="curvature"``: equal curvature at the zero, ``4 e ln2 / fwhm^2 = c``
      (c = :func:`zero_curvature`), so ``fwhm = sqrt(4 e ln2 / c)``;
    - ``match="diameter"``: equal peak-to-peak diameter, ``D_pp = fwhm / sqrt(ln2)``, so
      ``fwhm = D_pp sqrt(ln2)``.
    """
    o = dict(beam_opts or {})
    if match == "curvature":
        c = zero_curvature(**o)
        if not c > 0:
            raise ValueError("non-positive zero curvature %r" % (c,))
        return float(np.sqrt(4.0 * np.e * _LN2 / c))
    if match == "diameter":
        return float(peak_to_peak_diameter(**o) * np.sqrt(_LN2))
    raise ValueError("match must be 'curvature' or 'diameter', got %r" % (match,))


def make_vectorial_beam(rho_max=1500.0, d_rho=1.0, eps=0.0, mode="interp", grid_step=5.0,
                        fwhm_eps=None, **opt):
    """Vectorized beam ``f(x, y)`` with peak (ring maximum) 1, compatible with
    ``photons.make_model``.

    Parameters
    ----------
    rho_max : float
        Tabulation radius (nm); the beam is 0 for rho > rho_max (``mode="interp"``).
    d_rho : float
        Radial table step (nm), circular. Interpolation is linear in ``s = rho^2`` (np.interp),
        which keeps the exact ``I ~ c rho^2`` behaviour below the first node.
    eps : float
        Optional residual zero (>= 0): ``+ eps * exp(-4 ln2 rho^2/fwhm_eps^2)`` after the peak
        normalization (like ``beams.lg_donut`` with zero_model="gaussian"). ``fwhm_eps`` defaults
        to the diameter-matched LG fwhm (:func:`lg_equivalent_fwhm`, match="diameter").
    mode : "interp" (tabulated) or "exact" (evaluates the Bessel integrals at every call;
        slow for big arrays, exact for CRB computations; no rho_max cut).
    grid_step : float
        Cartesian grid step (nm, <= 5 recommended) of the 2-D table for "linear" (interpolated with
        ``scipy.interpolate.RegularGridInterpolator``, linear, 0 outside the square
        ``|x|, |y| <= rho_max``).
    opt : focal-field options (``z``, ``wavelength``, ``NA``, ``n``, ``filling``,
        ``polarization``, ``handedness``, ``charge``, ``n_theta``, ``pol_angle``).

    The returned function carries attributes ``kind="vectorial"``, ``opts``, ``I_max``, ``eps``,
    ``fwhm_eps``, ``mode``.
    """
    _split_opts(opt)
    opt.pop("method", None)
    eps = float(eps)
    if eps < 0:
        raise ValueError("eps must be >= 0, got %r" % (eps,))
    if mode not in ("interp", "exact"):
        raise ValueError("mode must be 'interp' or 'exact', got %r" % (mode,))
    rho_max = float(rho_max)
    if not rho_max > 0:
        raise ValueError("rho_max must be > 0")
    Imax = _max_search(opt)[0]
    if eps > 0:
        fe = lg_equivalent_fwhm(opt, match="diameter") if fwhm_eps is None else float(fwhm_eps)
        ae = 4.0 * _LN2 / fe ** 2
    else:
        fe, ae = fwhm_eps, 0.0
    circ = _is_circ(opt)

    if mode == "exact":
        def base(x, y):
            return intensity(x, y, **opt) / Imax
    elif circ:
        if not float(d_rho) > 0:
            raise ValueError("d_rho must be > 0")
        z = opt.get("z", 0.0)
        o = dict(opt)
        o.pop("z", None)
        r = np.arange(0.0, rho_max + 0.5 * float(d_rho), float(d_rho))
        s_tab = r * r
        I_tab = _radial_I(r, z=z, **o) / Imax

        def base(x, y):
            s2 = np.asarray(x, float) ** 2 + np.asarray(y, float) ** 2
            return np.interp(s2, s_tab, I_tab, right=0.0)
    else:
        if not float(grid_step) > 0:
            raise ValueError("grid_step must be > 0")
        g = np.arange(-rho_max, rho_max + 0.5 * grid_step, float(grid_step))
        X, Y = np.meshgrid(g, g, indexing="ij")
        T = intensity(X, Y, **opt) / Imax
        rgi = RegularGridInterpolator((g, g), T, method="linear", bounds_error=False,
                                      fill_value=0.0)

        def base(x, y):
            x = np.asarray(x, float)
            y = np.asarray(y, float)
            sh = np.broadcast(x, y).shape
            pts = np.stack([np.broadcast_to(x, sh).ravel(), np.broadcast_to(y, sh).ravel()], -1)
            return rgi(pts).reshape(sh)

    def f(x, y):
        out = base(x, y)
        if eps > 0:
            out = out + eps * np.exp(-ae * (np.asarray(x, float) ** 2 + np.asarray(y, float) ** 2))
        return out

    f.kind = "vectorial"
    f.opts = dict(opt)
    f.I_max = Imax
    f.eps = eps
    f.fwhm_eps = fe
    f.mode = mode
    return f


def compare_crb_vectorial_vs_lg(L_list=(50.0, 100.0, 150.0), N=100, sbr=None, **opt):
    """Central-limit CRB (fisher.crb_limit, nm) of the TCP with the vectorial donut vs the project
    LG donut matched by curvature and by diameter (and the default LG fwhm = 300 nm).

    Uses only public APIs: ``patterns.tcp_centers``, ``photons.make_model``, ``fisher.crb_limit``
    (and ``beams.make_beam`` for the LG). The vectorial beam is evaluated exactly
    (``make_vectorial_beam(mode="exact")``) to avoid interpolation error in the derivatives.

    Returns a dict of numpy arrays keyed ``L, crb_vectorial, crb_lg_curvature, crb_lg_diameter,
    crb_lg_300`` plus floats ``fwhm_curvature, fwhm_diameter, zero_depth`` and ``opts``.
    """
    from . import patterns, photons, fisher, beams
    _split_opts(opt)
    vb = make_vectorial_beam(mode="exact", **opt)
    f_c = lg_equivalent_fwhm(opt, match="curvature")
    f_d = lg_equivalent_fwhm(opt, match="diameter")
    lgs = {"crb_lg_curvature": beams.make_beam("donut", fwhm=f_c),
           "crb_lg_diameter": beams.make_beam("donut", fwhm=f_d),
           "crb_lg_300": beams.make_beam("donut", fwhm=300.0)}
    L_arr = np.asarray(L_list, dtype=float)
    out = {"L": L_arr, "crb_vectorial": np.empty(L_arr.size)}
    for key in lgs:
        out[key] = np.empty(L_arr.size)
    for i, L in enumerate(L_arr):
        c = patterns.tcp_centers(L)
        out["crb_vectorial"][i] = fisher.crb_limit(photons.make_model(c, vb, sbr=sbr), N)
        for key, b in lgs.items():
            out[key][i] = fisher.crb_limit(photons.make_model(c, b, sbr=sbr), N)
    out["fwhm_curvature"] = f_c
    out["fwhm_diameter"] = f_d
    out["zero_depth"] = zero_depth(**opt)
    out["N"] = N
    out["sbr"] = sbr
    out["opts"] = dict(opt)
    return out
