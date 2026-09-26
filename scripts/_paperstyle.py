# -*- coding: utf-8 -*-
"""Shared publication style for every scripts/fig_*.py (see equipo/.../reports/r03-pi.md).

Importing this module puts ``src/`` on ``sys.path`` (via ``_paperconfig``) and applies the
rcParams. API: COLORS, CYCLE, CMAP_SEQ, CMAP_DIV, COL1, COL2, apply(), figure(), panel_label(),
savefig(), report(), write_summary(), parse_args(), cached(), config_hash(), fig_dir(),
summary_dir().

Quick mode (``--quick``, set by ``parse_args`` or ``set_quick``): figures go to
``paper/figures/quick/``, summaries to ``data/quick/`` and caches to ``data/mc/*_quick.npz``;
the final products are never touched by a quick run.

Cache invalidation: ``cached`` stores in every npz a key = sha256 (16 hex) of the canonical
parameters of ``_paperconfig`` (every upper-case constant except paths), of the optional
``params`` passed by the caller and of the quick flag; a cache whose key differs (or is missing)
is recomputed.  Code changes are NOT detected: the final regeneration uses ``--no-cache``.
"""
import argparse
import hashlib
import json
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import _paperconfig  # noqa: E402  (puts src/ on sys.path)

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ---------------------------------------------------------------- colours (Okabe-Ito)
COLORS = {
    "lg": "#0072B2",
    "vec": "#D55E00",
    "opp": "#009E73",
    "lin": "#CC79A7",
    "cam": "#000000",
    "crb": "#000000",
    "mle": "#0072B2",
    "lms": "#E69F00",
    "mlms": "#009E73",
    "honest": "#0072B2",
    "naive": "#D55E00",
    "sbr": "#56B4E9",
    "extra": "#F0E442",
}
# Okabe-Ito order without the leading yellow-ish entries that print poorly
CYCLE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442"]
CMAP_SEQ = "viridis"
CMAP_DIV = "RdBu_r"

# ---------------------------------------------------------------- sizes (inches)
COL1 = 3.4
COL2 = 7.0


def apply():
    """Set the shared rcParams (called on import)."""
    from cycler import cycler
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 0.6,
        "axes.prop_cycle": cycler(color=CYCLE),
        "lines.linewidth": 1.2,
        "lines.markersize": 3.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.6,
        "ytick.minor.size": 1.6,
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "savefig.dpi": 300,
        "figure.dpi": 150,
        "image.cmap": CMAP_SEQ,
    })


apply()


def figure(width="single", aspect=0.75, nrows=1, ncols=1, **subplots_kw):
    """Create a figure of single (COL1) or double (COL2) column width; height = width*aspect.

    ``width`` may also be a number in inches. Returns ``(fig, axes)`` as plt.subplots does.
    """
    if width == "single":
        w = COL1
    elif width == "double":
        w = COL2
    else:
        w = float(width)
    subplots_kw.setdefault("constrained_layout", True)
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, w * aspect), **subplots_kw)
    return fig, axes


def panel_label(ax, letter, x=-0.02, y=1.02, **kw):
    """Bold '(a)' at the top-left of ``ax`` in axes coordinates."""
    letter = str(letter).strip("()")
    kw.setdefault("fontsize", 9)
    kw.setdefault("fontweight", "bold")
    return ax.text(x, y, "(%s)" % letter, transform=ax.transAxes, ha="right", va="bottom", **kw)


_STATE = {"quick": False}


def set_quick(flag):
    """Switch quick mode on/off (separate output directories, see module docstring)."""
    _STATE["quick"] = bool(flag)


def fig_dir(quick=None):
    """paper/figures (final) or paper/figures/quick (quick mode)."""
    q = _STATE["quick"] if quick is None else quick
    return os.path.join(_paperconfig.FIGDIR, _paperconfig.QUICK_DIRNAME) if q else _paperconfig.FIGDIR


def summary_dir(quick=None):
    """data (final) or data/quick (quick mode)."""
    q = _STATE["quick"] if quick is None else quick
    return os.path.join(_paperconfig.DATA, _paperconfig.QUICK_DIRNAME) if q else _paperconfig.DATA


def savefig(fig, name):
    """Write <fig_dir()>/<name>.pdf, close the figure, print and return the path."""
    d = fig_dir()
    os.makedirs(d, exist_ok=True)
    base = name[:-4] if name.endswith(".pdf") else name
    path = os.path.join(d, base + ".pdf")
    fig.savefig(path)
    plt.close(fig)
    print("wrote %s" % path)
    return path


_SUMMARY = {}


def _jsonable(v):
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist()
    return v


def report(key, value, unit="", se=None):
    """Print ``NUMBER key = value [± se] unit`` and store it for write_summary()."""
    value = _jsonable(value)
    se = _jsonable(se)
    txt = "NUMBER %s = %s" % (key, _fmt(value))
    if se is not None:
        txt += " +- %s" % _fmt(se)
    if unit:
        txt += " %s" % unit
    print(txt)
    entry = {"value": value, "unit": unit}
    if se is not None:
        entry["se"] = se
    _SUMMARY[key] = entry
    return value


def _fmt(v):
    if isinstance(v, float):
        return "%.6g" % v
    return str(v)


def write_summary(fig_id):
    """Dump the numbers reported so far to <summary_dir()>/<fig_id>_summary.json (quick runs
    also get ``"quick": true``)."""
    d = summary_dir()
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "%s_summary.json" % fig_id)
    out = dict(_SUMMARY)
    if _STATE["quick"]:
        out["quick"] = True
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("wrote %s" % path)
    return path


def parse_args(description, argv=None):
    """argparse with --quick and --no-cache (also switches quick mode on via set_quick)."""
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--quick", action="store_true", help="reduced Monte Carlo / grids; outputs to "
                   "paper/figures/quick, data/quick and data/mc/*_quick.npz (final files untouched)")
    p.add_argument("--no-cache", dest="no_cache", action="store_true", help="recompute, ignoring data/mc cache")
    a = p.parse_args(argv)
    set_quick(a.quick)
    return a


_NOT_PARAMS = ("ROOT", "SRC", "DATA", "MC", "FIGDIR", "QUICK_DIRNAME")
_CACHE_KEY = "_cache_key"


def _canon(v):
    if isinstance(v, dict):
        return {str(k): _canon(x) for k, x in sorted(v.items(), key=lambda kv: str(kv[0]))}
    if isinstance(v, (list, tuple)):
        return [_canon(x) for x in v]
    if isinstance(v, np.ndarray):
        return [_canon(x) for x in v.tolist()]
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, (bool, int, float, str)) or v is None:
        return v
    return repr(v)


def config_params():
    """Upper-case constants of _paperconfig (paths excluded), canonicalized."""
    return {k: _canon(getattr(_paperconfig, k)) for k in sorted(dir(_paperconfig))
            if k.isupper() and k not in _NOT_PARAMS}


def config_hash(params=None, quick=False):
    """16-hex sha256 of the _paperconfig parameters, the caller's ``params`` and ``quick``."""
    blob = json.dumps({"config": config_params(), "params": _canon(params), "quick": bool(quick)},
                      sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def cached(name, fn, quick=False, force=False, params=None):
    """Load data/mc/<name>[_quick].npz unless ``force`` or its stored parameter hash differs
    from ``config_hash(params, quick)``; otherwise run fn() (-> dict) and save it with the hash."""
    os.makedirs(_paperconfig.MC, exist_ok=True)
    fname = name + ("_quick" if quick else "") + ".npz"
    path = os.path.join(_paperconfig.MC, fname)
    key = config_hash(params, quick)
    if os.path.exists(path) and not force:
        with np.load(path, allow_pickle=False) as z:
            stored = str(z[_CACHE_KEY].item()) if _CACHE_KEY in z.files else None
            if stored == key:
                out = {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files
                       if k != _CACHE_KEY}
                print("cache: loaded %s (key %s)" % (path, key))
                return out
        print("cache: %s is stale (key %s != %s), recomputing" % (path, stored, key))
    t0 = time.time()
    out = fn()
    if _CACHE_KEY in out:
        raise KeyError("%r is reserved" % _CACHE_KEY)
    arrays = {k: np.asarray(v) for k, v in out.items()}
    arrays[_CACHE_KEY] = np.asarray(key)
    np.savez(path, **arrays)
    print("cache: computed %s in %.1f s and saved (key %s)" % (fname, time.time() - t0, key))
    return out
