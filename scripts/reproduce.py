# -*- coding: utf-8 -*-
"""Reproduce every product of the paper companion, in order (cross-platform).

Steps (each one must succeed before the next starts):

1. unit tests                ``python -m unittest discover -s tests``
2. (optional, ``--with-sweep``) iterative-MINFLUX sweep ``scripts/run_iterative_sweep.py``
   (otherwise the committed ``data/iterative_sweep.json`` is used)
3. all figures               ``scripts/make_all_figures.py --no-cache``
4. paper numbers             ``scripts/compute_paper_numbers.py``
5. provenance check          ``scripts/check_provenance.py``
6. acceptance test           ``python -m unittest tests.test_acceptance -v``
7. (optional) LaTeX          ``tectonic paper/main.tex`` if ``tectonic`` is on PATH (or given by
   ``--tectonic PATH`` / the ``TECTONIC`` environment variable); skipped otherwise.

``--quick`` runs steps 3-4 with reduced Monte Carlo into the separate quick outputs
(``paper/figures/quick/``, ``data/quick/``, ``paper/generated/quick/``); the final files are not
touched, so steps 5-7 then check the committed products.  ``--skip-tests`` skips step 1.
The full run takes about 25 minutes on a laptop (8 cores) (Monte Carlo of figures 5 and 8 and of
the misalignment numbers); ``--quick`` a few minutes.

Usage:  python scripts/reproduce.py [--quick] [--with-sweep] [--skip-tests] [--no-latex]
                                    [--tectonic PATH]
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def plan(quick=False, with_sweep=False, skip_tests=False, latex=None):
    """List of (label, argv) steps; ``latex`` is the tectonic executable or None."""
    py = sys.executable
    q = ["--quick"] if quick else []
    steps = []
    if not skip_tests:
        steps.append(("unit tests", [py, "-m", "unittest", "discover", "-s", "tests"]))
    if with_sweep:
        steps.append(("iterative sweep", [py, os.path.join("scripts", "run_iterative_sweep.py")] + q))
    steps += [
        ("figures", [py, os.path.join("scripts", "make_all_figures.py"), "--no-cache"] + q),
        ("paper numbers", [py, os.path.join("scripts", "compute_paper_numbers.py")] + q),
        ("provenance check", [py, os.path.join("scripts", "check_provenance.py")]),
        ("acceptance test", [py, "-m", "unittest", "tests.test_acceptance", "-v"]),
    ]
    if latex:
        steps.append(("LaTeX (tectonic)", [latex, os.path.join("paper", "main.tex")]))
    return steps


def find_tectonic(explicit=None):
    for cand in (explicit, os.environ.get("TECTONIC"), shutil.which("tectonic")):
        if cand and os.path.exists(cand):
            return cand
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="reduced MC into the quick outputs")
    ap.add_argument("--with-sweep", action="store_true", help="rerun the iterative sweep first")
    ap.add_argument("--skip-tests", action="store_true", help="skip the unit tests")
    ap.add_argument("--no-latex", action="store_true", help="never compile the manuscript")
    ap.add_argument("--tectonic", default=None, help="path of the tectonic executable")
    ap.add_argument("--dry-run", action="store_true", help="print the steps and exit")
    a = ap.parse_args(argv)
    latex = None if a.no_latex else find_tectonic(a.tectonic)
    steps = plan(a.quick, a.with_sweep, a.skip_tests, latex)
    if latex is None and not a.no_latex:
        print("note: tectonic not found (PATH / TECTONIC / --tectonic); LaTeX step skipped")
    env = dict(os.environ)
    if latex:
        env.setdefault("TECTONIC", latex)
    t_all = time.time()
    for i, (label, cmd) in enumerate(steps, 1):
        print("\n=== [%d/%d] %s: %s" % (i, len(steps), label, " ".join(cmd)), flush=True)
        if a.dry_run:
            continue
        t0 = time.time()
        rc = subprocess.call(cmd, cwd=ROOT, env=env)
        print("=== [%d/%d] %s: exit %d in %.0f s" % (i, len(steps), label, rc, time.time() - t0),
              flush=True)
        if rc != 0:
            print("\nFAILED at step '%s'" % label)
            return rc
    print("\nall %d steps OK in %.0f s" % (len(steps), time.time() - t_all))
    return 0


if __name__ == "__main__":
    sys.exit(main())
