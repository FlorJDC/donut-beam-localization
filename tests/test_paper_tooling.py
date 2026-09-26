# -*- coding: utf-8 -*-
"""Tests for the paper tooling (round 4): cache invalidation of _paperstyle.cached (D1),
separate --quick outputs (D2), the numbers.tex fallback (D4) and the shared constants of
_paperconfig (D5)."""
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

import numpy as np  # noqa: E402

import _paperconfig as C  # noqa: E402
import _paperstyle as S  # noqa: E402
import compute_paper_numbers as cpn  # noqa: E402
import make_all_figures as maf  # noqa: E402


class _TmpDirs(unittest.TestCase):
    """Redirect _paperconfig.MC / DATA / FIGDIR to a temporary tree."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="tooling_")
        self._saved = {k: getattr(C, k) for k in ("MC", "DATA", "FIGDIR")}
        C.MC = os.path.join(self.tmp, "data", "mc")
        C.DATA = os.path.join(self.tmp, "data")
        C.FIGDIR = os.path.join(self.tmp, "paper", "figures")
        S.set_quick(False)

    def tearDown(self):
        for k, v in self._saved.items():
            setattr(C, k, v)
        S.set_quick(False)
        shutil.rmtree(self.tmp, ignore_errors=True)


class TestCacheInvalidation(_TmpDirs):

    def _counter(self):
        calls = []

        def fn():
            calls.append(1)
            return {"x": np.arange(3.0), "n": len(calls)}
        return fn, calls

    def test_reuse_when_parameters_unchanged(self):
        fn, calls = self._counter()
        a = S.cached("t", fn, params={"a": 1})
        b = S.cached("t", fn, params={"a": 1})
        self.assertEqual(len(calls), 1)
        np.testing.assert_array_equal(a["x"], b["x"])
        self.assertNotIn("_cache_key", b)

    def test_recompute_when_caller_params_change(self):
        fn, calls = self._counter()
        S.cached("t", fn, params={"a": 1})
        out = S.cached("t", fn, params={"a": 2})
        self.assertEqual(len(calls), 2)
        self.assertEqual(out["n"], 2)

    def test_recompute_when_paperconfig_changes(self):
        """The observed r3 trap: MIS.n_patterns changed but the npz was reused."""
        fn, calls = self._counter()
        S.cached("t", fn)
        old = dict(C.MIS)
        try:
            C.MIS = dict(old, n_patterns=old["n_patterns"] + 1)
            S.cached("t", fn)
            self.assertEqual(len(calls), 2)
            S.cached("t", fn)                          # same new config: reused
            self.assertEqual(len(calls), 2)
        finally:
            C.MIS = old
        S.cached("t", fn)                              # back to the old config: recomputed
        self.assertEqual(len(calls), 3)

    def test_legacy_cache_without_key_is_recomputed(self):
        os.makedirs(C.MC)
        np.savez(os.path.join(C.MC, "t.npz"), x=np.zeros(2))
        fn, calls = self._counter()
        S.cached("t", fn)
        self.assertEqual(len(calls), 1)

    def test_force_and_quick_are_separate(self):
        fn, calls = self._counter()
        S.cached("t", fn)
        S.cached("t", fn, force=True)
        self.assertEqual(len(calls), 2)
        S.cached("t", fn, quick=True)
        self.assertEqual(len(calls), 3)
        self.assertTrue(os.path.exists(os.path.join(C.MC, "t_quick.npz")))

    def test_hash_ignores_paths(self):
        h = S.config_hash()
        C.DATA = os.path.join(self.tmp, "elsewhere")
        self.assertEqual(S.config_hash(), h)
        self.assertNotEqual(S.config_hash(quick=True), h)


class TestQuickOutputs(_TmpDirs):

    def test_savefig_and_summary_go_to_quick_dirs(self):
        import matplotlib.pyplot as plt
        S.set_quick(True)
        fig = plt.figure()
        path = S.savefig(fig, "figX_test")
        S._SUMMARY.clear()
        S.report("k", 1.0)
        spath = S.write_summary("figX")
        self.assertEqual(os.path.dirname(path), os.path.join(C.FIGDIR, "quick"))
        self.assertEqual(os.path.dirname(spath), os.path.join(C.DATA, "quick"))
        self.assertTrue(json.load(open(spath))["quick"])
        self.assertFalse(os.path.exists(os.path.join(C.FIGDIR, "figX_test.pdf")))
        self.assertFalse(os.path.exists(os.path.join(C.DATA, "figX_summary.json")))
        S.set_quick(False)
        fig = plt.figure()
        self.assertEqual(os.path.dirname(S.savefig(fig, "figX_test")), C.FIGDIR)
        self.assertNotIn("quick", json.load(open(S.write_summary("figX"))))
        S._SUMMARY.clear()

    def test_parse_args_sets_quick(self):
        a = S.parse_args("x", ["--quick"])
        self.assertTrue(a.quick)
        self.assertEqual(S.fig_dir(), os.path.join(C.FIGDIR, "quick"))
        S.parse_args("x", [])
        self.assertEqual(S.fig_dir(), C.FIGDIR)

    def test_make_all_figures_quick_pdf_path(self):
        s = os.path.join(SCRIPTS, "fig_8_misalignment.py")
        self.assertEqual(maf.pdf_for(s), os.path.join(maf.FIGDIR, "fig8_misalignment.pdf"))
        self.assertEqual(maf.pdf_for(s, True),
                         os.path.join(maf.FIGDIR, "quick", "fig8_misalignment.pdf"))

    def test_compute_paper_numbers_quick_writes_separate_files(self):
        final_json = os.path.join(self.tmp, "final", "paper_numbers.json")
        final_tex = os.path.join(self.tmp, "final", "numbers.tex")
        q_json = os.path.join(self.tmp, "quick", "paper_numbers.json")
        q_tex = os.path.join(self.tmp, "quick", "numbers.tex")
        saved = (cpn.NUMBERS_PATH, cpn.TEX_PATH, cpn.QUICK_NUMBERS_PATH, cpn.QUICK_TEX_PATH,
                 cpn.compute)
        fake = {"a_key": {"value": 1.5, "unit": "", "script": cpn.SCRIPT, "description": "d"}}
        seen = []
        try:
            cpn.NUMBERS_PATH, cpn.TEX_PATH = final_json, final_tex
            cpn.QUICK_NUMBERS_PATH, cpn.QUICK_TEX_PATH = q_json, q_tex
            cpn.compute = lambda quick=False: (seen.append(quick), dict(fake))[1]
            self.assertEqual(cpn.main(["--quick"]), 0)
        finally:
            (cpn.NUMBERS_PATH, cpn.TEX_PATH, cpn.QUICK_NUMBERS_PATH, cpn.QUICK_TEX_PATH,
             cpn.compute) = saved
        self.assertEqual(seen, [True])
        self.assertFalse(os.path.exists(final_json))
        self.assertFalse(os.path.exists(final_tex))
        d = json.load(open(q_json))
        self.assertIs(d["quick"], True)
        self.assertIn("pnum@a_key", io.open(q_tex, encoding="utf-8").read())
        self.assertNotIn("pnum@quick", io.open(q_tex, encoding="utf-8").read())

    def test_default_quick_paths_are_not_the_final_ones(self):
        self.assertNotEqual(os.path.abspath(cpn.QUICK_NUMBERS_PATH),
                            os.path.abspath(cpn.NUMBERS_PATH))
        self.assertNotEqual(os.path.abspath(cpn.QUICK_TEX_PATH), os.path.abspath(cpn.TEX_PATH))

    def test_fig1_quick_subprocess_leaves_final_outputs_untouched(self):
        """End to end on the real tree: fig_1 --quick must not modify the final PDF/summary."""
        pdf = os.path.join(ROOT, "paper", "figures", "fig1_schematic.pdf")
        summ = os.path.join(ROOT, "data", "fig1_summary.json")
        before = {p: (os.path.getmtime(p), open(p, "rb").read()) for p in (pdf, summ)
                  if os.path.exists(p)}
        proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fig_1_schematic.py"),
                               "--quick"], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        for p, (mt, content) in before.items():
            self.assertEqual(os.path.getmtime(p), mt, p)
            self.assertEqual(open(p, "rb").read(), content, p)
        self.assertTrue(os.path.exists(os.path.join(ROOT, "paper", "figures", "quick",
                                                    "fig1_schematic.pdf")))
        q = json.load(open(os.path.join(ROOT, "data", "quick", "fig1_summary.json")))
        self.assertIs(q["quick"], True)


class TestNumbersTexFallback(unittest.TestCase):

    NUM = {"crb_x_nm": {"value": 1.6051, "unit": "nm", "script": "s", "description": "d",
                        "se": 0.01}}

    def test_fallback_detokenizes_the_key(self):
        t = cpn.numbers_tex(self.NUM)
        self.assertIn(r"\textbf{??\detokenize{#1}??}", t)
        self.assertIn(r"\textbf{??se:\detokenize{#1}??}", t)
        self.assertNotIn(r"\textbf{??#1??}", t)

    def test_metadata_entries_are_skipped(self):
        t = cpn.numbers_tex(dict(self.NUM, quick=True))
        self.assertNotIn("pnum@quick", t)

    @staticmethod
    def _tectonic():
        exe = os.environ.get("TECTONIC") or shutil.which("tectonic")
        return exe if exe and os.path.exists(exe) else None

    def test_unknown_key_compiles_with_tectonic(self):
        exe = self._tectonic()
        if exe is None:
            self.skipTest("tectonic not available (set TECTONIC or put it on PATH)")
        tmp = tempfile.mkdtemp(prefix="pnumtex_")
        try:
            with io.open(os.path.join(tmp, "numbers.tex"), "w", encoding="utf-8") as fh:
                fh.write(cpn.numbers_tex(self.NUM))
            with io.open(os.path.join(tmp, "doc.tex"), "w", encoding="utf-8") as fh:
                fh.write("\\documentclass{article}\n\\input{numbers}\n\\begin{document}\n"
                         "Known \\pnum{crb_x_nm} $\\pm$ \\pnumse{crb_x_nm}; "
                         "unknown \\pnum{no_such_key} and \\pnumse{no_such_key}; "
                         "math $\\pnum{crb_x_nm}$.\n\\end{document}\n")
            proc = subprocess.run([exe, "doc.tex"], cwd=tmp, capture_output=True, text=True,
                                  timeout=600)
            out = proc.stdout + proc.stderr
            self.assertEqual(proc.returncode, 0, out)
            self.assertNotIn("Missing $ inserted", out)
            self.assertTrue(os.path.exists(os.path.join(tmp, "doc.pdf")))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestReproduce(unittest.TestCase):

    def test_step_order(self):
        import reproduce
        labels = [s[0] for s in reproduce.plan()]
        self.assertEqual(labels, ["unit tests", "figures", "paper numbers", "provenance check",
                                  "acceptance test"])
        full = reproduce.plan(with_sweep=True, latex="tectonic")
        self.assertEqual([s[0] for s in full][1], "iterative sweep")
        self.assertEqual(full[-1][0], "LaTeX (tectonic)")
        figs = dict(full)["figures"]
        self.assertIn("--no-cache", figs)
        self.assertNotIn("--quick", figs)

    def test_quick_flags(self):
        import reproduce
        steps = dict(reproduce.plan(quick=True, skip_tests=True))
        self.assertNotIn("unit tests", steps)
        self.assertIn("--quick", steps["figures"])
        self.assertIn("--quick", steps["paper numbers"])
        self.assertNotIn("--quick", steps["provenance check"])

    def test_dry_run_and_shell_wrapper(self):
        import reproduce
        self.assertEqual(reproduce.main(["--dry-run", "--no-latex"]), 0)
        sh = io.open(os.path.join(SCRIPTS, "reproduce.sh"), encoding="utf-8").read()
        self.assertIn("scripts/reproduce.py", sh)


class TestSharedConstants(unittest.TestCase):

    def test_no_duplicated_literals(self):
        """rb=40000, the MLE disk radii and eps=0.002 live only in _paperconfig."""
        for name in ("compute_paper_numbers.py", "fig_5_estimators.py", "fig_7_zero_depth.py"):
            src = io.open(os.path.join(SCRIPTS, name), encoding="utf-8").read()
            self.assertNotIn("40000", src, name)
            self.assertNotIn("search_radius=2.0 * L", src, name)
            self.assertNotIn("np.sqrt(0.002", src, name)
            self.assertNotIn("EPS_C, L_C = 0.002", src, name)
        self.assertEqual(C.N_REP_BIAS0, 40000)
        self.assertEqual(C.MLE_RADIUS_SWEEP_OVER_L, 2.0)
        self.assertEqual(C.MLE_RADIUS_CENTRE_OVER_L, 1.0)
        self.assertEqual(C.EPS_TRANSITION, 0.002)
        self.assertGreaterEqual(C.MIS_POP_N_PATTERNS, 4000)


if __name__ == "__main__":
    unittest.main()
