# -*- coding: utf-8 -*-
"""Tests for scripts/check_provenance.py (manuscript provenance wrapper).

Each test builds a minimal temporary project (paper/main.tex, sections, provenance.json,
data/paper_numbers.json, generated numbers.tex, structure/claims.json) and runs the wrapper with
the real agent-team checker.
"""
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

import check_provenance as cp  # noqa: E402
from compute_paper_numbers import numbers_tex  # noqa: E402

NUMBERS = {
    "crb_x_nm": {"value": 1.6051, "unit": "nm", "script": "scripts/s.py", "description": "d"},
    "eff": {"value": 0.99, "unit": "", "script": "scripts/s.py", "description": "d", "se": 0.005},
}


class _Tree(object):
    def __init__(self, body, provenance=None, claims=None, numbers=None, tex=None):
        self.root = tempfile.mkdtemp(prefix="prov_")
        for d in ("paper/sections", "paper/generated", "data", "structure", "scripts"):
            os.makedirs(os.path.join(self.root, d))
        with open(os.path.join(self.root, "scripts", "s.py"), "w") as fh:
            fh.write("# dummy\n")
        numbers = NUMBERS if numbers is None else numbers
        self.write("data/paper_numbers.json", json.dumps(numbers))
        self.write("paper/generated/numbers.tex", numbers_tex(numbers) if tex is None else tex)
        self.write("paper/provenance.json", json.dumps({} if provenance is None else provenance))
        self.write("structure/claims.json", json.dumps([] if claims is None else claims))
        self.write("paper/main.tex", "\\newcommand{\\src}[1]{}\n\\input{generated/numbers}\n"
                   "% \\src{commented_out_is_ignored}\n\\input{sections/intro}\n")
        self.write("paper/sections/intro.tex", body)

    def write(self, rel, text):
        with open(os.path.join(self.root, rel), "w", encoding="utf-8") as fh:
            fh.write(text)

    def run(self):
        return cp.run(self.root, cp.DEFAULT_CHECKER, verbose=False)

    def close(self):
        shutil.rmtree(self.root, ignore_errors=True)


GOOD_PROV = {"crb_x_nm": {"statement": "centre CRB", "type": "script",
                          "reproduce": "scripts/s.py", "number_keys": ["crb_x_nm"]},
             "eff": {"statement": "efficiency", "type": "script", "reproduce": "scripts/s.py",
                     "detail": "eff", "number_keys": ["eff"]}}
ALL_CLAIMS = [{"key": "c", "numbers": ["eff", "crb_x_nm"], "script": "scripts/s.py"}]


class TestCheckProvenance(unittest.TestCase):

    def _check(self, tree, ok, fragment=None):
        try:
            res, errors, _ = tree.run()
            self.assertEqual(res, ok, errors)
            if fragment is not None:
                self.assertTrue(any(fragment in e for e in errors), errors)
        finally:
            tree.close()

    def test_ok(self):
        self._check(_Tree("The CRB is \\pnum{crb_x_nm}\\src{crb_x_nm} nm "
                          "(\\pnum{eff}\\src{eff} $\\pm$ \\pnumse{eff}).\n", provenance=GOOD_PROV,
                          claims=ALL_CLAIMS), True)

    # ---------------------------------------------------------------- round 4 (D3 holes)
    def test_braceless_input_is_flattened(self):
        t = _Tree("\\input sections/other\n", provenance=GOOD_PROV, claims=ALL_CLAIMS)
        t.write("paper/sections/other.tex", "Hidden \\src{nowhere} and \\pnum{zzz}.\n")
        res, errors, _ = t.run()
        t.close()
        self.assertFalse(res)
        self.assertTrue(any("nowhere" in e for e in errors), errors)
        self.assertTrue(any("zzz" in e for e in errors), errors)

    def test_braceless_input_good_file_passes(self):
        t = _Tree("\\input sections/other \n", provenance=GOOD_PROV, claims=ALL_CLAIMS)
        t.write("paper/sections/other.tex", "Fine \\pnum{eff}\\src{eff}.\n")
        self._check(t, True)

    def test_braceless_input_missing_file_fails(self):
        self._check(_Tree("\\input sections/nothere\n", provenance=GOOD_PROV,
                          claims=ALL_CLAIMS), False, "flatten failed")

    def test_src_with_space_is_checked(self):
        self._check(_Tree("Claim \\src {nowhere}.\n", provenance=GOOD_PROV, claims=ALL_CLAIMS),
                    False, "nowhere")

    def test_src_and_pnum_with_space_resolve(self):
        self._check(_Tree("V \\pnum {eff}\\src {eff}.\n", provenance=GOOD_PROV,
                          claims=ALL_CLAIMS), True)

    def test_pnum_with_space_inside_braces_fails(self):
        self._check(_Tree("V \\pnum{ eff}\\src{eff}.\n", provenance=GOOD_PROV,
                          claims=ALL_CLAIMS), False, "whitespace")

    def test_macro_parameter_is_ignored(self):
        body = ("\\newcommand{\\nn}[1]{\\pnum{#1}}\n"
                "\\newcommand{\\mm}[2]{\\pnumse{##1}}\nV \\pnum{eff}\\src{eff}.\n")
        self._check(_Tree(body, provenance=GOOD_PROV, claims=ALL_CLAIMS), True)

    def test_pnum_without_src_fails(self):
        self._check(_Tree("Value \\pnum{eff} only.\n", provenance=GOOD_PROV, claims=ALL_CLAIMS),
                    False, "no provenance")

    def test_pnumse_without_src_fails(self):
        self._check(_Tree("Value \\pnumse{eff} only.\n", provenance=GOOD_PROV,
                          claims=ALL_CLAIMS), False, "no provenance")

    def test_pnum_backed_by_src_elsewhere_passes(self):
        """The rule is global: a \\src{e} anywhere in the tex whose entry lists k suffices."""
        self._check(_Tree("First \\src{eff}.\n\nLater \\pnum{eff} again \\pnumse{eff}.\n",
                          provenance=GOOD_PROV, claims=ALL_CLAIMS), True)

    def test_src_entry_not_covering_key_fails(self):
        self._check(_Tree("V \\pnum{eff}\\src{crb_x_nm}.\n", provenance=GOOD_PROV,
                          claims=ALL_CLAIMS), False, "no provenance")

    def test_detail_alone_covers_key(self):
        prov = {"e1": {"statement": "s", "type": "script", "reproduce": "scripts/s.py",
                       "detail": "eff"}}
        self._check(_Tree("V \\pnum{eff}\\src{e1}.\n", provenance=prov, claims=ALL_CLAIMS),
                    True)

    def test_comma_separated_src(self):
        self._check(_Tree("V \\pnum{eff} \\pnum{crb_x_nm}\\src{crb_x_nm,eff}.\n",
                          provenance=GOOD_PROV, claims=ALL_CLAIMS), True)

    def test_quick_numbers_fail(self):
        nums = dict(NUMBERS, quick=True)
        self._check(_Tree("x\n", numbers=nums, claims=ALL_CLAIMS), False, "--quick")

    def test_claim_with_empty_numbers_fails(self):
        self._check(_Tree("x\n", claims=ALL_CLAIMS + [{"key": "empty", "numbers": []}]),
                    False, "empty")

    def test_number_without_claim_fails(self):
        self._check(_Tree("x\n", claims=[{"key": "c", "numbers": ["eff"]}]), False,
                    "crb_x_nm")

    def test_normalize(self):
        self.assertEqual(cp.normalize("\\src {a} \\pnum\n{b} \\pnumse  {c} \\srcx {d}"),
                         "\\src{a} \\pnum{b} \\pnumse{c} \\srcx {d}")

    def test_unresolved_src_fails(self):
        self._check(_Tree("Claim \\src{nowhere}.\n", provenance=GOOD_PROV), False,
                    "agent-team check_provenance failed")

    def test_unknown_pnum_fails(self):
        self._check(_Tree("Value \\pnum{not_a_key}.\n"), False, "not_a_key")

    def test_pnumse_needs_se(self):
        self._check(_Tree("Value \\pnumse{crb_x_nm}.\n"), False, "pnumse")

    def test_provenance_number_keys_checked(self):
        prov = {"k": {"statement": "s", "type": "script", "reproduce": "scripts/s.py",
                      "number_keys": ["missing_key"]}}
        self._check(_Tree("Text \\src{k}.\n", provenance=prov), False, "missing_key")

    def test_claims_numbers_and_script_checked(self):
        self._check(_Tree("x\n", claims=[{"key": "c", "numbers": ["zzz"],
                                          "script": "scripts/s.py"}]), False, "zzz")
        self._check(_Tree("x\n", claims=[{"key": "c", "numbers": [],
                                          "script": "scripts/nope.py"}]), False, "nope.py")

    def test_numbers_tex_out_of_sync_fails(self):
        self._check(_Tree("x\n", tex="% stale\n"), False, "out of sync")

    def test_missing_input_fails(self):
        t = _Tree("\\input{sections/missing}\n")
        self._check(t, False, "flatten failed")

    def test_cli_exit_codes_on_real_project(self):
        """The real project must pass (acceptance test uses the same call)."""
        proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_provenance.py")],
                              cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertIn("all checks pass", proc.stdout)
        t = _Tree("Claim \\src{nowhere}.\n")
        try:
            proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_provenance.py"),
                                   "--root", t.root], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(proc.returncode, 1)
            self.assertNotIn("all checks pass", proc.stdout)
        finally:
            t.close()


class TestFlattenAndFormat(unittest.TestCase):

    def test_strip_comments_keeps_escaped_percent(self):
        self.assertEqual(cp.strip_comments("50\\% of x % comment\n% all"), "50\\% of x \n")

    def test_numbers_tex_format(self):
        from compute_paper_numbers import _fmt
        self.assertEqual(_fmt(3), "3")
        self.assertEqual(_fmt(1.60509557), "1.605")
        self.assertEqual(_fmt(12345.6), "12346")
        self.assertIn("\\times10^{-5}", _fmt(7.04e-5))
        self.assertEqual(_fmt(0.0), "0")
        t = numbers_tex(NUMBERS)
        self.assertIn("\\csname pnum@crb_x_nm\\endcsname{1.605}", t)
        self.assertIn("\\csname pnumse@eff\\endcsname{0.005}", t)
        self.assertNotIn("pnumse@crb_x_nm", t)


if __name__ == "__main__":
    unittest.main()
