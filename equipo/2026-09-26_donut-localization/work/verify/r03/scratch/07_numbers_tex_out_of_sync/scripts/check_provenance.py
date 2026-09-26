# -*- coding: utf-8 -*-
"""Provenance check of the manuscript (called by tests/test_acceptance.py).

1. Flattens ``paper/main.tex``: ``\\input{...}`` / ``\\include{...}`` are resolved recursively
   (relative to ``paper/``, ``.tex`` appended when missing; e.g. ``sections/*.tex`` and
   ``generated/numbers.tex``), LaTeX comments are stripped, and the result is written to a
   temporary file.
2. Runs ``agent-team/bin/check_provenance.py <flat> paper/provenance.json`` as a subprocess
   (every ``\\src{key}`` must resolve to a well-formed, reproducible registry entry).
3. Additionally checks that
   * every ``\\pnum{k}`` / ``\\pnumse{k}`` of the flattened tex is a key of
     ``data/paper_numbers.json`` (``\\pnumse`` also needs an ``se``);
   * every ``number_keys`` listed by an entry of ``paper/provenance.json`` exists there;
   * every ``numbers`` key of ``structure/claims.json`` exists and every ``script`` file exists;
   * ``paper/generated/numbers.tex`` is exactly what ``compute_paper_numbers.numbers_tex``
     produces from ``data/paper_numbers.json`` (in sync).

Prints ``all checks pass`` and exits 0 only if everything passes; otherwise prints every error
and exits 1.

Usage:  python scripts/check_provenance.py [--root DIR] [--checker PATH]
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

DEFAULT_ROOT = os.path.dirname(HERE)
DEFAULT_CHECKER = os.path.join(DEFAULT_ROOT, "agent-team", "bin", "check_provenance.py")

_INPUT_RE = re.compile(r"\\(input|include)\s*\{([^}]+)\}")
_PNUM_RE = re.compile(r"\\pnum\s*\{([^}]+)\}")
_PNUMSE_RE = re.compile(r"\\pnumse\s*\{([^}]+)\}")


def strip_comments(text):
    """Remove LaTeX comments (an unescaped % to the end of the line)."""
    out = []
    for line in text.splitlines():
        m = re.search(r"(?<!\\)%", line)
        out.append(line[:m.start()] if m else line)
    return "\n".join(out)


def flatten(path, base, _stack=None):
    """Recursively inline \\input / \\include of the tex file ``path`` (names relative to
    ``base``); comments stripped.  Raises FileNotFoundError / RecursionError on problems."""
    _stack = [] if _stack is None else _stack
    ap = os.path.abspath(path)
    if ap in _stack:
        raise RecursionError("circular \\input: %s" % " -> ".join(_stack + [ap]))
    with open(ap, encoding="utf-8") as fh:
        text = strip_comments(fh.read())

    def repl(m):
        name = m.group(2).strip()
        cand = os.path.join(base, name)
        if not os.path.exists(cand) and not cand.endswith(".tex"):
            cand += ".tex"
        if not os.path.exists(cand):
            raise FileNotFoundError("\\%s{%s} in %s: file not found" % (m.group(1), name, ap))
        return flatten(cand, base, _stack + [ap])

    return _INPUT_RE.sub(repl, text)


def run(root=DEFAULT_ROOT, checker=DEFAULT_CHECKER, verbose=True):
    """Run all checks; return (ok, list_of_errors, checker_stdout)."""
    errors = []
    paper = os.path.join(root, "paper")
    main_tex = os.path.join(paper, "main.tex")
    prov_path = os.path.join(paper, "provenance.json")
    numbers_path = os.path.join(root, "data", "paper_numbers.json")
    tex_numbers_path = os.path.join(paper, "generated", "numbers.tex")
    claims_path = os.path.join(root, "structure", "claims.json")

    numbers = {}
    try:
        with open(numbers_path, encoding="utf-8") as fh:
            numbers = json.load(fh)
    except (OSError, ValueError) as exc:
        errors.append("cannot read %s: %s" % (numbers_path, exc))

    # 1. flatten
    flat_text = None
    try:
        flat_text = flatten(main_tex, paper)
    except (OSError, RecursionError) as exc:
        errors.append("flatten failed: %s" % exc)

    # 2. agent-team checker
    out = ""
    if flat_text is not None:
        fd, flat = tempfile.mkstemp(suffix="_main_flat.tex", text=True)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(flat_text)
            proc = subprocess.run([sys.executable, checker, flat, prov_path], cwd=root,
                                  capture_output=True, text=True)
            out = proc.stdout + proc.stderr
            if proc.returncode != 0:
                errors.append("agent-team check_provenance failed (exit %d):\n%s"
                              % (proc.returncode, out.strip()))
        finally:
            os.remove(flat)

        # 3a. \pnum keys
        for k in sorted(set(_PNUM_RE.findall(flat_text))):
            if k.strip() not in numbers:
                errors.append("\\pnum{%s} is not a key of data/paper_numbers.json" % k)
        for k in sorted(set(_PNUMSE_RE.findall(flat_text))):
            e = numbers.get(k.strip())
            if not isinstance(e, dict) or "se" not in e:
                errors.append("\\pnumse{%s}: key missing or without 'se' in paper_numbers.json"
                              % k)

    # 3b. provenance number_keys
    prov = {}
    try:
        with open(prov_path, encoding="utf-8") as fh:
            prov = json.load(fh)
    except (OSError, ValueError) as exc:
        errors.append("cannot read %s: %s" % (prov_path, exc))
    if isinstance(prov, dict):
        for key, entry in prov.items():
            nk = entry.get("number_keys", []) if isinstance(entry, dict) else []
            if isinstance(nk, str):
                nk = [nk]
            for k in nk:
                if k not in numbers:
                    errors.append("provenance entry %r: number_keys %r not in paper_numbers.json"
                                  % (key, k))
    else:
        errors.append("paper/provenance.json must be a JSON object")

    # 3c. claims
    try:
        with open(claims_path, encoding="utf-8") as fh:
            claims = json.load(fh)
        if not isinstance(claims, list):
            errors.append("structure/claims.json must be a list")
            claims = []
    except (OSError, ValueError) as exc:
        errors.append("cannot read %s: %s" % (claims_path, exc))
        claims = []
    for c in claims:
        ck = c.get("key", "?") if isinstance(c, dict) else "?"
        if not isinstance(c, dict):
            errors.append("claim is not an object: %r" % (c,))
            continue
        for k in c.get("numbers", []):
            if k not in numbers:
                errors.append("claim %r: number %r not in paper_numbers.json" % (ck, k))
        scr = c.get("script")
        if scr and not os.path.exists(os.path.join(root, scr.split("::", 1)[0])):
            errors.append("claim %r: script %r does not exist" % (ck, scr))

    # 3d. numbers.tex in sync
    if numbers:
        from compute_paper_numbers import numbers_tex
        want = numbers_tex(numbers)
        try:
            with open(tex_numbers_path, encoding="utf-8", newline="") as fh:
                have = fh.read().replace("\r\n", "\n")
            if have != want:
                errors.append("paper/generated/numbers.tex is out of sync with "
                              "data/paper_numbers.json (rerun scripts/compute_paper_numbers.py)")
        except OSError as exc:
            errors.append("cannot read %s: %s" % (tex_numbers_path, exc))

    if verbose:
        if out:
            print(out.strip())
        for e in errors:
            print("ERROR: %s" % e)
        if not errors:
            print("claims: %d, numbers: %d, provenance entries: %d"
                  % (len(claims), len(numbers), len(prov) if isinstance(prov, dict) else 0))
            print("all checks pass")
    return (not errors), errors, out


def main(argv=None):
    ap = argparse.ArgumentParser(description="manuscript provenance check")
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--checker", default=DEFAULT_CHECKER)
    a = ap.parse_args(argv)
    ok, _, _ = run(os.path.abspath(a.root), os.path.abspath(a.checker))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
