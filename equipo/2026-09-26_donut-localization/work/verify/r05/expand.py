# -*- coding: utf-8 -*-
"""Print manuscript sections with \\pnum{k} replaced by <k=value> (verifier r05 reading aid)."""
import re, sys, io, os
root = r"C:\Users\BANGHO\Documents\GithubPRO\donut-beam-localization\paper"
nums = {}; ses = {}
PAT = re.compile(r"\\expandafter\\def\\csname pnum(se)?@(\S+)\\endcsname\{(.*)\}\s*$")
for line in io.open(os.path.join(root, "generated", "numbers.tex"), encoding="utf-8"):
    m = PAT.match(line)
    if m:
        (ses if m.group(1) else nums)[m.group(2)] = m.group(3)
def sub(t):
    t = re.sub(r"\\pnumse\{([^}]*)\}", lambda m: "<SE:%s>" % ses.get(m.group(1), "??"), t)
    t = re.sub(r"\\pnum\{([^}]*)\}", lambda m: "<%s=%s>" % (m.group(1), nums.get(m.group(1), "??")), t)
    t = re.sub(r"\\src\{[^}]*\}", "", t)
    return t
for f in sys.argv[1:]:
    p = os.path.join(root, "sections", f + ".tex")
    print("=" * 20, f)
    for i, l in enumerate(io.open(p, encoding="utf-8"), 1):
        sys.stdout.write("%d: %s" % (i, sub(l)))
