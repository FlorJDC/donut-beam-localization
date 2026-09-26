import sys, os
sys.path.insert(0, "tests"); sys.path.insert(0, "scripts")
from test_check_provenance import _Tree, GOOD_PROV
import check_provenance as cp
def go(name, body, **kw):
    t = _Tree(body, **kw)
    try:
        ok, errs, out = t.run()
        print("%-40s ok=%s %s" % (name, ok, [e[:110] for e in errs]))
    finally:
        t.close()
go("input without braces hides bad src", "\input sections/other\n")
t = None
# need other file: build manually
t = _Tree("\input sections/other\n")
t.write("paper/sections/other.tex", "Claim \src{nowhere} \pnum{zzz}\n")
print("input-no-braces:", t.run()[:2]); t.close()
go("escaped linebreak then comment", "a \\% \src{nowhere}\n", provenance=GOOD_PROV)
go("pnum with leading space", "v \pnum{ crb_x_nm}\n")
go("pnum without src", "v \pnum{crb_x_nm} nm\n")
prov = {"crb_x_nm": {"statement": "s", "type": "script", "reproduce": "sections/intro.tex"}}
go("paper-relative reproduce path", "v \src{crb_x_nm}\n", provenance=prov)
go("input with .tex ext", "\input{generated/numbers.tex}\n")
go("pnum inside verbatim-like \iffalse", "\iffalse \pnum{zzz}\fi\n")
