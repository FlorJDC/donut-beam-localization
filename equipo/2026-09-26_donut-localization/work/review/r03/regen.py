# Re-run compute_paper_numbers sections in memory (no file writes) and diff against the committed JSON.
import sys, json, time
sys.path.insert(0, "scripts")
import compute_paper_numbers as cpn
R = cpn.Registry()
t = time.time()
for name, fn in (("crb", lambda: cpn.crb_section(R)), ("camera", lambda: cpn.camera_section(R)),
                 ("eps", lambda: cpn.eps_section(R)), ("iterative", lambda: cpn.iterative_section(R, False)),
                 ("estimators", lambda: cpn.estimator_section(R, False)),
                 ("adaptive", lambda: cpn.adaptive_coverage_section(R, False)),
                 ("vectorial", lambda: cpn.vectorial_section(R))):
    fn(); print("# %s %.1f s" % (name, time.time() - t), flush=True)
pn = json.load(open("data/paper_numbers.json", encoding="utf-8"))
bad = 0
for k, e in R.items():
    v = pn.get(k)
    if v is None: print("MISSING in committed", k); bad += 1; continue
    if v["value"] != e["value"] or v.get("se") != e.get("se"):
        print("DIFF", k, v["value"], e["value"], v.get("se"), e.get("se")); bad += 1
print("regenerated", len(R), "keys; diffs", bad)
