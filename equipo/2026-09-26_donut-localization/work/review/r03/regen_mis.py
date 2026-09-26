import sys, json, time
sys.path.insert(0, "scripts")
import compute_paper_numbers as cpn
R = cpn.Registry(); t = time.time()
cpn.misalignment_section(R, False)
print("# misalignment %.1f s" % (time.time() - t))
pn = json.load(open("data/paper_numbers.json", encoding="utf-8"))
bad = 0
for k, e in R.items():
    v = pn.get(k)
    if v is None or v["value"] != e["value"] or v.get("se") != e.get("se"):
        print("DIFF", k, v and v["value"], e["value"]); bad += 1
print("regenerated", len(R), "keys; diffs", bad, "; committed total", len(pn))
