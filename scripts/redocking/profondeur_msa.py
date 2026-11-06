import json, glob, os, sys
R = "/gpfs/phherman/bioemu2_bench/redocking"
niv = {}
h = json.load(open(R + "/strate_B/homologie.json"))["systemes"]
for pid, v in h.items():
    niv[pid] = 3 if v["id90"]["entrees_avec_les_deux"] else 2 if v["id50"]["entrees_avec_les_deux"] else 1 if v["id30"]["entrees_avec_les_deux"] else 0
rows = []
for f in sorted(glob.glob(R + "/msa/*.json")):
    d = json.load(open(f)); i = os.path.basename(f)[:-5]
    if not isinstance(d, dict) or "chaines" not in d: continue
    ch = d["chaines"]; 
    rows.append((i, [ch[c]["non_apparie"] for c in sorted(ch)], [ch[c].get("apparie") for c in sorted(ch)]))
sel = sys.argv[1:] 
for i, u, p in rows:
    if sel and not any(i.startswith(s) or i in sel for s in sel): continue
    print(i, "niveau" , niv.get(i, "-"), "non apparie", u, "apparie", p)
