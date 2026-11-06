import json, glob, collections
R = "/gpfs/phherman/bioemu2_bench/redocking/sorties"
par = collections.defaultdict(dict)
for f in glob.glob(R + "/*/*/meta.json"):
    m = json.load(open(f))
    if m.get("statut") != "ok": continue
    for lag, v in m["lags"].items():
        par[m["id"]].setdefault(lag, []).append((v["s_par_structure"], m["tokens"], v["batch_size"], v["jobs_concurrents_gpu"]["max"]))
print("systeme tokens | s/structure lag 50 | s/structure lag 1000 | ratio 1000/50 (images, batch, concurrence max)")
rat = []
for sid, d in sorted(par.items(), key=lambda kv: list(kv[1].values())[0][0][1]):
    if "50" in d and "1000" in d:
        a, b = d["50"][0], d["1000"][0]
        rat.append(b[0] / a[0])
        print(f"{sid:6s} {a[1]:4d} | {a[0]:6.2f} (conc {a[3]}) | {b[0]:6.2f} (conc {b[3]}) | {b[0] / a[0]:.2f}")
if rat:
    import statistics
    print(f"ratio 1000/50 : mediane {statistics.median(rat):.2f}, min {min(rat):.2f}, max {max(rat):.2f} sur {len(rat)} systemes (conditions differentes : 24 images a lag 50, 15 images a lag 1000, concurrence GPU 3-4 aux deux)")
