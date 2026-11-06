#!/usr/bin/env python3
"""resume_plancher.py — resume du plancher du hasard par strate (relit resultats/plancher/*.json). Imprime ; ecrit resultats/plancher_par_strate.json. Les placements vides (F1 impossible : aucun placement sans choc a la distance native) sont COMPTES."""
import glob, json
import numpy as np
R = "/gpfs/phherman/bioemu2_bench/redocking"
P = [json.load(open(f)) for f in sorted(glob.glob(R + "/resultats/plancher/*.json"))]
out = {}
for s in ("A_corpus", "A_apparie", "B"):
    for k in ("F1", "F2"):
        sel = [p for p in P if p["strate"] == s]
        if not sel:
            continue
        vides = [p["id"] for p in sel if not p["plancher"][k]]
        ok = [p for p in sel if p["plancher"][k]]
        fn = np.concatenate([p["plancher"][k]["fnat"] for p in ok]) if ok else np.array([])
        dq = np.concatenate([p["plancher"][k]["dockq"] for p in ok]) if ok else np.array([])
        cl = np.concatenate([p["plancher"][k]["classe_capri"] for p in ok]) if ok else np.array([])
        p99 = [float(np.percentile(p["plancher"][k]["fnat"], 99)) for p in ok]
        mx = [float(max(p["plancher"][k]["fnat"])) for p in ok]
        o = {"systemes": len(sel), "placements_vides": len(vides), "ids_vides": vides, "n_placements": int(len(fn)), "fnat_p50_global": float(np.median(fn)) if len(fn) else None,
             "fnat_p99_global": float(np.percentile(fn, 99)) if len(fn) else None, "fnat_max_global": float(fn.max()) if len(fn) else None, "p99_par_systeme_max": max(p99) if p99 else None,
             "systemes_dont_p99_ge_0.1": int(sum(x >= 0.1 for x in p99)), "systemes_dont_max_ge_0.3": int(sum(x >= 0.3 for x in mx)), "dockq_max": float(dq.max()) if len(dq) else None,
             "n_ge_acceptable": int((cl != "incorrect").sum()) if len(cl) else 0, "n_ge_medium": int(np.isin(cl, ["medium", "high"]).sum()) if len(cl) else 0}
        out[f"{s}|{k}"] = o
        print(s, k, json.dumps(o, ensure_ascii=False))
json.dump(out, open(R + "/resultats/plancher_par_strate.json", "w"), indent=1)
