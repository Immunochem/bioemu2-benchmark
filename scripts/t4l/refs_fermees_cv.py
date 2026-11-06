import glob, json, os, sys, numpy as np, mdtraj as md
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C
T = "/gpfs/phherman/bioemu2_bench/t4l"; Z = f"{T}/donnees/zenodo/extrait/T4L-bechmarking/Closed-state"
out = {}
for lab in ("Exposed-closed", "Buried-closed"):
    a = [C.cv(md.load(f)) for f in sorted(glob.glob(f"{Z}/{lab}/*.pdb"))]
    for k in ("d1", "d2", "p"):
        v = np.array([r[k][0] for r in a]); out[f"{lab}_{k}_min_q25_med_q75_max"] = [round(float(x), 3) for x in np.percentile(v, [0, 25, 50, 75, 100])]
c = C.cv(md.load(f"{T}/x0/x0_5lzm.pdb")); out["cristal"] = {k: round(float(c[k][0]), 3) for k in ("d1", "d2", "p")}
print(json.dumps(out, indent=1)); json.dump(out, open(f"{T}/resultats/refs_fermees_cv.json", "w"), indent=1)
