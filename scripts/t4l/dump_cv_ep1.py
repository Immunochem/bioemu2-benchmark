"""dump_cv.py : exporte (d1, p) des 7 methodes de Bhakat (ses distances.npy et locked.pkl, lus par reproduire_v1.lire_pkl)
et de BioEmu2 (cv.npz, unite eq164_n500). Rien n'est recalcule : lecture seule, le trace ne depend que de fichiers sources."""
import json, sys, hashlib
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/t4l/scripts")
import reproduire_v1 as R
out, rej = {"methodes": {}, "rejets": []}, []
for nom, (dz, xtc, pdb, dan) in R.METHODES.items():
    if nom == "BioEmu":
        continue            # p non defini (809 atomes, pas de cycle de Phe4)
    try:
        d = np.load(f"{dan}/distances.npy")
        p = R.lire_pkl(f"{dan}/locked.pkl")
        assert d.shape[0] == p.shape[0], (d.shape, p.shape)
        out["methodes"][nom] = {"d1": [round(float(x), 5) for x in d[:, 0]], "p": [round(float(x), 5) for x in p],
                                "sha_distances": hashlib.sha256(open(f"{dan}/distances.npy","rb").read()).hexdigest()[:16]}
    except Exception as e:
        out["rejets"].append({"methode": nom, "erreur": repr(e)})
z = np.load("/gpfs/phherman/bioemu2_bench/t4l/unites/eq164_n500/cv.npz")
out["bioemu2"] = {"d1": [round(float(x), 5) for x in z["d1"]], "p": [round(float(x), 5) for x in z["p"]],
                  "physique": [bool(x) for x in z["physique"]]}
out["n_rejets"] = len(out["rejets"])
json.dump(out, sys.stdout)
