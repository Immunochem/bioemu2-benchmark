#!/usr/bin/env python3
"""verifier_enligne.py — les CV calculees EN LIGNE (cv.npz, sur les coordonnees en memoire) == les CV recalculees HORS LIGNE depuis brut/ (xtc + cif, relus) ?
Borne ecrite avant lecture : le .xtc quantifie les coordonnees a 0,001 nm (precision par defaut) -> erreur de coordonnee <= 0,0005 nm par axe ->
erreur de distance <= 2*sqrt(3)*0,0005 = 0,0018 nm. Le critere « 1e-4 nm » de la decision du 2026-11-06 ne tenait pas compte de cette quantification :
corrige a 2e-3 nm AVANT de lire le resultat. Pour p (rapport de deux produits scalaires), la borne est empirique : on la rapporte, on ne la decrete pas."""
import json, os, sys
import numpy as np, mdtraj as md
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C
T = "/gpfs/phherman/bioemu2_bench/t4l"
racine = sys.argv[1]
bad, out = 0, {}
for u in sorted(os.listdir(f"{T}/{racine}")):
    d = f"{T}/{racine}/{u}"
    if not os.path.isfile(f"{d}/cv.npz"):
        continue
    z = np.load(f"{d}/cv.npz")
    try:
        r = C.cv(md.load_xtc(f"{d}/brut/samples.xtc", top=f"{d}/brut/topology.cif"))
    except Exception as e:
        out[u] = {"ERREUR": repr(e)[:300]}; bad += 1; continue
    e = {"n": int(len(z["d1"])), "max_d1_nm": float(np.abs(r["d1"] - z["d1"]).max()), "max_d2_nm": float(np.abs(r["d2"] - z["d2"]).max()),
         "max_p": float(np.abs(r["p"] - z["p"]).max()), "etats_differents_en_ligne_hors_ligne": int((C.etat(r["d1"], r["p"]) != C.etat(z["d1"], z["p"])).sum())}
    e["ok_borne_2e-3_nm"] = bool(e["max_d1_nm"] < 2e-3 and e["max_d2_nm"] < 2e-3)
    bad += (not e["ok_borne_2e-3_nm"]); out[u] = e
    print(u, json.dumps(e))
print(f"UNITES EN ECART : {bad} sur {len(out)}")
json.dump(out, open(f"{T}/resultats/verif_enligne_{racine}.json", "w"), indent=1)
sys.exit(1 if bad else 0)
