#!/usr/bin/env python3
"""proxy_cb.py — BioEmu v1 n'a pas le cycle de Phe4 (topology squelette + CB, 809 atomes) : p n'y est PAS defini.
Un PROXY (le CB de Phe4 a la place du centre du cycle) est calibre sur les structures completes de l'auteur, puis applique a BioEmu v1.
Calibration : 2197 structures GenAI + 221 structures fermees (MD de l'auteur) + 5LZM. Critere ecrit avant : accord de signe >= 90 %.
Le proxy n'est PAS sa definition : tout ce qui en sort est etiquete « proxy »."""
import json, os, sys
import numpy as np, mdtraj as md
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C
T = "/gpfs/phherman/bioemu2_bench/t4l"
Z = f"{T}/donnees/zenodo/extrait/T4L-bechmarking"
M = {"AlphaFlow": ("AlphaFlow", "combined.xtc", "prot-masses.pdb"), "ESMFlow": ("ESMFlow", "combined.xtc", "prot-masses.pdb"), "AF-cluster": ("AF-cluster", "combined.xtc", "prot-masses.pdb"),
     "ConforFold": ("ConforFold", "combined.xtc", "prot-masses.pdb"), "ConfRover": ("ConfRover", "T4L.xtc", "T4L.pdb"), "rMSA-AF2": ("rMSA-AF2", "combined.xtc", "prot-masses.pdb"),
     "rMSA-AF2-recyc": ("rMSA-AF2-recyc", "combined.xtc", "prot-masses.pdb")}
rej, pr, pc = [], [], []
par = {}
for nom, (d, x, t) in M.items():
    r = C.cv(md.load(f"{Z}/Genai-models/{d}/{x}", top=f"{Z}/Genai-models/{d}/{t}"))
    pr.append(r["p"]); pc.append(r["p_cb_proxy"]); par[nom] = (r["p"], r["p_cb_proxy"])
import glob
for lab in ("Exposed-closed", "Buried-closed"):
    a, b = [], []
    for f in sorted(glob.glob(f"{Z}/Closed-state/{lab}/*.pdb")):
        try:
            r = C.cv(md.load(f)); a.append(r["p"][0]); b.append(r["p_cb_proxy"][0])
        except Exception as e:
            rej.append({"fichier": f, "raison": repr(e)[:200]})
    pr.append(np.array(a)); pc.append(np.array(b)); par[lab] = (np.array(a), np.array(b))
P, Q = np.concatenate(pr), np.concatenate(pc)
ok = np.isfinite(P) & np.isfinite(Q)
acc = float(((P > 0) == (Q > 0))[ok].mean())
res = {"n_calibration": int(ok.sum()), "accord_de_signe": acc, "pearson": float(np.corrcoef(P[ok], Q[ok])[0, 1]),
       "pente_moindres_carres": float(np.polyfit(P[ok], Q[ok], 1)[0]), "critere_ecrit_avant_accord_ge_0_90": bool(acc >= 0.90),
       "par_jeu": {k: {"n": int(len(a)), "accord_de_signe": float(((a > 0) == (b > 0)).mean())} for k, (a, b) in par.items()},
       "confusion_ring_vs_proxy": {"ring_pos_proxy_pos": int(((P > 0) & (Q > 0))[ok].sum()), "ring_pos_proxy_neg": int(((P > 0) & (Q < 0))[ok].sum()),
                                   "ring_neg_proxy_pos": int(((P < 0) & (Q > 0))[ok].sum()), "ring_neg_proxy_neg": int(((P < 0) & (Q < 0))[ok].sum())}}
# BioEmu v1
tr = md.load(f"{Z}/Genai-models/BioEmu/samples.xtc", top=f"{Z}/Genai-models/BioEmu/topology.pdb")
r = C.cv(tr)
v1 = r["p_cb_proxy"]
res["bioemu_v1_proxy"] = {"n": int(v1.size), "n_cb_nan": int((~np.isfinite(v1)).sum()), "p_gt_0": int((v1 > 0).sum()), "p_lt_0": int((v1 < 0).sum()),
                          "p_min_med_max": [float(np.nanmin(v1)), float(np.nanmedian(v1)), float(np.nanmax(v1))], "d1_ferme": int((r["d1"] < 2.5).sum()),
                          "exposed_open_pct": float(100 * (v1 > 0).mean()), "buried_open_pct": float(100 * (v1 < 0).mean())}
res["rejets"] = rej; res["n_rejets"] = len(rej)
json.dump(res, open(f"{T}/resultats/proxy_cb.json", "w"), indent=1)
print(json.dumps(res, indent=1))
