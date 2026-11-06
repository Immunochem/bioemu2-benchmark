#!/usr/bin/env python3
"""rmsd_cristal.py — DISPERSION structurale : RMSD Ca (162 residus, appariement ordinal, superposition optimale) de chaque structure au cristal 5LZM.
Ensembles : BioEmu v1 (499), les sept GenAI de l'auteur, et chaque unite BioEmu2 de `unites`. Rend mediane et quantiles 5/95/max, par ensemble."""
import json, os, sys, numpy as np, mdtraj as md
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C
T = "/gpfs/phherman/bioemu2_bench/t4l"; Z = f"{T}/donnees/zenodo/extrait/T4L-bechmarking/Genai-models"
N = 162
def ca(tr):
    res = list(list(tr.topology.chains)[0].residues)
    C.indices(tr.topology)                                        # garde ordinal (noms de residus)
    return tr.atom_slice([a.index for r in res[:N] for a in r.atoms if a.name == "CA"])
ref = ca(md.load(f"{T}/x0/x0_5lzm.pdb"))
jeux = {"BioEmu_v1": (f"{Z}/BioEmu/samples.xtc", f"{Z}/BioEmu/topology.pdb")}
for nom, (d, x, t) in {"AlphaFlow": ("AlphaFlow", "combined.xtc", "prot-masses.pdb"), "ESMFlow": ("ESMFlow", "combined.xtc", "prot-masses.pdb"), "AF-cluster": ("AF-cluster", "combined.xtc", "prot-masses.pdb"),
                       "ConforFold": ("ConforFold", "combined.xtc", "prot-masses.pdb"), "ConfRover": ("ConfRover", "T4L.xtc", "T4L.pdb"), "rMSA-AF2": ("rMSA-AF2", "combined.xtc", "prot-masses.pdb"),
                       "rMSA-AF2-recyc": ("rMSA-AF2-recyc", "combined.xtc", "prot-masses.pdb")}.items():
    jeux[nom] = (f"{Z}/{d}/{x}", f"{Z}/{d}/{t}")
for rac in sys.argv[1:]:
    for u in sorted(os.listdir(f"{T}/{rac}")):
        if os.path.isfile(f"{T}/{rac}/{u}/brut/samples.xtc"):
            jeux[f"BioEmu2:{u}"] = (f"{T}/{rac}/{u}/brut/samples.xtc", f"{T}/{rac}/{u}/brut/topology.cif")
out, rej = {}, []
for k, (x, t) in jeux.items():
    try:
        tr = ca(md.load_xtc(x, top=t)) if x.endswith(".xtc") and not t.endswith(".cif") else ca(md.load_xtc(x, top=t))
        r = md.rmsd(tr, ref, frame=0) * 10.0
        out[k] = {"n": int(tr.n_frames), "rmsd_cristal_A_min_5_50_95_max": [round(float(v), 2) for v in np.percentile(r, [0, 5, 50, 95, 100])]}
        print(k, out[k], flush=True)
    except Exception as e:
        rej.append({"jeu": k, "raison": repr(e)[:200]}); print("REJET", k, repr(e)[:200], flush=True)
out["rejets"] = rej
nom = "rmsd_cristal.json" if sys.argv[1:] == ["unites"] else "rmsd_cristal_" + "_".join(sys.argv[1:]) + ".json"   # ne jamais ecraser le fichier de la production
json.dump(out, open(f"{T}/resultats/{nom}", "w"), indent=1)
print("REJETS :", len(rej))
