#!/usr/bin/env python3
"""recouvrement.py — U2 (equilibrium) contre U3 lag 1000 ns : les deux ensembles occupent-ils la meme region ?
Pour chaque systeme : RMSD Calpha (superposition optimale sur TOUS les Calpha du complexe, appariement ORDINAL des residus car les deux ensembles partagent la topologie) entre
  - toutes les paires internes de U2 (25*24/2), toutes les paires internes de U3-1000 (idem), et les 25*25 paires CROISEES U2 x U3-1000 ;
  - rapport croise / moyenne des deux internes (1 = meme distribution vue par ses distances ; > 1 = ensembles separes).
Ce n est PAS une metrique DynaPIN : c est un descripteur de recouvrement, declare comme tel. Ecrit generation/comparaison/recouvrement.json."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L
import mdtraj as md, numpy as np

OUT = L.G / "comparaison"
OUT.mkdir(exist_ok=True)


def charger(unite):
    d = L.G / "sorties" / unite
    return md.load_xtc(str(d / "samples.xtc"), top=str(d / "topology.cif"))


res = {}
for pdb in sorted(L.corpus_28()):
    a, b = charger(f"{pdb}__U2__eq"), charger(f"{pdb}__U3__eq_lag1000ns")
    ca = a.topology.select("name CA")
    cb = b.topology.select("name CA")
    if len(ca) != len(cb):
        res[pdb] = {"statut": "ecarte", "raison": f"{len(ca)} CA contre {len(cb)}"}
        continue
    tous = md.join([a.atom_slice(ca), b.atom_slice(ca)])
    n = a.n_frames
    M = np.zeros((2 * n, 2 * n))
    for i in range(2 * n):
        M[i] = md.rmsd(tous, tous, frame=i) * 10.0           # A
    iu = np.triu_indices(n, k=1)
    int_a, int_b = M[:n, :n][iu], M[n:, n:][iu]
    cross = M[:n, n:].ravel()
    res[pdb] = {"statut": "ok", "n": int(n), "rmsd_interne_U2_moy_A": round(float(int_a.mean()), 3), "rmsd_interne_U3_1000_moy_A": round(float(int_b.mean()), 3),
                "rmsd_croise_moy_A": round(float(cross.mean()), 3),
                "rapport_croise_sur_interne_moyen": round(float(cross.mean() / np.mean([int_a.mean(), int_b.mean()])), 3),
                "rmsd_interne_U2_mediane_A": round(float(np.median(int_a)), 3), "rmsd_interne_U3_1000_mediane_A": round(float(np.median(int_b)), 3), "rmsd_croise_mediane_A": round(float(np.median(cross)), 3)}
    print(f"RECOUVREMENT {pdb} interne U2 {res[pdb]['rmsd_interne_U2_moy_A']} · interne U3-1000 {res[pdb]['rmsd_interne_U3_1000_moy_A']} · croise {res[pdb]['rmsd_croise_moy_A']} · rapport {res[pdb]['rapport_croise_sur_interne_moyen']}", flush=True)
(OUT / "recouvrement.json").write_text(json.dumps(res, indent=1) + "\n")
ko = [p for p, r in res.items() if r["statut"] != "ok"]
print(f"RECOUVREMENT_FIN {len(res)} systemes, ecartes {ko} (rejets = {len(ko)})")
sys.exit(0 if not ko else 1)
