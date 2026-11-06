#!/usr/bin/env python3
"""controle_mesure.py — CONTROLE POSITIF DE LA MESURE. Mes fnat / iRMSD / lRMSD / DockQ (redock_lib.Reference.mesurer) contre ceux de DynaPIN, STRUCTURE PAR STRUCTURE, sur les piles du pilote 6CWG :
  reference = 6CWG/x0/x0_md0.pdb (image 0 de la MD), echantillons = 6CWG/prop/<unite>/brut/samples.xtc, table DynaPIN = hydrogenes/out/R6_all/piles/cristal_vsMD0__<unite>/dynapin/tables/qc_trajectory_metrics.csv
  (image 0 de la table = la reference elle-meme : DockQ 1, fnat 1 ; les echantillons sont les lignes suivantes, dans l ordre).
Si l ecart depasse les seuils, le script sort en 1 et DIT lequel : une mesure qui ne retrouve pas l outil de reference ne decide rien."""
import csv, json, os, sys
from pathlib import Path
os.environ["CUDA_VISIBLE_DEVICES"] = ""
sys.dont_write_bytecode = True
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
sys.path.insert(0, str(R / "scripts"))
import numpy as np, mdtraj as md
import redock_lib as RL
import redock_bioemu as RB

W = Path("/gpfs/phherman/bioemu2_bench/6CWG")
H = Path("/gpfs/phherman/bioemu2_bench/hydrogenes/out/R6_all/piles")
info = RB.systeme_depuis_x0(W / "x0/x0_md0.pdb")
Rf = RL.Reference(info["ref"])
print(f"reference 6CWG x0_md0 : {info['tokens']} residus, d_COM {Rf.d_nat:.2f} A, {Rf.n_contacts_nat} couples de residus en contact, recepteur = chaine {Rf.rec}", flush=True)
seuils = {"fnat": 0.02, "irmsd_A": 0.05, "lrmsd_A": 0.10, "dockq": 0.02}
cols = {"fnat": "fnat", "irmsd_A": "iRMSD", "lrmsd_A": "lRMSD", "dockq": "DockQ"}
pire = {k: 0.0 for k in seuils}
bilan = []
for pile, unite in (("cristal_vsMD0__cristal_lag0001", "cristal_lag0001"), ("cristal_vsMD0__cristal_lag0050", "cristal_lag0050"), ("cristal_vsMD0__cristal_lag0200", "cristal_lag0200"),
                    ("cristal_vsMD0__cristal_lag1000", "cristal_lag1000"), ("md0_self__md0_lag0050", "md0_lag0050")):
    d = W / "prop" / unite / "brut"
    tr = md.load_xtc(str(d / "samples.xtc"), top=str(d / "topology.cif"))
    m = Rf.mesurer(tr)
    rows = list(csv.DictReader(open(H / f"{pile}/dynapin/tables/qc_trajectory_metrics.csv")))
    rows.sort(key=lambda r: int(float(r["Frame"])))
    z = rows[0]
    garde0 = abs(float(z["DockQ"]) - 1) < 1e-6 and abs(float(z["fnat"]) - 1) < 1e-6
    rows = rows[1:]
    if len(rows) != tr.n_frames:
        print(f"{unite} : {len(rows)} lignes DynaPIN pour {tr.n_frames} images -> ARRET (appariement ordinal impossible)"); sys.exit(2)
    ligne = {"unite": unite, "n": tr.n_frames, "image0_dynapin_ok": garde0}
    for k, c in cols.items():
        dd = np.array(m[k]) - np.array([float(r[c]) for r in rows])
        ligne[k] = {"mediane_moi": round(float(np.median(m[k])), 4), "mediane_dynapin": round(float(np.median([float(r[c]) for r in rows])), 4), "ecart_max_abs": round(float(np.abs(dd).max()), 4),
                    "ecart_median": round(float(np.median(dd)), 4)}
        pire[k] = max(pire[k], float(np.abs(dd).max()))
    bilan.append(ligne)
    print(unite, json.dumps({k: ligne[k] for k in cols}), "image0 DynaPIN =", "ok" if garde0 else "NON", flush=True)
ok = all(pire[k] < seuils[k] for k in seuils)
(R / "resultats").mkdir(exist_ok=True)
(R / "resultats/controle_mesure.json").write_text(json.dumps({"seuils": seuils, "pire_ecart_abs": pire, "ok": ok, "unites": bilan}, indent=1) + "\n")
print("PIRE ECART |moi - DynaPIN| :", {k: round(v, 4) for k, v in pire.items()}, "seuils", seuils, "->", "OK" if ok else "HORS SEUIL")
sys.exit(0 if ok else 1)
