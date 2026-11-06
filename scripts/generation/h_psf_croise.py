#!/usr/bin/env python3
"""h_psf_croise.py — l etat de protonation des histidines est LU dans le PSF de la MD (h_core.variantes_depuis_psf, le meme code que l outil
de reconstruction) et CROISE avec prep_fix/curation/protonation_his.tsv (releve independant de la curation). Rien n est code en dur :
un desaccord de nombre ou d etat par systeme est un REJET. Ecrit hydrogenes/psf_croise.json. A lancer avec le python de l env evomd (openmm)."""
import csv, json, os, sys
from collections import Counter
from pathlib import Path
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/hydrogenes/scripts")
import h_core as H

G = Path("/gpfs/phherman/bioemu2_bench/generation")
B = Path("/gpfs/phherman/bioemu2_bench")
STD = {"HSD": "HID", "HSE": "HIE", "HSP": "HIP"}
tsv = {}
for r in csv.DictReader(open(B / "prep_fix/curation/protonation_his.tsv"), delimiter="\t"):
    if r["replique"] == "rep1":
        tsv.setdefault(r["pdb"], []).append(STD.get(r["etat_deduit_des_H"], r["etat_deduit_des_H"]))
ponts = json.loads((G / "ponts/reference.json").read_text())["systemes"]
HR = os.environ.get("H_RACINE", "hydrogenes")
taches = [l.split("\t") for l in (G / HR / "taches.tsv").read_text().splitlines() if l.strip()]
out, rej = {}, []
for pdb, liste, psf in taches:
    ref = [l for l in open(liste) if l.strip()][0].strip()
    V = H.variantes_depuis_psf(psf, ref)
    his = Counter(v for v in V.values() if v.startswith("HI"))
    t = Counter(tsv.get(pdb, []))
    cys = Counter(v for v in V.values() if v.startswith("CY"))
    ok_cys = cys.get("CYX", 0) == 2 * ponts[pdb]["n_ponts"] and cys.get("CYS", 0) == ponts[pdb]["n_cys_libres"]
    ok = dict(his) == dict(t) and ok_cys
    out[pdb] = {"his_psf": dict(his), "his_tsv_curation": dict(t), "cys_psf": dict(Counter(v for v in V.values() if v.startswith("CY"))), "accord_his": dict(his) == dict(t), "accord_cys_avec_ponts_reference": ok_cys, "n_ponts_reference": ponts[pdb]["n_ponts"], "n_cys_libres_reference": ponts[pdb]["n_cys_libres"], "accord": ok}
    if not ok:
        rej.append(pdb)
    print(f"PSF_CROISE {pdb} psf={dict(his)} tsv={dict(t)} cys={out[pdb]['cys_psf']} {'OK' if ok else 'DESACCORD'}", flush=True)
(G / HR / "psf_croise.json").write_text(json.dumps({"systemes": out, "desaccords": rej}, indent=1) + "\n")
print(f"PSF_CROISE_BILAN {len(out)} systemes, desaccords = {len(rej)} {rej}")
sys.exit(0 if not rej else 1)
