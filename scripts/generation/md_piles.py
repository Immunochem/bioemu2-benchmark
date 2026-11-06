#!/usr/bin/env python3
"""md_piles.py — piles de REFERENCE MD pour DynaPIN : 28 systemes x 3 repliques, stride 5 (200 images de 1000), atomes avec H du champ de force.
DERIVE de 6CWG/scripts/preparer.py:72-99 (memes choix : xtc + structure.pdb de la MD, stride 5, controle de chaine cassee par CA–CA max, save_pdb mdtraj, image 0 = reference).
ECARTS DECLARES : (1) chaines renommees A, B sequentiellement (DynaPIN `-ch A,B` ; 2AJF porte A,C) ; (2) un systeme dont le xtc n a pas 1000 images, ou dont une chaine est
cassee, est un REJET (entree au manifeste), jamais un saut ; (3) garde_H (0 H, topologie constante) rejoue sur chaque pile ; (4) sortie sous generation/md_ref/piles/<PDB>__MD__rep<r>/pile.pdb.
Ecrit generation/md_ref/preparation.json. Python du venv bioemu2 (mdtraj)."""
import json, os, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/hydrogenes/orig")
import gen_lib as L
import garde_H
import mdtraj as md, numpy as np

C = L.BENCH / "corpus/systemes"
PILES = L.G / "md_ref/piles"
PILES.mkdir(parents=True, exist_ok=True)
bilan, rej = {"piles": 0}, []
t0 = time.time()
RENUM = "--renumeroter" in sys.argv
SEL = [a for a in sys.argv[1:] if not a.startswith("--")]
for pdb in sorted(L.corpus_28()):
    if SEL and pdb not in SEL:
        continue
    for rep in (1, 2, 3):
        u = f"{pdb}__MD__rep{rep}"
        try:
            d = C / pdb / f"rep{rep}"
            t = md.load(str(d / "trajectory.xtc"), top=str(d / "structure.pdb"))
            if t.n_frames != 1000:
                raise RuntimeError(f"{t.n_frames} images, 1000 attendues")
            t = t[::5]
            chs = list(t.topology.chains)
            if len(chs) != 2:
                raise RuntimeError(f"{len(chs)} chaines (DynaPIN -ch A,B en attend 2)")
            for i, ch in enumerate(chs):
                ch.chain_id = chr(ord("A") + i)
            if RENUM:                                  # numerotation ORDINALE 1..n par chaine (contrat d entree §10 regle 6) : DSSP echoue sur les codes d insertion (3RJQ, 4M3K, 5HGG)
                for ch in t.topology.chains:
                    for o, r in enumerate(ch.residues, start=1):
                        r.resSeq = o
            ca = [a.index for a in t.topology.atoms if a.name == "CA"]
            dca = np.linalg.norm(t.xyz[:, ca[1:]] - t.xyz[:, ca[:-1]], axis=2)
            same = np.array([t.topology.atom(ca[i + 1]).residue.chain.index == t.topology.atom(ca[i]).residue.chain.index for i in range(len(ca) - 1)])
            if float(dca[:, same].max()) > 4.5:
                raise RuntimeError(f"chaine cassee (CA-CA max {float(dca[:, same].max()):.2f} A) : recollage PBC a faire")
            p = PILES / u
            p.mkdir(parents=True, exist_ok=True)
            t.save_pdb(str(p / "pile.pdb"))
            g = garde_H.verifier(str(p / "pile.pdb"))
            if g["code"] != 0:
                raise RuntimeError(f"garde_H code {g['code']} {g['differences']}")
            L.ajouter_manifeste({"type": "pile_MD", "unite": u, "systeme": pdb, "usage": "MD", "regime": f"rep{rep}_stride5", "statut": "ok", "n_modeles": g["n_modeles"], "H_min": g["H_min"],
                                 "atomes_par_modele": g["atomes_par_modele"], "ca_ca_max_A": round(float(dca[:, same].max()), 2), "pile": str(p), "numerotation": "ordinale_1_n_par_chaine" if RENUM else "auteur_(resSeq_de_la_MD)", "jobid": os.environ.get("SLURM_JOB_ID")})
            bilan["piles"] += 1
            print(f"MD_PILE {u} ok {g['n_modeles']} modeles, {g['atomes_par_modele']} atomes", flush=True)
        except Exception as e:
            rej.append({"unite": u, "raison": f"{type(e).__name__}: {str(e)[:250]}"})
            L.ajouter_manifeste({"type": "rejet", "unite": u, "systeme": pdb, "usage": "MD", "statut": "rejet_pile_MD", "motif": f"{type(e).__name__}: {str(e)[:250]}"})
            print(f"MD_PILE REJET {u} : {e}", flush=True)
bilan.update(rejets=rej, n_rejets=len(rej), duree_s=round(time.time() - t0, 1))
(L.G / "md_ref").mkdir(exist_ok=True)
(L.G / ("md_ref/preparation_reprise.json" if SEL else "md_ref/preparation.json")).write_text(json.dumps(bilan, indent=1) + "\n")
print("MD_PILES_FIN", json.dumps({k: v for k, v in bilan.items() if k != "rejets"}), f"(rejets = {len(rej)})")
sys.exit(0 if not rej else 1)
