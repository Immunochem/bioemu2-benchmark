#!/usr/bin/env python3
"""md_contre_x0.py — DockQ / fnat / iRMSD de CHAQUE image de la MD de reference contre le MEME
x0 cure que les piles BioEmu2, pour rendre le balayage de lag comparable.

POURQUOI : DynaPIN note chaque image contre l image 0 de SA pile. Les piles BioEmu2 portent le
x0 cure en image 0 (h_preparer.py:3) ; les piles MD portent leur propre 1re image
(md_piles.py:32-48). Les deux series ne sont pas sur la meme echelle.

AUCUNE REIMPLEMENTATION : on appelle redock_lib.Reference(...).mesurer(...), le chemin fige que
le re-docking utilise et qui retrouve DynaPIN a 0,0005 pres. Stride 5 identique a md_piles.py,
donc 200 images par replique, 0,5 ns par image.

🔴 Appariement ORDINAL (index de chaine, ordinal du residu dans la chaine, nom d atome).
🔴 Tout desaccord d atomes lourds est un REJET inscrit au bilan, jamais un saut silencieux.

ECART DECLARE, un seul : la topologie du champ de force porte `OXT` (oxygene terminal du
carboxyle) a la fin de chaque chaine, le x0 cure ne l a pas. Mesure du 2026-11-06 : l ecart
vaut EXACTEMENT 2 atomes (un par chaine) sur 6CWG, 1AY7 et 1GCQ. On retire donc les `OXT`
du cote MD, on les COMPTE, et on mesure leur distance minimale a l autre chaine : un OXT a
moins de 5 A de l autre chaine serait un atome d interface, et son retrait deplacerait fnat.
Tout autre desaccord reste un REJET.
"""
import json, os, sys, time
from pathlib import Path
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/generation/scripts")
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/redocking/scripts")
import gen_lib as L
import redock_lib as RL
import mdtraj as md
import numpy as np

C = L.BENCH / "corpus/systemes"
SORTIE = L.G / "md_contre_x0"
SORTIE.mkdir(parents=True, exist_ok=True)
NS_PAR_IMAGE = 0.5                      # stride 5 sur 1000 images de 100 ns (md_piles.py:2)
SEL = [a for a in sys.argv[1:] if not a.startswith("--")]


def cles(top):
    """cle ORDINALE -> index d atome, atomes lourds seuls. Copie de h_preparer.cles."""
    out = {}
    for ch in top.chains:
        for o, r in enumerate(ch.residues, start=1):
            for at in r.atoms:
                if at.element.symbol != "H":
                    out[(ch.index, o, at.name)] = at.index
    return out


bilan = {"systemes": 0, "repliques": 0, "images": 0, "rejets": [], "ns_par_image": NS_PAR_IMAGE}
lignes, t0 = [], time.time()
n_oxt, oxt_dmin = 0, []

for pdb in sorted(L.corpus_28()):
    if SEL and pdb not in SEL:
        continue
    try:
        x0 = md.load(str(L.chemin_x0(pdb, 1)))
        L.verifier_x0(pdb, 1)
        kref = cles(x0.topology)
        ordre = sorted(kref, key=lambda c: kref[c])
        i_ref = [kref[c] for c in ordre]
        x0h = x0.atom_slice(i_ref)
        R = RL.Reference(x0h)
    except Exception as e:
        bilan["rejets"].append({"systeme": pdb, "etape": "x0", "raison": f"{type(e).__name__}: {str(e)[:200]}"})
        print(f"REJET {pdb} x0 : {e}", flush=True)
        continue

    nrep = 0
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
                raise RuntimeError(f"{len(chs)} chaines, 2 attendues")
            for i, ch in enumerate(chs):
                ch.chain_id = chr(ord("A") + i)
            for ch in t.topology.chains:
                for o, r in enumerate(ch.residues, start=1):
                    r.resSeq = o
            ke = cles(t.topology)
            oxt = sorted(k for k in set(ke) - set(kref) if k[2] == "OXT")
            # un OXT d interface deplacerait fnat : on mesure avant de retirer
            if oxt:
                X = t.xyz[0] * 10.0
                ch_at = {}
                for ch in t.topology.chains:
                    ch_at[ch.index] = [a.index for a in ch.atoms if a.element.symbol != "H"]
                for k in oxt:
                    autre = [i for ci, ids in ch_at.items() if ci != k[0] for i in ids]
                    d = float(np.min(np.linalg.norm(X[autre] - X[ke[k]], axis=1)))
                    oxt_dmin.append({"unite": u, "chaine": k[0], "residu_ordinal": k[1], "dmin_A": round(d, 2)})
                for k in oxt:
                    ke.pop(k)
                n_oxt += len(oxt)
            if set(ke) != set(kref):
                raise RuntimeError(f"{len(set(ke) ^ set(kref))} atomes lourds en desaccord ordinal avec le x0 "
                                   f"apres retrait de {len(oxt)} OXT (x0 {len(kref)}, MD {len(ke)})")
            th = t.atom_slice([ke[c] for c in ordre])
            m = R.mesurer(th)
            for i in range(th.n_frames):
                lignes.append({"systeme": pdb, "replique": f"rep{rep}", "image": i,
                               "ns": round(i * NS_PAR_IMAGE, 2), "dockq": m["dockq"][i],
                               "fnat": m["fnat"][i], "irmsd_A": m["irmsd_A"][i], "lrmsd_A": m["lrmsd_A"][i]})
            bilan["images"] += th.n_frames
            nrep += 1
            print(f"OK {u} : {th.n_frames} images, DockQ {m['dockq'][0]:.3f} -> {m['dockq'][-1]:.3f} "
                  f"({time.time()-t0:.0f} s cumulees)", flush=True)
        except Exception as e:
            bilan["rejets"].append({"unite": u, "systeme": pdb, "etape": "mesure",
                                    "raison": f"{type(e).__name__}: {str(e)[:200]}"})
            print(f"REJET {u} : {e}", flush=True)
    bilan["repliques"] += nrep
    if nrep:
        bilan["systemes"] += 1

bilan["oxt_retires"] = n_oxt
bilan["oxt_dmin_min_A"] = round(min((o["dmin_A"] for o in oxt_dmin), default=-1.0), 2)
bilan["oxt_a_moins_de_5A"] = sum(1 for o in oxt_dmin if o["dmin_A"] < 5.0)
bilan["oxt_detail"] = oxt_dmin
bilan["n_rejets"] = len(bilan["rejets"])
bilan["duree_s"] = round(time.time() - t0, 1)
suff = "_smoke" if SEL else ""
(SORTIE / f"mesures{suff}.json").write_text(json.dumps(lignes) + "\n")
(SORTIE / f"bilan{suff}.json").write_text(json.dumps(bilan, indent=1) + "\n")
print("FIN", json.dumps({k: v for k, v in bilan.items() if k not in ("rejets", "oxt_detail")}),
      f"(REJETS = {bilan['n_rejets']} · OXT retires = {n_oxt}, dont {bilan['oxt_a_moins_de_5A']} a moins de 5 A de l autre chaine)")
sys.exit(0 if bilan["n_rejets"] == 0 else 1)
