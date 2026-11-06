#!/usr/bin/env python3
"""h_diag_ss.py — avant/pendant l etage H : combien d images feraient echouer l outil de reconstruction a cause des ponts ?
OpenMM (Topology.createDisulfideBonds, seuil 0,3 nm) INFERE une liaison S–S entre deux SG a < 3,0 A, quelle que soit la declaration de la MD.
  cas A : deux cysteines LIBRES (CYS dans le PSF) dont les SG sont a < 3,0 A  -> liaison inferee, variante CYS : « No template … 1 S atom too many »
  cas B : un pont DECLARE (CYX dans le PSF) dont les SG sont a > 3,0 A        -> aucune liaison inferee, variante CYX : template absent
Compte par systeme / unite / image, sur les xtc de generation/sorties (U3). Ecrit hydrogenes/diag_ss.json. CPU seul."""
import json, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L
import mdtraj as md, numpy as np

HR, US = os.environ.get("H_RACINE", "hydrogenes"), os.environ.get("H_USAGE", "U3")
ref = json.loads((L.G / "ponts/reference.json").read_text())["systemes"]
out, tot = {}, {"A": 0, "B": 0, "images": 0, "images_touchees": 0}
for pdb in sorted(L.corpus_28()):
    x0 = md.load(str(L.chemin_x0(pdb, 1)))
    ch = {chr(ord("A") + i): c for i, c in enumerate(x0.topology.chains)}
    sg = {}
    for cl, c in ch.items():
        for o, r in enumerate(c.residues, start=1):
            if r.name == "CYS":
                sg[(cl, o)] = [a.index for a in r.atoms if a.name == "SG"][0]
    ponts = {(p["chaine_i"], p["ordinal_i"], p["chaine_j"], p["ordinal_j"]) for p in ref[pdb]["ponts"]}
    lies = {(a, b) for (a, i, b, j) in ponts} and set()
    lies = set()
    for (a, i, b, j) in ponts:
        lies |= {(a, i), (b, j)}
    libres = [k for k in sg if k not in lies]
    d = {"n_cys": len(sg), "n_ponts": len(ponts), "n_libres": len(libres), "unites": {}}
    for ud in sorted((L.G / "sorties").glob(f"{pdb}__{US}__*")):
        tr = md.load_xtc(str(ud / "samples.xtc"), top=str(ud / "topology.cif"))
        # indices SG de l echantillon par ordinal
        ch2 = {chr(ord("A") + i): c for i, c in enumerate(tr.topology.chains)}
        sg2 = {(cl, o): [a.index for a in r.atoms if a.name == "SG"][0] for cl, c in ch2.items() for o, r in enumerate(c.residues, start=1) if r.name == "CYS"}
        A, B, touche = [], [], set()
        for f in range(tr.n_frames):
            for a in range(len(libres)):
                for b in range(a + 1, len(libres)):
                    dist = float(np.linalg.norm(tr.xyz[f, sg2[libres[a]]] - tr.xyz[f, sg2[libres[b]]]) * 10)
                    if dist < 3.0:
                        A.append((f, libres[a], libres[b], round(dist, 2))); touche.add(f)
            for (a, i, b, j) in ponts:
                dist = float(np.linalg.norm(tr.xyz[f, sg2[(a, i)]] - tr.xyz[f, sg2[(b, j)]]) * 10)
                if dist > 3.0:
                    B.append((f, (a, i), (b, j), round(dist, 2))); touche.add(f)
        d["unites"][ud.name] = {"n_images": tr.n_frames, "cas_A_libres_proches": A, "cas_B_ponts_rompus": B, "images_touchees": sorted(touche)}
        tot["A"] += len(A); tot["B"] += len(B); tot["images"] += tr.n_frames; tot["images_touchees"] += len(touche)
    out[pdb] = d
    tp = sum(len(u["images_touchees"]) for u in d["unites"].values())
    print(f"DIAG_SS {pdb} cys={d['n_cys']} ponts={d['n_ponts']} libres={d['n_libres']} images touchees {tp}/{sum(u['n_images'] for u in d['unites'].values())} "
          f"(A {sum(len(u['cas_A_libres_proches']) for u in d['unites'].values())}, B {sum(len(u['cas_B_ponts_rompus']) for u in d['unites'].values())})", flush=True)
(L.G / HR / "diag_ss.json").write_text(json.dumps({"total": tot, "systemes": out}, ensure_ascii=False) + "\n")
print("DIAG_SS_TOTAL", json.dumps(tot))
