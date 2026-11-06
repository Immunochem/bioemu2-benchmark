"""J3c — `structure.pdb` est-il l'image 0 du xtc ? (L0a : oui sur 8 repliques pilotes sur 9 ; ici sur les 90.)

Pour chaque replique : ecart brut max (nm), ecart de translation du centroide (A), RMSD atomes lourds apres
superposition de structure.pdb contre l'image 0, et image du xtc la plus proche de structure.pdb (RMSD min sur 1000).
"""
import json
import os
import sys

import mdtraj as md
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import CACHE, ROOT, jload, now_iso  # noqa: E402

sel = jload(ROOT + os.environ.get("CORPUS_SELECTION", "/selection.json"))["selection"]
SUF = os.environ.get("CORPUS_SUFFIX", "")
out = open(CACHE + "/verif_structure_vs_image0%s.jsonl" % SUF, "w", encoding="utf-8")
n = n_ident = 0
for s in sel:
    base = "%s/systemes/%s" % (ROOT, s["pdb"])
    for rep in sorted(x for x in os.listdir(base) if x.startswith("rep")):
        d = "%s/%s" % (base, rep)
        p = md.load(d + "/structure.pdb")
        t = md.load(d + "/trajectory.xtc", top=d + "/structure.pdb")
        heavy = [a.index for a in t.topology.atoms if a.element is not None and a.element.symbol != "H"]
        raw = float(np.abs(t.xyz[0] - p.xyz[0]).max())
        shift = float(np.linalg.norm(t.xyz[0][heavy].mean(axis=0) - p.xyz[0][heavy].mean(axis=0)) * 10)
        r0 = float(md.rmsd(t[0], p, atom_indices=heavy)[0] * 10)
        rall = md.rmsd(t, p, atom_indices=heavy) * 10
        i_min = int(np.argmin(rall))
        rec = {"pdb": s["pdb"], "replique": int(rep[3:]), "ecart_brut_max_nm": round(raw, 4), "translation_centroide_A": round(shift, 2),
               "rmsd_apres_superposition_image0_A": round(r0, 3), "image_la_plus_proche": i_min,
               "rmsd_image_la_plus_proche_A": round(float(rall[i_min]), 3), "date": now_iso()}
        out.write(json.dumps(rec) + "\n")
        print("S0 %s %s ecart_brut_max_nm=%.4f translation_A=%.2f rmsd_sup_image0_A=%.3f plus_proche=image%d (%.3f A)" % (
            s["pdb"], rep, raw, shift, r0, i_min, rall[i_min]), flush=True)
        n += 1
        n_ident += (raw < 0.002)
print("BILAN replicas=%d  structure.pdb==image0 (ecart brut < 0,002 nm)=%d  autres=%d" % (n, n_ident, n - n_ident))
print("FINI s3c %s" % now_iso())
