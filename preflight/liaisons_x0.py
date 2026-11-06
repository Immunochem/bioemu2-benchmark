#!/usr/bin/env python3
"""liaisons_x0.py — combien de x0 CURES portent une jonction C–N ETIREE (1,6 a 2,9 A), acceptee par la fenetre 2,9 A ?

curer_entree.py refuse une rupture (> 2,9 A) mais ne COMPTE pas a part les jonctions etirees qu'il accepte (boucles
reconstruites des structures MD : 2,01 / 2,11 / 2,40 / 2,55 A mesurees sur le socle, RAPPORT_PREP.md §3). Ce script relit
les 178 x0 ecrits (curation/x0/*.pdb), avec `residus_de` de preparer_cristal (importe via curer_entree), et releve la plus
longue liaison C–N de chaque chaine. Lecture seule. Ne corrige rien.
"""
import collections
import glob
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
os.environ["CUDA_VISIBLE_DEVICES"] = ""
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/prep_fix")
import curer_entree as ce
import mdtraj as md

prep = ce.prep
X0 = Path("/gpfs/phherman/bioemu2_bench/prep_fix/curation/x0")


def main():
    lignes, par_sys = [], collections.defaultdict(list)
    for f in sorted(glob.glob(str(X0 / "*.pdb"))):
        nom = Path(f).stem
        pdb, kind = nom.split("__")
        t = md.load(f)
        top = t.topology
        lettres = list(dict.fromkeys(prep.chaine_id(top, c) for c in top.chains))
        etire, vmax = [], 0.0
        for L in lettres:
            res = prep.residus_de(t, L)
            rs, d = res["residus"], res["d_suivant_nm"]
            for k in range(len(rs) - 1):
                dd = d[k] * 10
                vmax = max(vmax, dd)
                if dd > 1.6:
                    etire.append((L, rs[k].resSeq, round(dd, 3)))
        par_sys[pdb].append((kind, vmax, etire))
        lignes.append((pdb, kind, vmax, etire))
    print(f"x0 relus : {len(lignes)} ; avec au moins une liaison C–N > 1,6 A : {sum(1 for l in lignes if l[3])} ; > 1,9 A : {sum(1 for l in lignes if l[2] > 1.9)} ; max global : {max(l[2] for l in lignes):.3f} A")
    print(f"{'systeme':<8} {'x0':>3} {'max C–N (A)':>12}  x0 avec jonction > 1,6 A (kind: nb ; plus longue)")
    for pdb, v in par_sys.items():
        et = [(k, len(e), max(x[2] for x in e)) for k, m, e in v if e]
        print(f"{pdb:<8} {len(v):>3} {max(m for k, m, e in v):>12.3f}  {et if et else '-'}")
    tous = sorted((l for l in lignes if l[3]), key=lambda l: -l[2])
    print("\nTOUTES les jonctions > 1,6 A (chaine du x0, ordinal, d en A), x0 les plus etirees d'abord :")
    for pdb, kind, vmax, e in tous[:40]:
        print(f"  {pdb}__{kind}: {e}")
    print("\nhistogramme des liaisons C–N les plus longues par x0 (A) :")
    bins = [(0, 1.5), (1.5, 1.6), (1.6, 1.9), (1.9, 2.4), (2.4, 2.9), (2.9, 99)]
    for lo, hi in bins:
        print(f"  [{lo} ; {hi}) : {sum(1 for l in lignes if lo <= l[2] < hi)} x0")


if __name__ == "__main__":
    main()
