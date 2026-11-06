#!/usr/bin/env python3
"""inspecter_5lzm.py — que contient 5LZM, atome par atome, avant toute curation ? Puis les CV du cristal (garde ordinal).
Ne modifie rien : lit donnees/pdb5lzm/5LZM.pdb, ecrit resultats/5lzm_inspection.json."""
import collections, json, os, sys
import numpy as np, mdtraj as md
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C
T = "/gpfs/phherman/bioemu2_bench/t4l"
f = f"{T}/donnees/pdb5lzm/5LZM.pdb"
lignes = open(f).read().splitlines()
atm = [l for l in lignes if l.startswith("ATOM")]
het = [l for l in lignes if l.startswith("HETATM")]
alt = [l for l in atm + het if l[16] != " "]
link = [l for l in lignes if l.startswith(("LINK", "SSBOND", "CONECT"))]
occ = collections.Counter(l[54:60].strip() for l in atm)
print("ATOM", len(atm), "HETATM", len(het), "altLoc", len(alt), "LINK/SSBOND/CONECT", len(link))
for l in link[:12]: print("  ", l)
print("occupations :", dict(occ))
# atomes lourds attendus par residu (topologie standard) vs presents
ATTENDU = {"ALA":5,"ARG":11,"ASN":8,"ASP":8,"CYS":6,"GLN":9,"GLU":9,"GLY":4,"HIS":10,"ILE":8,"LEU":8,"LYS":9,"MET":8,"PHE":11,"PRO":7,"SER":6,"THR":7,"TRP":14,"TYR":12,"VAL":7}
par_res = collections.OrderedDict()
for l in atm:
    k = (l[21], int(l[22:26]), l[26], l[17:20])
    par_res.setdefault(k, []).append(l[12:16].strip())
manq = []
for (ch, n, ic, nom), at in par_res.items():
    att = ATTENDU[nom] + (1 if n == 164 else 0)   # OXT du C-terminal
    if len(at) != att: manq.append((n, nom, len(at), att))
print("residus ATOM :", len(par_res), "| residus a nombre d'atomes lourds != attendu :", len(manq))
for m in manq: print("   resSeq %d %s : %d atomes presents, %d attendus" % m)
print("codes d'insertion :", sum(1 for k in par_res if k[2] != " "))
tr = md.load(f)    # mdtraj : topologie des ATOM + HETATM
print("mdtraj :", tr.n_atoms, "atomes,", tr.n_residues, "residus,", tr.n_chains, "chaine(s)")
prot = tr.atom_slice(tr.topology.select("protein"))
r = C.cv(prot)
print("CV 5LZM : d1 = %.4f nm, d2 = %.4f nm, p = %.4f" % (r["d1"][0], r["d2"][0], r["p"][0]))
print("resSeq des ordinaux :", r["rapport"]["resSeq_ordinaux"], "| cycle Phe4 :", r["rapport"]["atomes_cycle_phe4"], "| diffs sequence vs papier :", r["rapport"]["diffs_sequence_vs_papier"])
json.dump({"n_ATOM": len(atm), "n_HETATM": len(het), "n_altloc": len(alt), "HETATM_par_residu": dict(collections.Counter(l[17:20] for l in het)),
           "residus_atomes_manquants": manq, "d1_nm": float(r["d1"][0]), "d2_nm": float(r["d2"][0]), "p": float(r["p"][0])},
          open(f"{T}/resultats/5lzm_inspection.json", "w"), indent=1)
