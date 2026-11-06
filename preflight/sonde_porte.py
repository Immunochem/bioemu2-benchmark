#!/usr/bin/env python3
"""sonde_porte.py — la porte CPU est-elle faisable ? Temps, GPU absent, topologie du contrat.
Aucun modele charge, aucun CUDA. x0 de reference : 6CWG/x0/x0_cristal.pdb (LECTURE SEULE ; fait
passer la porte au J2 du pilote)."""
import os, sys, time, json
os.environ["CUDA_VISIBLE_DEVICES"] = ""
t0 = time.time()
import mdtraj as md
import torch
from bioemu2 import System, map_cif_atom_order
print("import bioemu2 + torch : %.1f s ; torch %s ; cuda dispo = %s" % (time.time() - t0, torch.__version__, torch.cuda.is_available()))
import bioemu2
print("bioemu2 :", os.path.dirname(bioemu2.__file__))
X0 = "/gpfs/phherman/bioemu2_bench/6CWG/x0/x0_cristal.pdb"
ref = md.load(X0)
seqs = list(ref.topology.to_fasta())
t0 = time.time(); system = System.from_description(proteins=seqs); print("from_description : %.1f s ; chain_ids %s" % (time.time() - t0, system.chain_ids))
t0 = time.time(); top = system.topology; print("system.topology : %.1f s ; %d atomes, %d residus, %d chaines" % (time.time() - t0, top.n_atoms, top.n_residues, top.n_chains))
for ch in top.chains:
    rs = list(ch.residues)
    print("  chaine", getattr(ch, "chain_id", None), "n_res", len(rs))
    for r in (rs[0], rs[1], rs[-1]):
        print("    ", r.name, r.resSeq, [a.name for a in r.atoms])
ile = [r for r in top.residues if r.name == "ILE"][:1]
for r in ile: print("  ILE", r.resSeq, [a.name for a in r.atoms])
t0 = time.time(); x = map_cif_atom_order(system, X0); print("map_cif_atom_order : %.1f s ; %d atomes mappes" % (time.time() - t0, x.n_atoms))
print("elements non H dans le contrat :", sorted({a.element.symbol for a in top.atoms}))
print("noms d atomes contrat : OXT present ?", any(a.name == "OXT" for a in top.atoms), "; H presents ?", any(a.element.symbol == "H" for a in top.atoms))
