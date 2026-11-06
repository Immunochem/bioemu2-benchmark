#!/usr/bin/env python3
"""redock_bioemu.py — passerelle vers bioemu2 (CPU ou GPU) : du x0 CURE au System, a ses ponts disulfure DECLARES et a la reference dans l ORDRE d atomes du System.
Ponts : paires SG-SG < 2,5 A dans le x0 cristallin (meme regle que 6CWG/scripts/propager.py::systeme_ss, docs/system.md « Declare every disulfide »), intra OU inter-chaine ;
un pont INTER-chaine est signale (il fixerait l interface par construction)."""
import sys
import numpy as np
import mdtraj as md
from bioemu2 import Atom, Bond, Protein, System, map_cif_atom_order

sys.dont_write_bytecode = True


def systeme_depuis_x0(x0):
    ref = md.load(str(x0))
    seqs = list(ref.topology.to_fasta())
    nu = System.from_description(proteins=seqs)
    if list(nu.chain_ids)[:2] != ["A", "B"]:
        raise ValueError(f"identifiants de chaine inattendus {list(nu.chain_ids)}")
    top = ref.topology
    sg = []
    for ch in top.chains:
        for o, r in enumerate(ch.residues, start=1):
            if r.name == "CYS":
                for a in r.atoms:
                    if a.name == "SG":
                        sg.append((ch.index, o, a.index))
    ponts, inter = [], 0
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            if np.linalg.norm(ref.xyz[0][sg[i][2]] - ref.xyz[0][sg[j][2]]) * 10.0 < 2.5:
                ponts.append((sg[i][:2], sg[j][:2]))
                inter += sg[i][0] != sg[j][0]
    if not ponts:
        systeme = nu
    else:
        bonds = [Bond(atom1=Atom(name="SG", res_name="CYS", res_id=a[1], asym_id="AB"[a[0]]), atom2=Atom(name="SG", res_name="CYS", res_id=b[1], asym_id="AB"[b[0]])) for a, b in ponts]
        systeme = System([Protein(sequence=e.sequence, asym_ids=list(e.asym_ids)) for e in nu.entities], bonds=bonds)
    ref_sys = map_cif_atom_order(systeme, str(x0))
    return {"nu": nu, "systeme": systeme, "sequences": seqs, "ponts": [f"{'AB'[a[0]]}{a[1]}-{'AB'[b[0]]}{b[1]}" for a, b in ponts], "n_ponts_inter_chaine": int(inter), "ref": ref_sys,
            "tokens": sum(len(s) for s in seqs)}
