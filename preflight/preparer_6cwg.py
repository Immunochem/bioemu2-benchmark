#!/usr/bin/env python3
"""preparer.py — 6CWG : x0 (cristal ; image 0 de la replique 1) en atomes lourds, et les 3 repliques MD (stride 5).
Sortie : x0/x0_cristal.pdb, x0/x0_md0.pdb, x0/prep.json, md/rep{1,2,3}/pile.pdb (avec H du champ de force).
Aucun calcul de verite terrain : on ne fait que mettre en forme ; DynaPIN fera les mesures."""
import json, sys
from pathlib import Path
import mdtraj as md
import numpy as np

W = Path("/gpfs/phherman/bioemu2_bench/6CWG")
CRIST = "/gpfs/phherman/evomd-data/Trajectory_Inputs/Trajectory_Inputs/6CWG/6CWG_input.pdb"
DYN = Path("/gpfs/phherman/bioemu2_bench/dynabench/6CWG")


def el(a):
    if a.element is None:
        sys.exit(f"ARRET : element inconnu pour {a}")
    return a.element.symbol


def cles(top):
    out = {}
    for ch in top.chains:
        for o, r in enumerate(ch.residues, start=1):
            for at in r.atoms:
                if el(at) != "H":
                    out[(ch.index, o, at.name)] = at.index
    return out


def main():
    (W / "x0").mkdir(parents=True, exist_ok=True)
    cr = md.load(CRIST)
    m0 = md.load(str(DYN / "rep1/trajectory.xtc"), top=str(DYN / "rep1/structure.pdb"))[0]
    info = {"cristal_fichier": CRIST, "md0_fichier": str(DYN / "rep1/trajectory.xtc") + " image 0"}
    sc, sm = list(cr.topology.to_fasta()), list(m0.topology.to_fasta())
    info["sequences_identiques"] = sc == sm
    info["longueurs"] = [len(s) for s in sc]
    info["noms_residus_cristal"] = sorted({r.name for r in cr.topology.residues})
    info["noms_residus_md0"] = sorted({r.name for r in m0.topology.residues})
    if sc != sm:
        sys.exit("ARRET : les sequences du cristal et de l image 0 different.")
    kc, km = cles(cr.topology), cles(m0.topology)
    communes = sorted(set(kc) & set(km), key=lambda k: kc[k])
    info.update(n_lourds_cristal=len(kc), n_lourds_md0=len(km), n_lourds_communs=len(communes),
                seulement_cristal=sorted(map(str, set(kc) - set(km)))[:10],
                seulement_md0=sorted(map(str, set(km) - set(kc)))[:10],
                n_H_cristal=cr.n_atoms - len(kc), n_H_md0=m0.n_atoms - len(km))
    if len(communes) < 0.99 * len(kc):
        sys.exit(f"ARRET : {len(communes)} atomes lourds communs sur {len(kc)}.")
    x_c = cr.atom_slice([kc[k] for k in communes])
    x_m = m0.atom_slice([km[k] for k in communes])
    x_c.save_pdb(str(W / "x0/x0_cristal.pdb"))
    x_m.save_pdb(str(W / "x0/x0_md0.pdb"))
    # ecart cristal <-> image 0 (indication ; la mesure de reference est celle de DynaPIN, au J4)
    ca = x_c.topology.select("name CA")
    d = {}
    xc, xm = x_c.xyz[0] * 10, x_m.xyz[0] * 10
    def rmsd_sup(idx_fit, idx_eval):
        a, b = xc[idx_fit], xm[idx_fit]
        ca_, cb_ = a - a.mean(0), b - b.mean(0)
        u, s, vt = np.linalg.svd(cb_.T @ ca_)
        dd = np.sign(np.linalg.det(u @ vt)); u[:, -1] *= dd
        r = u @ vt
        mc = a.mean(0); mm = b.mean(0)
        e = (xm[idx_eval] - mm) @ r - (xc[idx_eval] - mc)
        return float(np.sqrt((e ** 2).sum(1).mean()))
    chA = np.array([i for i in ca if x_c.topology.atom(i).residue.chain.index == 0])
    chB = np.array([i for i in ca if x_c.topology.atom(i).residue.chain.index == 1])
    info["ecart_x0_A_Ca_apres_superposition_sur_A"] = round(rmsd_sup(chA, chA), 3)
    info["ecart_x0_B_Ca_apres_superposition_sur_A"] = round(rmsd_sup(chA, chB), 3)
    info["ecart_x0_B_Ca_apres_superposition_sur_B"] = round(rmsd_sup(chB, chB), 3)
    info["ecart_x0_tous_lourds_apres_superposition_globale"] = round(rmsd_sup(np.arange(len(xc)), np.arange(len(xc))), 3)
    # repliques MD, stride 5 (comme l analyse DynaPIN existante), avec H
    info["md"] = {}
    for r in (1, 2, 3):
        d_ = DYN / f"rep{r}"
        t = md.load(str(d_ / "trajectory.xtc"), top=str(d_ / "structure.pdb"))
        if t.n_frames != 1000:
            sys.exit(f"ARRET : rep{r} a {t.n_frames} images, 1000 attendues.")
        t = t[::5]
        top = t.topology
        pairs = []
        for ch in top.chains:
            c = [a.index for a in ch.atoms if a.name == "CA"]
            pairs += list(zip(c[:-1], c[1:]))
        dca = md.compute_distances(t, np.array(pairs)) * 10
        # contact A-B : plus petite distance CA(A)-CA(B)
        a_ = [a.index for a in top.atoms if a.name == "CA" and a.residue.chain.index == 0]
        b_ = [a.index for a in top.atoms if a.name == "CA" and a.residue.chain.index == 1]
        dmin = md.compute_contacts  # (non utilise)
        pp = np.array([(i, j) for i in a_ for j in b_])
        dab = (md.compute_distances(t, pp) * 10).min(1)
        info["md"][f"rep{r}"] = {"images_stride5": int(t.n_frames), "max_CA_CA_A": round(float(dca.max()), 2),
                                 "min_CA_AB_max_sur_images": round(float(dab.max()), 2)}
        if dca.max() > 4.5:
            sys.exit(f"ARRET : rep{r} contient une chaine cassee (CA-CA max {dca.max():.2f} A) : recollage PBC a faire.")
        (W / f"md/rep{r}").mkdir(parents=True, exist_ok=True)
        t.save_pdb(str(W / f"md/rep{r}/pile.pdb"))
    (W / "x0/prep.json").write_text(json.dumps(info, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(info, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
