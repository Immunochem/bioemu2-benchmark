#!/usr/bin/env python3
"""q2_metriques.py — mesures GEOMETRIQUES cheap (atomes lourds, CPU) pour U4/Q2 : distance des centres de masse des deux chaines (Calpha), contacts residu-residu
inter-chaines (CAPRI : tout couple d atomes lourds < 5 A), fnat contre une reference, derive intrinseque de chaque chaine. APPARIEMENT ORDINAL partout (jamais resSeq).
Ce n est PAS DynaPIN : c est un crible. Controle positif obligatoire (`--controle`) : son fnat doit retrouver celui de DynaPIN structure par structure sur les piles du pilote."""
import argparse, csv, json, sys
from pathlib import Path
import mdtraj as md
import numpy as np

SEUIL_A = 5.0


def residus_par_chaine(top):
    """liste, par chaine, des indices d atomes lourds de chaque residu (ordinal dans la chaine)."""
    out = []
    for ch in top.chains:
        out.append([[a.index for a in r.atoms if a.element.symbol != "H"] for r in ch.residues])
    return out


def d_com_ca(traj, ia, ib):
    ca = traj.xyz[:, ia].mean(axis=1); cb = traj.xyz[:, ib].mean(axis=1)
    return np.linalg.norm(ca - cb, axis=1) * 10.0


def indices_chaines(top):
    """indices de TOUS les atomes de chaque chaine (BioEmu2 ne sort que des lourds) : le centre de masse du potentiel de steering est le centre de ces atomes,
    pour que le guidage translate la chaine ENTIERE (un centre des seuls Calpha deplacait les Calpha de 9 A et laissait le reste en arriere)."""
    return [[a.index for a in ch.atoms] for ch in top.chains]


def d_com_atomes(traj, ia, ib):
    return np.linalg.norm(traj.xyz[:, ia].mean(axis=1) - traj.xyz[:, ib].mean(axis=1), axis=1) * 10.0


def indices_ca(top):
    return [[a.index for a in ch.atoms if a.name == "CA"] for ch in top.chains]


def contacts_residus(traj, res, ch_a=0, ch_b=1, seuil=SEUIL_A):
    """pour chaque image : ensemble {(i_res_A, j_res_B)} de couples en contact (atomes lourds < seuil A), ordinaux dans leur chaine."""
    ra, rb = res[ch_a], res[ch_b]
    ia = np.concatenate([np.array(r, dtype=int) for r in ra]); ib = np.concatenate([np.array(r, dtype=int) for r in rb])
    oa = np.concatenate([np.full(len(r), k) for k, r in enumerate(ra)]); ob = np.concatenate([np.full(len(r), k) for k, r in enumerate(rb)])
    out = []
    from scipy.spatial import cKDTree
    for f in range(traj.n_frames):
        xa = traj.xyz[f, ia] * 10.0; xb = traj.xyz[f, ib] * 10.0
        ta, tb = cKDTree(xa), cKDTree(xb)
        pairs = ta.query_ball_tree(tb, seuil)
        s = set()
        for i, lst in enumerate(pairs):
            for j in lst:
                s.add((int(oa[i]), int(ob[j])))
        out.append(s)
    return out


def residus_interface(ref, res, seuil=10.0):
    """CAPRI/DockQ : residus de chaque chaine ayant un atome lourd a moins de `seuil` A d un atome lourd de l AUTRE chaine, dans la reference. Renvoie {chaine: [ordinaux]}."""
    from scipy.spatial import cKDTree
    out = {}
    ia = np.concatenate([np.array(r, dtype=int) for r in res[0]]); ib = np.concatenate([np.array(r, dtype=int) for r in res[1]])
    oa = np.concatenate([np.full(len(r), k) for k, r in enumerate(res[0])]); ob = np.concatenate([np.full(len(r), k) for k, r in enumerate(res[1])])
    xa, xb = ref.xyz[0, ia] * 10.0, ref.xyz[0, ib] * 10.0
    ta, tb = cKDTree(xa), cKDTree(xb)
    pa = ta.query_ball_tree(tb, seuil)
    out[0] = sorted({int(oa[i]) for i, l in enumerate(pa) if l})
    pb = tb.query_ball_tree(ta, seuil)
    out[1] = sorted({int(ob[i]) for i, l in enumerate(pb) if l})
    return out


def indices_squelette(top, residus_par_chaine, noms=("N", "CA", "C", "O")):
    """indices des atomes de squelette des residus donnes {chaine: [ordinaux]}, dans l ORDRE CANONIQUE (chaine, ordinal, nom) : deux topologies dont l ordre d atomes
    differe (ordre PDB contre ordre BioEmu2) donnent deux listes de MEME sens, jamais deux listes triees par indice."""
    idx = []
    for ci in sorted(residus_par_chaine):
        residues = list(list(top.chains)[ci].residues)
        for o in residus_par_chaine[ci]:
            par_nom = {at.name: at.index for at in residues[o].atoms}
            for n in noms:
                if n in par_nom:
                    idx.append(par_nom[n])
                else:
                    raise ValueError(f"residu d ordinal {o} de la chaine {ci} sans atome {n}")
    return idx


def kabsch_rmsd(X, Y):
    """RMSD apres superposition optimale, X [n_frames, n, 3], Y [n, 3] (memes unites)."""
    Yc = Y - Y.mean(axis=0)
    out = np.empty(len(X))
    for f in range(len(X)):
        Xc = X[f] - X[f].mean(axis=0)
        U, S, Vt = np.linalg.svd(Xc.T @ Yc)
        d = np.sign(np.linalg.det(U @ Vt))
        D = np.diag([1.0, 1.0, d])
        R = U @ D @ Vt
        out[f] = np.sqrt(np.mean(np.sum((Xc @ R - Yc) ** 2, axis=1)))
    return out


def irmsd(traj, ref, res_interface, noms=("N", "CA", "C", "O")):
    """RMSD (A) des atomes de squelette des residus d interface de la REFERENCE (appariement ordinal), superposes sur eux-memes (definition CAPRI)."""
    it = indices_squelette(traj.topology, res_interface, noms)
    ir = indices_squelette(ref.topology, res_interface, noms)
    return kabsch_rmsd(traj.xyz[:, it] * 10.0, ref.xyz[0, ir] * 10.0)


def fnat(contacts_modele, contacts_ref):
    nat = contacts_ref
    return [len(c & nat) / len(nat) if nat else float("nan") for c in contacts_modele]


def derive_chaine(traj, ref, idx_ca):
    return (md.rmsd(traj, ref, atom_indices=idx_ca) * 10.0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--controle", action="store_true", help="CONTROLE POSITIF : fnat de ce module contre celui de DynaPIN sur les piles du pilote (cristal_lag0050, ref image 0 MD)")
    a = ap.parse_args()
    if a.controle:
        W = Path("/gpfs/phherman/bioemu2_bench/6CWG")
        ref = md.load(str(W / "x0/x0_md0.pdb"))
        d = W / "prop/cristal_lag0050/brut"
        tr = md.load_xtc(str(d / "samples.xtc"), top=str(d / "topology.cif"))
        resR, resT = residus_par_chaine(ref.topology), residus_par_chaine(tr.topology)
        assert [len(c) for c in resR] == [len(c) for c in resT], "chaines de longueurs differentes"
        cr = contacts_residus(ref, resR)[0]
        ct = contacts_residus(tr, resT)
        fm = np.array(fnat(ct, cr))
        rows = list(csv.DictReader(open("/gpfs/phherman/bioemu2_bench/hydrogenes/out/R6_all/piles/cristal_vsMD0__cristal_lag0050/dynapin/tables/qc_trajectory_metrics.csv")))
        rows.sort(key=lambda r: int(float(r["Frame"])))
        fd = np.array([float(r["fnat"]) for r in rows[1:]])
        print(f"contacts natifs de la reference (image 0 MD) : {len(cr)} couples residu-residu < {SEUIL_A} A")
        print(f"fnat ce module : mediane {np.median(fm):.3f} · DynaPIN : mediane {np.median(fd):.3f}")
        dd = fm - fd
        print(f"ecart par structure : mediane {np.median(dd):+.4f}, min {dd.min():+.4f}, max {dd.max():+.4f}, |ecart| max {np.abs(dd).max():.4f} sur {len(dd)} structures")
        ri = residus_interface(ref, resR)
        im = irmsd(tr, ref, ri)
        idd = np.array([float(r["iRMSD"]) for r in rows[1:]])
        di = im - idd
        print(f"interface (10 A) de la reference : {len(ri[0])} residus (chaine A) + {len(ri[1])} (chaine B)")
        print(f"iRMSD ce module : mediane {np.median(im):.3f} · DynaPIN : mediane {np.median(idd):.3f} · ecart par structure : mediane {np.median(di):+.4f}, |ecart| max {np.abs(di).max():.4f}")
        ia, ib = indices_ca(tr.topology)[0], indices_ca(tr.topology)[1]
        print(f"d_COM Calpha : reference {d_com_ca(ref, indices_ca(ref.topology)[0], indices_ca(ref.topology)[1])[0]:.2f} A ; pilote cristal_lag0050 mediane {np.median(d_com_ca(tr, ia, ib)):.2f} A (min {d_com_ca(tr, ia, ib).min():.2f}, max {d_com_ca(tr, ia, ib).max():.2f})")
        sys.exit(0 if (np.abs(dd).max() < 0.02 and np.abs(di).max() < 0.05) else 1)
