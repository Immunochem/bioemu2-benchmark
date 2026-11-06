#!/usr/bin/env python3
"""redock_lib.py — leurres, plancher du hasard et mesures CAPRI du redocking BioEmu2.

MESURES : on IMPORTE `q2_metriques_vendu.py` (copie de steering/scripts/q2_metriques.py, sha256 epingle ci-dessous ; contacts residu-residu < 5 A sur atomes lourds, interface a 10 A,
iRMSD sur le squelette N-CA-C-O de l interface de la REFERENCE, Kabsch, APPARIEMENT ORDINAL partout). On y AJOUTE seulement ce qu il n a pas : lRMSD, DockQ, classes CAPRI.
  - fnat    = |contacts_modele ∩ contacts_ref| / |contacts_ref|  (couples de RESIDUS ordinaux, atomes lourds < 5 A)
  - iRMSD   = RMSD du squelette des residus d interface de la reference (atome lourd < 10 A de l autre chaine), superposes sur eux-memes
  - lRMSD   = on superpose le squelette du RECEPTEUR (la plus grande chaine) sur la reference, RMSD du squelette du LIGAND (la plus petite) sans nouvel ajustement
  - DockQ   = (fnat + 1/(1+(iRMSD/1,5)^2) + 1/(1+(lRMSD/8,5)^2)) / 3  (Basu & Wallner 2016) ; classes DockQ : <0,23 incorrect, <0,49 acceptable, <0,80 medium, sinon high
  - classe CAPRI (criteres d origine) : high fnat>=0,5 et (lRMSD<=1 ou iRMSD<=1) ; medium fnat>=0,3 et (lRMSD<=5 ou iRMSD<=2) ; acceptable fnat>=0,1 et (lRMSD<=10 ou iRMSD<=4) ; sinon incorrect
  - d_COM   = distance entre les CENTRES GEOMETRIQUES des atomes lourds des deux chaines (convention du smoke : « d_com_atomes »), en A
LEURRES : la plus PETITE chaine est deplacee (la grande reste en place). Distance = d_COM native + delta le long de l axe natif (centre fixe -> centre mobile).
  T  : translation pure (orientation relative native conservee : information donnee gratuitement)
  R  : translation + rotation uniforme aleatoire de la chaine mobile autour de son centre ; sans choc (distance atomique inter-chaines >= 3,0 A) ou rejetee
  RD : rotation uniforme aleatoire + direction aleatoire (uniforme sur la sphere), meme distance d_COM : aucune orientation ni direction native
  C0 : le cristal lui-meme (delta 0) — CONTROLE PLAFOND : ce que le modele fait d une pose deja juste
PLANCHER DU HASARD (sans modele) : F1 = placement aleatoire (rotation + direction) a d_COM = native, sans choc ; F2 = placement aleatoire amene AU CONTACT (distance atomique min 3,5 A).
"""
import hashlib, importlib.util, sys, zlib
from pathlib import Path
import numpy as np
import mdtraj as md
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

sys.dont_write_bytecode = True
RACINE = Path("/gpfs/phherman/bioemu2_bench/redocking")
SHA_Q2 = "a691cbb93ffbf6f4341372491d0df08ceb1ad3e4042015112e945ebef78efa75"


def _charger_q2():
    p = RACINE / "scripts/q2_metriques_vendu.py"
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    if sha != SHA_Q2:
        raise SystemExit(f"ARRET : {p} a change (sha256 {sha[:16]}, epingle {SHA_Q2[:16]}) — les mesures importent ce module, elles ne suivent pas un changement en silence.")
    sp = importlib.util.spec_from_file_location("q2_metriques_vendu", p)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


Q = _charger_q2()
DMIN_SANS_CHOC_A = 3.0
DMIN_CONTACT_A = 3.5


def seed_de(*parts):
    return zlib.crc32(":".join(str(p) for p in parts).encode()) % (2 ** 31)


# ───────────────────────── mesures ─────────────────────────
def kabsch(P, T):
    """rotation et translation tels que (P - moyenne(P)) @ R + moyenne(T) approche T. P, T : (n, 3)."""
    Pm, Tm = P.mean(axis=0), T.mean(axis=0)
    U, S, Vt = np.linalg.svd((P - Pm).T @ (T - Tm))
    d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1.0, 1.0, d]) @ Vt, Pm, Tm


def classe_capri(f, ir, lr):
    if f >= 0.5 and (lr <= 1.0 or ir <= 1.0):
        return "high"
    if f >= 0.3 and (lr <= 5.0 or ir <= 2.0):
        return "medium"
    if f >= 0.1 and (lr <= 10.0 or ir <= 4.0):
        return "acceptable"
    return "incorrect"


def classe_dockq(dq):
    return "incorrect" if dq < 0.23 else "acceptable" if dq < 0.49 else "medium" if dq < 0.80 else "high"


class Reference:
    """Structure de reference (cristal) a 2 chaines, dans l ORDRE d atomes du System BioEmu2 (sortie de map_cif_atom_order)."""

    def __init__(self, ref):
        self.ref, top = ref, ref.topology
        if top.n_chains != 2:
            raise ValueError(f"2 chaines attendues, {top.n_chains} trouvees")
        self.res = Q.residus_par_chaine(top)
        self.contacts = Q.contacts_residus(ref, self.res)[0]
        self.ri = Q.residus_interface(ref, self.res)
        self.idx = Q.indices_chaines(top)[:2]
        self.ca = Q.indices_ca(top)[:2]
        n0, n1 = len(self.res[0]), len(self.res[1])
        self.rec = 0 if n0 >= n1 else 1
        self.lig = 1 - self.rec            # la plus petite chaine = ligand = chaine DEPLACEE dans les leurres
        sq = lambda c: Q.indices_squelette(top, {c: list(range(len(self.res[c])))})
        self.bb_rec, self.bb_lig = sq(self.rec), sq(self.lig)
        self.d_nat = float(Q.d_com_atomes(ref, *self.idx)[0])
        self.n_contacts_nat = len(self.contacts)
        self.xyz_A = ref.xyz[0] * 10.0

    def lrmsd(self, traj):
        X = traj.xyz * 10.0
        out = np.empty(traj.n_frames)
        Rr, Rl = self.xyz_A[self.bb_rec], self.xyz_A[self.bb_lig]
        for f in range(traj.n_frames):
            Rot, Pm, Tm = kabsch(X[f][self.bb_rec], Rr)
            L = (X[f][self.bb_lig] - Pm) @ Rot + Tm
            out[f] = np.sqrt(np.mean(np.sum((L - Rl) ** 2, axis=1)))
        return out

    def mesurer(self, traj):
        """traj : md.Trajectory de MEME topologie que la reference (memes atomes dans le meme ordre). Rend un dict de listes (une valeur par image)."""
        top = traj.topology
        if [a.name for a in top.atoms] != [a.name for a in self.ref.topology.atoms]:
            raise ValueError("ordre/noms d atomes differents de la reference : appariement ordinal impossible")
        ct = Q.contacts_residus(traj, Q.residus_par_chaine(top))
        fn = np.array(Q.fnat(ct, self.contacts))
        ir = Q.irmsd(traj, self.ref, self.ri)
        lr = self.lrmsd(traj)
        dq = (fn + 1.0 / (1.0 + (ir / 1.5) ** 2) + 1.0 / (1.0 + (lr / 8.5) ** 2)) / 3.0
        X = traj.xyz * 10.0
        dmin = [float(cKDTree(X[f][self.idx[1]]).query(X[f][self.idx[0]])[0].min()) for f in range(traj.n_frames)]
        return {"fnat": fn.round(4).tolist(), "irmsd_A": ir.round(3).tolist(), "lrmsd_A": lr.round(3).tolist(), "dockq": dq.round(4).tolist(),
                "classe_capri": [classe_capri(a, b, c) for a, b, c in zip(fn, ir, lr)], "classe_dockq": [classe_dockq(x) for x in dq],
                "d_com_A": Q.d_com_atomes(traj, *self.idx).round(3).tolist(), "n_contacts": [len(c) for c in ct], "dmin_inter_A": [round(x, 3) for x in dmin],
                "derive_A_ca_A": Q.derive_chaine(traj, self.ref, self.ca[0]).round(3).tolist(), "derive_B_ca_A": Q.derive_chaine(traj, self.ref, self.ca[1]).round(3).tolist()}


# ───────────────────────── leurres ─────────────────────────
def _axe(Rf):
    X = Rf.xyz_A
    cf, cm = X[Rf.idx[Rf.rec]].mean(axis=0), X[Rf.idx[Rf.lig]].mean(axis=0)
    u = (cm - cf) / np.linalg.norm(cm - cf)
    return cf, cm, u


def _traj(Rf, X_A):
    return md.Trajectory(xyz=(X_A / 10.0).astype(np.float32)[None] if X_A.ndim == 2 else (X_A / 10.0).astype(np.float32), topology=Rf.ref.topology)


def dmin_inter(Rf, X_A):
    return float(cKDTree(X_A[Rf.idx[Rf.rec]]).query(X_A[Rf.idx[Rf.lig]])[0].min())


def leurre(Rf, arm, delta, rng, max_essais=500):
    """-> (X en A [n_atomes, 3] ou None, dict d audit). arm in C0, T, R, RD. Le dict porte n_essais, dmin_inter_A, d_com_A, accepte."""
    X = Rf.xyz_A.copy()
    cf, cm, u = _axe(Rf)
    m = Rf.idx[Rf.lig]
    if arm == "C0":
        return X, {"n_essais": 1, "accepte": True, "dmin_inter_A": dmin_inter(Rf, X), "d_com_A": Rf.d_nat}
    if arm == "T":
        X[m] += u * delta
        return X, {"n_essais": 1, "accepte": True, "dmin_inter_A": dmin_inter(Rf, X), "d_com_A": Rf.d_nat + delta}
    for k in range(1, max_essais + 1):
        Rot = Rotation.random(random_state=rng).as_matrix()
        if arm == "R":
            cible = cm + u * delta
        elif arm == "RD":
            v = rng.normal(size=3); v /= np.linalg.norm(v)
            cible = cf + v * (Rf.d_nat + delta)
        else:
            raise ValueError(arm)
        Y = X.copy()
        Y[m] = (X[m] - cm) @ Rot.T + cible
        d = dmin_inter(Rf, Y)
        if d >= DMIN_SANS_CHOC_A:
            dc = float(np.linalg.norm(Y[m].mean(axis=0) - cf))
            return Y, {"n_essais": k, "accepte": True, "dmin_inter_A": d, "d_com_A": dc}
    return None, {"n_essais": max_essais, "accepte": False, "dmin_inter_A": None, "d_com_A": None}


def placement_au_contact(Rf, rng, pas=2.0):
    """F2 : rotation aleatoire, direction aleatoire, la chaine mobile (centre sur la droite centre fixe -> direction) est amenee jusqu a une distance atomique min de 3,5 A."""
    X = Rf.xyz_A
    cf, cm = X[Rf.idx[Rf.rec]].mean(axis=0), X[Rf.idx[Rf.lig]].mean(axis=0)
    m = Rf.idx[Rf.lig]
    Rot = Rotation.random(random_state=rng).as_matrix()
    v = rng.normal(size=3); v /= np.linalg.norm(v)
    arbre = cKDTree(X[Rf.idx[Rf.rec]])
    base = (X[m] - cm) @ Rot.T
    f = lambda t: float(arbre.query(base + cf + v * t)[0].min())
    t_prec = Rf.d_nat + 120.0
    if f(t_prec) < DMIN_CONTACT_A:
        return None
    t = t_prec
    while t > 0:
        t -= pas
        if f(t) < DMIN_CONTACT_A:
            lo, hi = t, t + pas            # f(lo) < 3,5 <= f(hi)
            for _ in range(14):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f(mid) < DMIN_CONTACT_A else (lo, mid)
            Y = X.copy(); Y[m] = base + cf + v * hi
            return Y
    return None


def plancher(Rf, n1, n2, graine):
    """-> {'F1': mesures, 'F2': mesures, 'audit': ...}. Memes metriques que le modele, pour la meme reference."""
    rng = np.random.default_rng(graine)
    X1, essais1 = [], 0
    while len(X1) < n1 and essais1 < 40 * n1:
        essais1 += 1
        Y, au = leurre(Rf, "RD", 0.0, rng, max_essais=1)
        if Y is not None:
            X1.append(Y)
    X2, essais2 = [], 0
    while len(X2) < n2 and essais2 < 10 * n2:
        essais2 += 1
        Y = placement_au_contact(Rf, rng)
        if Y is not None:
            X2.append(Y)
    out = {"audit": {"F1_acceptes": len(X1), "F1_essais": essais1, "F2_acceptes": len(X2), "F2_essais": essais2}}
    for nom, XL in (("F1", X1), ("F2", X2)):
        out[nom] = Rf.mesurer(_traj(Rf, np.stack(XL))) if XL else None
    return out
