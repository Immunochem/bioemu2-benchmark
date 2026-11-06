#!/usr/bin/env python3
"""redock_cond.py — les CONDITIONS du redocking, leur verification, et la mesure du temoin negatif. Importe redock_lib (leurres, mesures). Meme code CPU et GPU :
le job GPU et la verification CPU appellent `construire_leurres` avec les memes graines, donc les memes x0.

CONDITIONS POSITIVES (par systeme, K repetitions chacune) :
  C0_d0   cristal tel quel (controle plafond)            T_d10 T_d20 T_d30   translation pure, d_COM = native + delta
  R_d10 R_d20 R_d30   translation + rotation aleatoire   RD_d20             rotation + direction aleatoires, d_COM = native + 20
CONDITIONS NEGATIVES (pas d orientation native : T et C0 n existent pas) : R_d10 R_d20 R_d30, RD_d20.
"""
import sys
import numpy as np
import mdtraj as md
from scipy.spatial import cKDTree
import redock_lib as RL
from redock_lib import Q

sys.dont_write_bytecode = True
COND_POS = [("C0", 0.0), ("T", 10.0), ("T", 20.0), ("T", 30.0), ("R", 10.0), ("R", 20.0), ("R", 30.0), ("RD", 20.0)]
COND_NEG = [("R", 10.0), ("R", 20.0), ("R", 30.0), ("RD", 20.0)]
K_DEFAUT = 3


def nom_cond(arm, delta):
    return f"{arm}_d{int(delta)}"


def verifier_leurre(Rf, X):
    """Verifications d un leurre X (A) contre la reference : memes atomes ; chaque chaine INTACTE (distances internes inchangees, liaisons C-N conservees) ;
    distance atomique minimale entre chaines. Rend un dict de constats chiffres ; `ok` = intacts ET fini."""
    ref = Rf.xyz_A
    fini = bool(np.isfinite(X).all()) and X.shape == ref.shape
    out = {"fini": fini}
    if not fini:
        out["ok"] = False
        return out
    ecart = []
    for c in (0, 1):
        i = Rf.idx[c]
        # chaine intacte : superposition rigide sur la reference -> RMSD ~ 0 (la chaine n est que deplacee)
        Rot, Pm, Tm = RL.kabsch(X[i], ref[i])
        ecart.append(float(np.sqrt(np.mean(np.sum(((X[i] - Pm) @ Rot + Tm - ref[i]) ** 2, axis=1)))))
    out["rmsd_interne_chaine_A_A"] = [round(e, 5) for e in ecart]
    # liaisons C(i)-N(i+1) : memes distances que dans le cristal (aucune rupture introduite)
    top = Rf.ref.topology
    rupt = 0
    for ch in top.chains:
        res = list(ch.residues)
        for a, b in zip(res, res[1:]):
            ca = next((x.index for x in a.atoms if x.name == "C"), None); nb = next((x.index for x in b.atoms if x.name == "N"), None)
            if ca is None or nb is None:
                continue
            d1, d0 = np.linalg.norm(X[ca] - X[nb]), np.linalg.norm(ref[ca] - ref[nb])
            if abs(d1 - d0) > 1e-3:
                rupt += 1
    out["liaisons_CN_modifiees"] = rupt
    out["dmin_inter_A"] = round(RL.dmin_inter(Rf, X), 3)
    out["d_com_A"] = round(float(np.linalg.norm(X[Rf.idx[0]].mean(axis=0) - X[Rf.idx[1]].mean(axis=0))), 3)
    out["ok"] = max(ecart) < 1e-3 and rupt == 0
    return out


def construire_leurres(Rf, ident, conditions, K=K_DEFAUT):
    """-> (liste de X en A pour les leurres ACCEPTES, audit[list de dict, un par (condition, repetition), accepte ou non])."""
    frames, audit = [], []
    for arm, delta in conditions:
        for rep in range(K):
            seed = RL.seed_de(ident, arm, delta, rep)
            X, au = RL.leurre(Rf, arm, delta, np.random.default_rng(seed))
            ent = {"cond": nom_cond(arm, delta), "arm": arm, "delta_A": delta, "rep": rep, "seed": seed, **au}
            if X is not None:
                v = verifier_leurre(Rf, X)
                ent["verif"] = v
                ent["accepte"] = bool(au["accepte"] and v["ok"])
                if ent["accepte"]:
                    ent["indice_image"] = len(frames)
                    frames.append(X)
            audit.append(ent)
    return frames, audit


def traj_de(Rf, frames):
    return md.Trajectory(xyz=(np.stack(frames) / 10.0).astype(np.float32), topology=Rf.ref.topology)


# ───────────────────────── temoin negatif ─────────────────────────
def mesurer_neg(Rf, traj, site_chaine, site_residus):
    """Temoin negatif : pas de bonne reponse, donc PAS de fnat. On mesure s il y a une interface (contacts, d_COM, distance min) et OU elle est :
    `site_residus` = residus (ordinaux de la chaine `site_chaine`) qui, dans le cristal d ORIGINE de cette chaine, touchent son vrai partenaire.
    rappel_site = part du site natif touchee par l etrangere ; precision_site = part des residus touches qui sont dans le site ; base = |site| / N_residus (au hasard)."""
    top = traj.topology
    ct = Q.contacts_residus(traj, Q.residus_par_chaine(top))
    S = set(site_residus)
    n_res = len(Rf.res[site_chaine])
    X = traj.xyz * 10.0
    out = {"n_contacts": [], "n_residus_touches_chaine_site": [], "rappel_site": [], "precision_site": [], "d_com_A": Q.d_com_atomes(traj, *Rf.idx).round(3).tolist(), "dmin_inter_A": [],
           "base_hasard_precision": len(S) / n_res, "taille_site": len(S), "n_residus_chaine_site": n_res}
    for f, c in enumerate(ct):
        touches = {p[site_chaine] for p in c}
        out["n_contacts"].append(len(c))
        out["n_residus_touches_chaine_site"].append(len(touches))
        out["rappel_site"].append(round(len(S & touches) / len(S), 4) if S else None)
        out["precision_site"].append(round(len(S & touches) / len(touches), 4) if touches else None)
        out["dmin_inter_A"].append(round(float(cKDTree(X[f][Rf.idx[1]]).query(X[f][Rf.idx[0]])[0].min()), 3))
    out["derive_A_ca_A"] = Q.derive_chaine(traj, Rf.ref, Rf.ca[0]).round(3).tolist()
    out["derive_B_ca_A"] = Q.derive_chaine(traj, Rf.ref, Rf.ca[1]).round(3).tolist()
    return out
