#!/usr/bin/env python3
"""cv_t4l.py — les variables collectives de Bhakat (JCIM 2026, 10.1021/acs.jcim.6c02044) appliquees a une trajectoire mdtraj.

DEFINITIONS (papier, section Methods, eq. 1 ; code de l'auteur Codes/locking_cord.py et Codes/distance.py) :
  d1 = distance Ca Ser44 - Ile150            (nm)   d1 < 2,5 nm = FERME, d1 > 2,5 nm = OUVERT
  d2 = distance Ca Glu22 - Gln141            (nm)
  p  = (d_60,67 . d_67,4) / |d_60,67|^2      d_60,67 = Ca67 - Ca60 ; d_67,4 = centre(CG,CD1,CD2,CE1,CE2,CZ de Phe4) - Ca67
       p > 0 : Phe4 EXPOSE · p < 0 : Phe4 ENFOUI

NUMEROTATION — APPARIEMENT ORDINAL, JAMAIS PAR resSeq.
  Le residu « n » (1-based) est le n-ieme residu de la chaine, quelle que soit la valeur de resSeq dans le fichier.
  Raison mesuree : le topology.pdb de BioEmu v1 est numerote de 0 (MET A 0 ... LEU A 163) ; le `locking_cord.py` de
  l'auteur selectionne « residue N » = resSeq, ce qui l'a oblige a ecrire une variante a part pour BioEmu (res 59/66/3).
  Ici : on compte les residus, puis on VERIFIE leur nom. Un nom inattendu arrete le calcul (ValueError), il ne le corrige pas.

UNITES : mdtraj rend des nm. Aucune conversion ici. Le seuil 2,5 est en nm.

CE QUI N'EST PAS CACHE : p n'est defini que si les 6 atomes du cycle de Phe4 existent. Un topology « squelette + CB »
(BioEmu v1) rend p = NaN et la trajectoire est DECLAREE sans p ; ce n'est jamais un p = 0 ni un « expose » par defaut.
"""
import numpy as np
import mdtraj as md

SEQ_PAPIER = ("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCNGVITK"
              "DEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRCALINMVFQMGETGVAGFTNSLRM"
              "LQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL")
ATTENDU = {4: "PHE", 22: "GLU", 44: "SER", 60: "LYS", 67: "PHE", 141: "GLN", 150: "ILE"}
CYCLE = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
SEUIL_D1_NM = 2.5
UN = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
      "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
      "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H", "MSE": "M"}


def residus(top):
    """Liste ordonnee des residus de la SEULE chaine. Plusieurs chaines = arret."""
    if top.n_chains != 1:
        raise ValueError(f"{top.n_chains} chaines ; T4L attendu en une seule chaine")
    return list(list(top.chains)[0].residues)


def sequence(top):
    return "".join(UN.get(r.name, "X") for r in residus(top))


def indices(top):
    """Indices d'atomes par ordinal (1-based) -> verifies par nom de residu. Rend un dict + un rapport."""
    res = residus(top)
    if len(res) not in (162, 164):      # 164 = sequence du papier ; 162 = residus MODELISES de 5LZM (Asn163, Leu164 : REMARK 465)
        raise ValueError(f"{len(res)} residus, 162 ou 164 attendus")
    for n, nom in ATTENDU.items():
        if res[n - 1].name != nom:
            raise ValueError(f"ordinal {n} = {res[n-1].name} (resSeq {res[n-1].resSeq}), {nom} attendu : numerotation non transferable")
    ca = {}
    for n in (22, 44, 60, 67, 141, 150):
        a = [at.index for at in res[n - 1].atoms if at.name == "CA"]
        if len(a) != 1:
            raise ValueError(f"CA du residu ordinal {n} : {len(a)} atome(s)")
        ca[n] = a[0]
    cyc = [at.index for at in res[3].atoms if at.name in CYCLE]
    cb = [at.index for at in res[3].atoms if at.name == "CB"]
    noms = sorted(at.name for at in res[3].atoms if at.name in CYCLE)
    rapport = {"resSeq_ordinaux": {n: res[n - 1].resSeq for n in ATTENDU}, "atomes_cycle_phe4": noms,
               "p_defini": len(cyc) == 6, "cb_phe4": len(cb) == 1, "n_residus": len(res),
               "diffs_sequence_vs_papier": [(i + 1, SEQ_PAPIER[i], s) for i, s in enumerate(sequence(top)) if s != SEQ_PAPIER[i]]}
    rapport["_cb"] = cb[0] if len(cb) == 1 else None
    return ca, (cyc if len(cyc) == 6 else None), rapport


def cv(traj):
    """-> dict d1, d2 (nm), p (NaN si le cycle manque), et le rapport d'indices."""
    ca, cyc, rap = indices(traj.topology)
    xyz = traj.xyz.astype(np.float64)
    d1 = np.linalg.norm(xyz[:, ca[44]] - xyz[:, ca[150]], axis=1)
    d2 = np.linalg.norm(xyz[:, ca[22]] - xyz[:, ca[141]], axis=1)
    if cyc is None:
        p = np.full(traj.n_frames, np.nan)
    else:
        a = xyz[:, ca[67]] - xyz[:, ca[60]]
        b = xyz[:, cyc].mean(axis=1) - xyz[:, ca[67]]
        p = (a * b).sum(axis=1) / (a * a).sum(axis=1)
    # PROXY, PAS LA DEFINITION DE L'AUTEUR : meme formule avec le CB de Phe4 a la place du centre du cycle (pour un topology squelette+CB comme BioEmu v1).
    cb = rap.pop("_cb")
    if cb is None:
        p_cb = np.full(traj.n_frames, np.nan)
    else:
        a = xyz[:, ca[67]] - xyz[:, ca[60]]
        b = xyz[:, cb] - xyz[:, ca[67]]
        p_cb = (a * b).sum(axis=1) / (a * a).sum(axis=1)
    return {"d1": d1, "d2": d2, "p": p, "p_cb_proxy": p_cb, "rapport": rap}


def etat(d1, p):
    """0 exposed/open · 1 exposed/closed · 2 buried/open · 3 buried/closed · -1 INDEFINI (p NaN ou p == 0).
    Les seuils sont les siens : d1 < 2,5 nm ferme ; p > 0 expose ; p < 0 enfoui. p == 0 n'est dans AUCUN etat chez lui
    (masques lock > 0 et lock < 0) : ici il est compte a part, jamais range en silence."""
    d1 = np.asarray(d1); p = np.asarray(p)
    out = np.full(d1.shape, -1, dtype=int)
    ferme = d1 < SEUIL_D1_NM
    ok = np.isfinite(p)
    out[ok & (p > 0) & ~ferme] = 0
    out[ok & (p > 0) & ferme] = 1
    out[ok & (p < 0) & ~ferme] = 2
    out[ok & (p < 0) & ferme] = 3
    return out


NOMS = ["exposed/open", "exposed/closed", "buried/open", "buried/closed"]


def populations(d1, p):
    e = etat(d1, p)
    n = len(e)
    c = [int((e == k).sum()) for k in range(4)]
    return {"n_total": n, "n_indefini": int((e == -1).sum()), "comptes": dict(zip(NOMS, c)),
            "pourcent_sur_definis": {k: (100.0 * v / max(1, sum(c))) for k, v in zip(NOMS, c)}}
