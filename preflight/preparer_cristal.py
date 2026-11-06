#!/usr/bin/env python3
"""preparer_cristal.py — un cristal AVIDbase -> ce qu il faut pour conditionner BioEmu2.

═══ POURQUOI CE SCRIPT EXISTE ═══
BioEmu2 genere depuis des SEQUENCES. Le cristal, lui, a des residus manquants
(densite absente) et une numerotation qui saute. Si on donne a BioEmu2 la sequence
complete pendant que DynaPIN compare au cristal, tout ce qui suit un gap glisse d un
cran et le DockQ est faux SANS AVERTISSEMENT. C est le piege « un index n est pas
transferable », deja paye sur ce projet.

⇒ On extrait la sequence OBSERVEE, on refuse les cas ou la numerotation n est pas
   exploitable, et on ecrit l index de correspondance pour que la comparaison soit
   verifiable apres coup.

═══ QU EST-CE QU UN RESIDU DE LA CHAINE ? (corrige le 2026-11-04) ═══
Ni « ATOM vs HETATM », ni le NOM du residu : la CONNECTIVITE PEPTIDIQUE, c est-a-dire
la distance C(i)–N(i+1) entre residus consecutifs dans le fichier, dans la plage
covalente [LIEN_CN_MIN_NM, LIEN_CN_MAX_NM].
  - un residu d acide amine (standard OU modifie : MSE, SEP, TPO, PTR, MLY, CSO...) lie
    en peptide a un voisin APPARTIENT a la chaine ;
  - un acide amine LIBRE (aucune liaison C–N avec un voisin ; additif de cristallisation,
    souvent un HETATM portant un OXT) N APPARTIENT PAS a la chaine, quel que soit son
    numero ;
  - un gap (« rupture ») est une rupture de cette liaison entre deux membres consecutifs.
    Un saut de numerotation resSeq SANS rupture (numerotation IMGT, auteur) n est PAS un gap.
Deux defauts mesures disparaissent : `MSE` (HETATM en ligne dans la chaine) etait exclu par
le filtre sur le nom et creait un faux gap ; un `SER 501` libre etait compte comme residu
et creait un faux gap de 84 residus. Le nom ne sert plus qu a choisir la LETTRE (table
canonique wwPDB/CCD de Biopython), jamais a decider de l appartenance.

═══ CONVENTIONS, ET LEUR SOURCE ═══
  contact d interface : 5,0 A entre atomes LOURDS  -> convention CAPRI/DockQ
                        (PRD-etage-MD.md §3.1, item 2 : « 5,0 A obligatoire pour tout
                        f_nat ou DockQ publie, sinon le chiffre n est comparable a rien »)
  borne de contact BioEmu2 : distance Ca + 1,0 A arrondie au 0,5 A superieur
                        -> recette de examples/sample_protein_protein_constraints.py
  🔴 Le conditionnement BioEmu2 BIAISE, il ne contraint pas (idem, verbatim :
     « Contacts bias predictions; they do not enforce hard limits »).
  perimetre : <= 1030 residus observes (cible + binder). Exposant VRAM BioEmu2 MESURE
     0,98 (227 tokens -> 2,55 Gio ; 720 -> 7,92), soit ~11,3 Gio a 1030 tokens, dans une
     slice 1g.24gb de 23,62 Gio. 🔴 Cet exposant est celui de BioEmu2 : ne jamais y
     transporter le 2,96 de RF3 (une contrainte appartient a son moteur).

Usage : preparer_cristal.py <cristal.pdb> <chaine_cible> <chaine_binder> [--sortie X.json]
"""
import argparse
import json
import math
import sys

import mdtraj as md
import numpy as np
from Bio import __version__ as BIOPYTHON_VERSION
# Table canonique nom de composant -> lettre du parent standard, GENEREE par Biopython
# depuis le dictionnaire de composants wwPDB (« component->parent relationships »).
# Elle ne decide JAMAIS de l appartenance a la chaine, seulement de la lettre.
from Bio.Data.PDBData import protein_letters_3to1_extended as TABLE_CCD

# 5,0 A en nanometres — mdtraj travaille en nm, et une unite non declaree est une
# erreur qui se lit comme un resultat.
SEUIL_CONTACT_NM = 0.50
MARGE_CA_NM = 0.10  # +1,0 A sur la distance Ca, recette BioEmu2
PALIER_NM = 0.05    # arrondi au 0,5 A superieur

# Liaison peptidique C(i)–N(i+1), en nm.
# 🔴 PROVENANCE DU DEFAUT LIEN_CN_MAX_NM = 0,29 nm (2,9 A) — mesure du 2026-11-04,
#    prep_fix/resultats/socle_creux.out (tests/socle_creux.py), jamais postulee :
#  POPULATION A — 12 structures MD du socle DynaBench (`<PDB>_input.pdb`), 24 chaines, 6 731
#    paires consecutives : 6 726 paires entre 1,143 et 1,6 A ; 4 jonctions de boucles
#    RECONSTRUITES a 2,014 / 2,115 / 2,403 / 2,550 A (API RCSB : positions d entite adjacentes,
#    aucun residu manquant entre a et b) ; 1 vraie coupure a 12,480 A (1ZM4, 18 residus d entite
#    manquants entre les deux, selon l API). Creux de A : (2,550 ; 12,480), rapport 4,9.
#  POPULATION B — 79 cristaux RCSB (X-ray), 162 chaines, 49 335 paires : liaisons 1,227–1,444 A,
#    AUCUNE paire entre 1,5 et 3,3 A, 117 ruptures de 3,301 a 66,8 A (105 confirmees par l API :
#    residus d entite manquants entre les deux ; 12 ou l API reste muette, C–N de 4,7 a 24,4 A).
#  Defaut = milieu de (2,550 ; 3,301) = 2,93 -> 2,9 A : au-dessus de la plus longue jonction
#    continue de A, au-dessous de la plus courte rupture de B.
#  🔴 ET NON le milieu du creux de A seul (~7,5 A) : il lierait 24 des 105 ruptures confirmees de
#    B (d <= 7,55 A), c est-a-dire de vrais gaps. Le creux de A n a qu UNE coupure de ce cote.
#  LIMITES : 4 jonctions d un cote, 1 rupture minimale de l autre ; marges 0,35 A (2,9 - 2,55) et
#    0,40 A (3,301 - 2,9). Calibre sur MD DynaBench + X-ray ; NON verifie sur cryo-EM ni sur une
#    autre preparation MD. A l ancien defaut 1,9 A, 3 des 12 complexes du socle (1BKD, 2AJF,
#    2OT3) sortaient avec des gaps qu un residu reconstruit n est pas.
LIEN_CN_MIN_NM = 0.10
LIEN_CN_MAX_NM = 0.29   # 2,9 A : population A (6 731 paires MD) + population B (49 335 paires cristaux), ci-dessus

# Perimetre : residus observes, cible + binder. Voir l en-tete pour la mesure VRAM.
MAX_RESIDUS_DEFAUT = 1030

# Les 20 standard : sert UNIQUEMENT a dire quels residus de la chaine sont « modifies ».
TROIS_VERS_UN = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def lettre_de(nom):
    """Lettre du residu : le parent standard (table CCD), sinon `X`.
    `X` = residu LIE en peptide mais sans parent standard connu de la table (UNK, SEC,
    PYL, acides amines D...). Il est compte et SIGNALE, jamais ignore."""
    return TABLE_CCD.get(nom, TROIS_VERS_UN.get(nom, "X"))


def dorsale(r):
    """{'N': idx, 'CA': idx, 'C': idx} (index d atomes) ou None si un des trois manque."""
    d = {}
    for at in r.atoms:
        if at.name in ("N", "CA", "C") and at.name not in d:
            d[at.name] = at.index
    return d if len(d) == 3 else None


def distance_CN_nm(xyz, dorsale_i, dorsale_j):
    """Distance C(i) -> N(j) en nm. LA mesure de connectivite peptidique."""
    return float(np.linalg.norm(xyz[dorsale_i["C"]] - xyz[dorsale_j["N"]]))


def est_lie(d_nm):
    """Vrai si la distance C–N est dans la plage covalente (None = pas de voisin)."""
    return d_nm is not None and LIEN_CN_MIN_NM <= d_nm <= LIEN_CN_MAX_NM


def residus_de(traj, lettre_chaine):
    """Les residus de la CHAINE de lettre donnee, dans l ordre du fichier.

    Appartient a la chaine tout residu d acide amine (N, CA, C presents) LIE en peptide a
    un voisin. Retourne un dict :
      residus        liste de residus mdtraj (les membres), ordre du fichier
      d_suivant_nm   distance C(i)->N(i+1) vers le membre suivant (None pour le dernier)
      exclus         residus d acide amine LIBRES (non lies), hors de l etendue de la chaine
    """
    top = traj.topology
    xyz = traj.xyz[0]
    cibles = [c for c in top.chains if chaine_id(top, c) == lettre_chaine]
    if not cibles:
        dispo = sorted({chaine_id(top, c) for c in top.chains})
        sys.exit(f"🔴 ARRET : chaine « {lettre_chaine} » absente. Presentes : {dispo}")

    # candidats : tout residu qui a une dorsale N-CA-C (standard, modifie ou libre)
    candidats, sans_dorsale = [], []
    for ch in cibles:
        for r in ch.residues:
            d = dorsale(r)
            if d is not None:
                candidats.append((r, d))
            elif r.name in TABLE_CCD or r.name in TROIS_VERS_UN or r.name == "UNK":
                sans_dorsale.append(r)
    if sans_dorsale:
        liste = ", ".join(f"{r.name} {r.resSeq}" for r in sans_dorsale[:8])
        sys.exit(f"🔴 ARRET : chaine « {lettre_chaine} » : {len(sans_dorsale)} residu(s) "
                 f"d acide amine SANS dorsale N-CA-C complete ({liste}). La liaison "
                 f"peptidique ne peut pas etre mesuree : les ignorer decalerait l ordinal "
                 f"en silence.")

    # liaison C(i)–N(i+1) entre candidats consecutifs, ordre du fichier
    n = len(candidats)
    d_suiv = [distance_CN_nm(xyz, candidats[i][1], candidats[i + 1][1])
              for i in range(n - 1)] + [None]
    lie = [est_lie(d_suiv[i]) or (i > 0 and est_lie(d_suiv[i - 1])) for i in range(n)]
    idx_lies = [i for i in range(n) if lie[i]]
    if not idx_lies:
        sys.exit(f"🔴 ARRET : chaine « {lettre_chaine} » sans AUCUN residu lie en peptide "
                 f"(C–N dans [{LIEN_CN_MIN_NM * 10:.1f} ; {LIEN_CN_MAX_NM * 10:.1f}] A) : "
                 f"{n} acide(s) amine(s) libre(s), pas une chaine.")
    premier, dernier = idx_lies[0], idx_lies[-1]

    membres, exclus, internes = [], [], []
    for i, (r, d) in enumerate(candidats):
        if lie[i]:
            membres.append(i)
        elif premier < i < dernier:
            internes.append(r)
        else:
            exclus.append({
                "residu": r.name, "resSeq": int(r.resSeq),
                "n_atomes": r.n_atoms,
                "OXT": any(at.name == "OXT" for at in r.atoms),
                "position": "apres la chaine" if i > dernier else "avant la chaine",
                "d_C_N_voisin_angstrom": [
                    None if i == 0 else round(distance_CN_nm(
                        xyz, candidats[i - 1][1], d) * 10, 3),
                    None if i == n - 1 else round(d_suiv[i] * 10, 3)],
                "raison": "dorsale N-CA-C mais AUCUNE liaison peptidique C–N avec un voisin : "
                          "acide amine libre (ou ligand), pas un residu de la chaine",
            })
    if internes:
        liste = ", ".join(f"{r.name} {r.resSeq}" for r in internes[:8])
        sys.exit(f"🔴 ARRET : chaine « {lettre_chaine} » : {len(internes)} residu(s) "
                 f"d acide amine ISOLE(S) a l INTERIEUR de la chaine ({liste}), sans liaison "
                 f"peptidique avec leurs voisins. Ambigu : residu de la chaine entre deux "
                 f"gaps, ou acide amine libre ? Ni l un ni l autre n est tranche en silence.")

    # d_suivant entre MEMBRES (consecutifs dans la liste des candidats : aucun interne)
    out_res = [candidats[i][0] for i in membres]
    out_d = [d_suiv[i] for i in membres]
    out_d[-1] = None
    return {"residus": out_res, "d_suivant_nm": out_d, "exclus": exclus}


def chaine_id(top, chaine):
    """mdtraj perd les identifiants de chaine du PDB selon la version : on prend
    `chain_id` quand il existe, sinon on retombe sur l index alphabetique — et on le
    DIT, parce qu un repli silencieux sur A/B est exactement comment on compare deux
    chaines differentes sans le savoir."""
    cid = getattr(chaine, "chain_id", None)
    if cid:
        return cid
    return chr(ord("A") + chaine.index)


def codes_insertion(chemin_pdb):
    """Les codes d insertion lus dans le PDB BRUT, par chaine.

    🔴 mdtraj ECRASE le code d insertion : 111, 111A et 111B deviennent trois residus
    de `resSeq` 111. Les VHH sont numerotes IMGT et le CDR3 porte justement des
    insertions (111.x / 112.x) — donc `resSeq` n est PAS une identite exploitable ici.
    On lit le fichier pour pouvoir le DIRE, et la correspondance cristal <-> echantillon
    BioEmu2 se fait sur la position ORDINALE, seule grandeur que les deux partagent.
    """
    par_chaine = {}
    with open(chemin_pdb) as f:
        for ligne in f:
            if not ligne.startswith(("ATOM  ", "HETATM")):
                continue
            code = ligne[26:27]
            if code != " ":
                par_chaine.setdefault(ligne[21:22], set()).add(
                    f"{ligne[22:26].strip()}{code}")
    return {c: sorted(v) for c, v in par_chaine.items()}


def garde_numerotation(residus, d_suivant_nm, etiquette, insertions):
    """Decrit la numerotation et refuse ce qui rendrait la correspondance ambigue.

    Ce qui est TOLERE et documente : les doublons de `resSeq` dus aux codes
    d insertion, et les gaps de densite. Ce qui est REFUSE : une numerotation qui
    DECROIT, parce qu alors l ordre du fichier ne suit pas la sequence et l ordinal
    cesse d etre une identite.

    UN GAP est une RUPTURE de la liaison peptidique C–N entre deux membres consecutifs
    (`d_suivant_nm`), PAS un saut de `resSeq` : un saut de numerotation avec liaison intacte
    (numerotation IMGT, numerotation d auteur) est compte a part et n est pas un gap.
    """
    numeros = [r.resSeq for r in residus]
    if any(b < a for a, b in zip(numeros, numeros[1:])):
        sys.exit(f"🔴 ARRET ({etiquette}) : numerotation DECROISSANTE. L ordre du "
                 f"fichier ne suit pas la sequence, l ordinal n est plus une identite.")
    doublons = len(numeros) - len(set(numeros))
    ruptures = [k for k in range(len(residus) - 1) if not est_lie(d_suivant_nm[k])]
    sauts_sans_rupture = [k for k in range(len(residus) - 1)
                          if numeros[k + 1] > numeros[k] + 1 and est_lie(d_suivant_nm[k])]
    liens = [d_suivant_nm[k] * 10 for k in range(len(residus) - 1)
             if est_lie(d_suivant_nm[k])]
    return {
        "n_residus_observes": len(residus),
        "identite": "position ORDINALE 1..N dans la sequence observee",
        "premier_resSeq": numeros[0],
        "dernier_resSeq": numeros[-1],
        "n_resSeq_en_double": doublons,
        "codes_insertion": insertions,
        "n_gaps": len(ruptures),
        "gaps": [{"apres": numeros[k], "avant": numeros[k + 1],
                  "apres_ordinal": k + 1, "avant_ordinal": k + 2,
                  "manquants": max(numeros[k + 1] - numeros[k] - 1, 0),
                  "d_C_N_angstrom": round(d_suivant_nm[k] * 10, 3)} for k in ruptures],
        "n_sauts_resSeq_sans_rupture": len(sauts_sans_rupture),
        "liaison_peptidique_mesuree_angstrom": {
            "n_liaisons": len(liens),
            "min": round(min(liens), 3) if liens else None,
            "max": round(max(liens), 3) if liens else None,
        },
    }


def arrondi_sup(x_nm):
    return math.ceil(x_nm / PALIER_NM) * PALIER_NM


def main():
    global LIEN_CN_MAX_NM
    ap = argparse.ArgumentParser()
    ap.add_argument("cristal")
    ap.add_argument("chaine_cible", help="la chaine antigene (colonne `cible` du manifeste)")
    ap.add_argument("chaine_binder", help="la chaine VHH (colonne `binder` du manifeste)")
    ap.add_argument("--sortie", default=None)
    ap.add_argument("--refuser-gaps", action="store_true",
                    help="refuse un cristal a trous de densite. A utiliser pour "
                         "constituer le lot de benchmark : un gap oblige a donner une "
                         "chaine amputee a BioEmu2, variable qu on ne veut pas melanger "
                         "a ce qu on mesure.")
    ap.add_argument("--max-contacts", type=int, default=3,
                    help="nombre de contacts de conditionnement. 3 = la recette de "
                         "l exemple officiel. Plus de contacts = plus de biais, et le "
                         "biais est ce que le benchmark doit mesurer, pas presupposer.")
    ap.add_argument("--max-residus", type=int, default=MAX_RESIDUS_DEFAUT,
                    help=f"perimetre : refuse un complexe de plus de N residus observes "
                         f"(cible + binder). Defaut {MAX_RESIDUS_DEFAUT} = limite MESUREE "
                         f"de BioEmu2 dans une slice 1g.24gb (exposant VRAM 0,98).")
    ap.add_argument("--lien-cn-max-angstrom", type=float, default=LIEN_CN_MAX_NM * 10,
                    help=f"borne haute de la liaison peptidique C(i)–N(i+1), en A (defaut "
                         f"{LIEN_CN_MAX_NM * 10:.1f}, calibre sur 12 structures MD du socle + 79 "
                         f"cristaux : voir le commentaire de LIEN_CN_MAX_NM). A 1.9, les jonctions "
                         f"de boucles reconstruites (2,0–2,55 A) comptent comme gaps.")
    a = ap.parse_args()
    LIEN_CN_MAX_NM = a.lien_cn_max_angstrom / 10     # nm : une unite non declaree est une erreur

    traj = md.load(a.cristal)
    if traj.n_frames != 1:
        sys.exit(f"🔴 ARRET : {traj.n_frames} modeles dans le cristal, attendu 1.")

    chaine_c = residus_de(traj, a.chaine_cible)
    chaine_b = residus_de(traj, a.chaine_binder)
    res_cible, res_binder = chaine_c["residus"], chaine_b["residus"]

    total = len(res_cible) + len(res_binder)
    if total > a.max_residus:
        sys.exit(f"🔴 ARRET : {total} residus observes (cible {len(res_cible)} + binder "
                 f"{len(res_binder)}) > perimetre {a.max_residus}. Au-dela, la VRAM "
                 f"BioEmu2 (exposant 0,98 mesure) n est plus garantie dans une slice "
                 f"1g.24gb. --max-residus pour passer outre, en connaissance de cause.")

    ins = codes_insertion(a.cristal)
    # ── l ordinal, 1-indexe, est LA cle de correspondance avec un echantillon BioEmu2
    ordinal = {}
    for lot in (res_cible, res_binder):
        for i, r in enumerate(lot, start=1):
            ordinal[r.index] = i

    info = {
        "cristal": a.cristal,
        "chaine_cible": a.chaine_cible,
        "chaine_binder": a.chaine_binder,
        "seuil_contact_angstrom": SEUIL_CONTACT_NM * 10,
        "convention_contact": "CAPRI/DockQ, atomes lourds, 5,0 A (PRD §3.1 item 2)",
        "appartenance_chaine": {
            "critere": "connectivite peptidique C(i)-N(i+1) entre residus consecutifs",
            "fenetre_angstrom": [LIEN_CN_MIN_NM * 10, LIEN_CN_MAX_NM * 10],
            "table_nom_vers_lettre": f"Bio.Data.PDBData.protein_letters_3to1_extended "
                                     f"(biopython {BIOPYTHON_VERSION}, {len(TABLE_CCD)} entrees)",
        },
        "perimetre": {"max_residus": a.max_residus, "total_observe": total},
        "cible": garde_numerotation(res_cible, chaine_c["d_suivant_nm"],
                                    f"cible {a.chaine_cible}", ins.get(a.chaine_cible, [])),
        "binder": garde_numerotation(res_binder, chaine_b["d_suivant_nm"],
                                     f"binder {a.chaine_binder}", ins.get(a.chaine_binder, [])),
    }
    if a.refuser_gaps and (info["cible"]["n_gaps"] or info["binder"]["n_gaps"]):
        sys.exit(f"🔴 ARRET : gaps de densite "
                 f"(cible {info['cible']['n_gaps']}, binder {info['binder']['n_gaps']}) "
                 f"et --refuser-gaps demande. Un gap oblige a donner a BioEmu2 une "
                 f"chaine amputee : c est une variable de plus dans le benchmark.")
    for role, lot, chaine in (("cible", res_cible, chaine_c), ("binder", res_binder, chaine_b)):
        info[role]["sequence_observee"] = "".join(lettre_de(r.name) for r in lot)
        info[role]["acides_amines_modifies"] = [
            {"residu": r.name, "resSeq": int(r.resSeq), "ordinal": i,
             "lettre": lettre_de(r.name)}
            for i, r in enumerate(lot, start=1) if r.name not in TROIS_VERS_UN]
        info[role]["residus_exclus_non_lies"] = chaine["exclus"]

    # ── contacts d interface, atomes LOURDS uniquement
    xyz = traj.xyz[0]
    lourds = lambda r: [at.index for at in r.atoms if at.element.symbol != "H"]
    paires = []
    for rc in res_cible:
        ic = lourds(rc)
        if not ic:
            continue
        pc = xyz[ic]
        for rb in res_binder:
            ib = lourds(rb)
            if not ib:
                continue
            d = np.linalg.norm(pc[:, None, :] - xyz[ib][None, :, :], axis=-1)
            dmin = float(d.min())
            if dmin <= SEUIL_CONTACT_NM:
                paires.append((dmin, rc, rb))

    if not paires:
        sys.exit(f"🔴 ARRET : AUCUN contact a {SEUIL_CONTACT_NM * 10:.1f} A entre "
                 f"« {a.chaine_cible} » et « {a.chaine_binder} ». Ce PDB n est pas un "
                 f"complexe, ou les chaines sont mal nommees.")

    info["n_paires_interface"] = len(paires)
    # en ORDINAL, puisque c est la seule cle que le cristal et un echantillon partagent
    info["residus_interface_cible"] = sorted({ordinal[rc.index] for _, rc, _ in paires})
    info["residus_interface_binder"] = sorted({ordinal[rb.index] for _, _, rb in paires})

    # ── conditionnement : les contacts les plus serres, bornes sur la distance Ca
    ca = {}
    for r in res_cible + res_binder:
        sel = [at.index for at in r.atoms if at.name == "CA"]
        if sel:
            ca[r.index] = xyz[sel[0]]

    contacts = []
    for dmin, rc, rb in sorted(paires, key=lambda p: p[0]):
        if rc.index not in ca or rb.index not in ca:
            continue  # pas de Ca -> inutilisable pour une borne Ca, on passe
        d_ca = float(np.linalg.norm(ca[rc.index] - ca[rb.index]))
        contacts.append({
            "cible_ordinal": ordinal[rc.index],
            "cible_resSeq": int(rc.resSeq),
            "cible_residu": rc.name,
            "binder_ordinal": ordinal[rb.index],
            "binder_resSeq": int(rb.resSeq),
            "binder_residu": rb.name,
            "d_min_lourds_angstrom": round(dmin * 10, 3),
            "d_ca_angstrom": round(d_ca * 10, 3),
            "borne_angstrom": round(arrondi_sup(d_ca + MARGE_CA_NM) * 10, 1),
        })
        if len(contacts) >= a.max_contacts:
            break
    info["contacts_conditionnement"] = contacts

    texte = json.dumps(info, indent=2, ensure_ascii=False)
    if a.sortie:
        with open(a.sortie, "w") as f:
            f.write(texte + "\n")
        print(f"ecrit : {a.sortie}")
    else:
        print(texte)

    # ── le resume qui se lit sans ouvrir le JSON
    for role, lettre in (("cible ", a.chaine_cible), ("binder", a.chaine_binder)):
        d = info[role.strip()]
        print(f"{role} {lettre} : {d['n_residus_observes']} residus observes "
              f"· {d['n_gaps']} gap(s) · resSeq {d['premier_resSeq']}-{d['dernier_resSeq']}"
              + (f" · {d['n_resSeq_en_double']} resSeq en double "
                 f"(insertions {','.join(d['codes_insertion'])})"
                 if d["n_resSeq_en_double"] else "")
              + (f" · {len(d['acides_amines_modifies'])} modifie(s) lie(s) en peptide "
                 f"({','.join(sorted({m['residu'] for m in d['acides_amines_modifies']}))})"
                 if d["acides_amines_modifies"] else "")
              + (f" · {len(d['residus_exclus_non_lies'])} acide(s) amine(s) LIBRE(S) exclu(s) "
                 f"({','.join(e['residu'] + str(e['resSeq']) for e in d['residus_exclus_non_lies'])})"
                 if d["residus_exclus_non_lies"] else ""), file=sys.stderr)
    print(f"interface : {info['n_paires_interface']} paires a 5,0 A · "
          f"{len(info['residus_interface_cible'])} residus cible · "
          f"{len(info['residus_interface_binder'])} residus binder", file=sys.stderr)
    if info["cible"]["n_gaps"] or info["binder"]["n_gaps"]:
        print("🔸 Des gaps existent : la sequence donnee a BioEmu2 doit etre la sequence "
              "OBSERVEE ci-dessus, pas la sequence complete de la PDB.", file=sys.stderr)
    if info["cible"]["n_resSeq_en_double"] or info["binder"]["n_resSeq_en_double"]:
        print("🔸 Codes d insertion presents (numerotation IMGT du CDR3). La "
              "correspondance avec un echantillon BioEmu2 se fait sur l ORDINAL : "
              "`resSeq` seul designerait plusieurs residus.", file=sys.stderr)
    n_x = sum(info[r]["sequence_observee"].count("X") for r in ("cible", "binder"))
    if n_x:
        print(f"🔸 {n_x} residu(s) LIE(S) en peptide sans parent standard dans la table "
              f"CCD : lettre `X` dans la sequence observee (voir `acides_amines_modifies`). "
              f"A trancher avant de donner cette sequence a BioEmu2.", file=sys.stderr)


if __name__ == "__main__":
    main()
