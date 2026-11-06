#!/usr/bin/env python3
"""curer_entree.py — une structure de MD DynaBench -> un x0 que BioEmu2 accepte, ou un REFUS bruyant.

═══ POURQUOI ═══
Le contrat d'entree de BioEmu2 se valide ENTIEREMENT sur CPU : `System.from_description` et
`map_cif_atom_order` ne chargent aucun modele et n'ont pas besoin de GPU. Trois jobs GPU ont
ete perdus un par un sur des blocages que cette porte voit en 0,2 s : noms d'ATOMES CHARMM
(ILE `CD` au lieu de `CD1`), nom de RESIDU non standard (`DDE`), identifiants de CHAINE non
sequentiels (`A`,`C`). Ce script LIT d'abord tout ce qui peut bloquer, CLASSE chaque blocage,
corrige UNIQUEMENT ce qui preserve l'objet (renommer, renumeroter), REFUSE ce qui le changerait.

═══ CRITERE : LA COMMENSURABILITE ═══
L'objet donne a BioEmu2 et l'objet que la MD a simule portent LES MEMES RESIDUS.
  - renommer un atome (`CD`->`CD1`), un residu de protonation (`HSD`->`HIS`), une chaine
    (`C`->`B`), renumeroter : preservent l'objet -> applique, DECLARE, par table EXPLICITE ;
  - retirer l'adduit d'un `DDE`, combler une rupture de chaine, supprimer un residu, effacer un
    ion ou un ligand : changent l'objet -> REFUS, avec sa classe et le message litteral.
Un renommage HORS de la table explicite fait refuser le systeme, meme si mdtraj le fait en
silence : ce que mdtraj convertit implicitement est COMPARE a la table, atome par atome.

═══ CE QUI EST IMPORTE, PAS REECRIT ═══
  prep_fix/preparer_cristal.py (sha256 epingle ci-dessous) : residus_de, garde_numerotation,
      lettre_de, dorsale, chaine_id, TABLE_CCD, TROIS_VERS_UN  (connectivite C–N, ruptures)
  6CWG/scripts/preparer.py (sha256 epingle ci-dessous) : cles, el  (cles d'atomes lourds)
  bioemu2 : System.from_description, map_cif_atom_order  (la porte ; CPU ; aucun modele)
La liste des atomes attendus par residu est LUE dans la topologie de BioEmu2 elle-meme, jamais
recopiee.

Usage : curer_entree.py <structure.pdb> <x0_sortie.pdb> [--rapport r.json]
Code de sortie : 0 = PRET (la porte passe sur le x0 ecrit) ; 3 = REFUSE (classes + messages) ;
                 2 = usage.
"""
import argparse
import collections
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True          # ne jamais ecrire un __pycache__ dans 6CWG/ (autre lot)
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # le contrat se valide sur CPU ; aucun GPU, jamais

import mdtraj as md
import numpy as np
from scipy.spatial import cKDTree

# DECLARED DEVIATION (public copy, 2026-11-06): the two pinned dependencies are published
# next to this file instead of at their cluster paths. The pinned sha256 are UNCHANGED and
# still verified at load time, so the functions this script imports are the ones that ran.
RACINE = Path(__file__).resolve().parent
PREP_CRISTAL = RACINE / "preparer_cristal.py"
SHA_PREP_CRISTAL = "8ead19bb668b4c952c6c455de21db95d78e808ed9e2e48df8387d176707f4924"
PREP_6CWG = RACINE / "preparer_6cwg.py"
SHA_PREP_6CWG = "f3650e47d4cfa9c6bd4a9279441228bb2453e15a22a7e563d06a050c5633542d"


def charger(chemin, nom, sha_attendu):
    sha = hashlib.sha256(Path(chemin).read_bytes()).hexdigest()
    if sha != sha_attendu:
        sys.exit(f"🔴 ARRET : {chemin} a change (sha256 {sha[:16]}…, epingle {sha_attendu[:16]}…). "
                 f"La curation importe ces fonctions : elle ne suit pas un changement en silence.")
    spec = importlib.util.spec_from_file_location(nom, chemin)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


prep = charger(PREP_CRISTAL, "preparer_cristal", SHA_PREP_CRISTAL)
prep6 = charger(PREP_6CWG, "preparer_6cwg", SHA_PREP_6CWG)

# ══════════ LA TABLE EXPLICITE ══════════
# Renommage de RESIDUS : etats de protonation de l'histidine. Memes atomes lourds, meme residu.
RENOMMAGE_RESIDUS = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS"}
# Renommage d'ATOMES LOURDS, convention CHARMM -> convention PDB (celle de BioEmu2) :
RENOMMAGE_ATOMES_PAR_RESIDU = {("ILE", "CD"): "CD1"}        # CHARMM nomme CD le delta de l'isoleucine
RENOMMAGE_ATOMES_TERMINAL_C = {"OT1": "O", "OT2": "OXT"}    # oxygenes du carboxyle C-terminal CHARMM
# Atome ECARTE par contrat (absent du modele de BioEmu2, ignore explicitement par l'amont,
# `bioemu2/structure/ref_struct.py:172` ignorable_atoms=("OXT", "OP3", "NT")) :
ATOMES_ECARTES = ("OXT",)
SEUIL_SUPERPOSITION_A = 0.5     # deux atomes LOURDS a moins de 0,5 A ne sont pas des atomes distincts


def contrat_standard():
    """{nom_residu: [atomes lourds]} lu dans la topologie de BioEmu2 : LE contrat, pas une copie."""
    from bioemu2 import System
    top = System.from_description(proteins=["ACDEFGHIKLMNPQRSTVWY"]).topology
    return {r.name: [a.name for a in r.atoms] for r in top.residues}


# ══════════ LECTURE BRUTE : ce que mdtraj et biotite masquent ══════════
def scanner_brut(chemin):
    """Colonnes fixes du PDB : altLoc, code d'insertion, occupation, element, TER, MODEL, chaines."""
    s = {"n_atom": 0, "n_hetatm": 0, "n_ter": 0, "n_model": 0, "n_altloc": 0, "n_icode": 0,
         "n_occ_lt1": 0, "n_occ_0": 0, "n_element_vide": 0, "n_nonfini": 0}
    runs, residus, cle_prec = [], [], None
    for l in open(chemin):
        r = l[:6]
        if r == "TER   ":
            s["n_ter"] += 1
        elif r == "MODEL ":
            s["n_model"] += 1
        if r not in ("ATOM  ", "HETATM"):
            continue
        s["n_atom" if r == "ATOM  " else "n_hetatm"] += 1
        chain, seq, ic, nom_res = l[21], l[22:26].strip(), l[26], l[17:20].strip()
        if l[16] != " ":
            s["n_altloc"] += 1
        if ic != " ":
            s["n_icode"] += 1
        try:
            occ = float(l[54:60])
            if occ < 1.0:
                s["n_occ_lt1"] += 1
            if occ == 0.0:
                s["n_occ_0"] += 1
            xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            if not all(np.isfinite(xyz)):
                s["n_nonfini"] += 1
        except ValueError:
            s["n_nonfini"] += 1
        if not l[76:78].strip():
            s["n_element_vide"] += 1
        if not runs or runs[-1] != chain:
            runs.append(chain)
        cle = (chain, seq, ic, nom_res)
        if cle != cle_prec:
            residus.append(cle)
            cle_prec = cle
    s["chaines_runs"] = runs
    s["chaines_discontigues"] = sorted({c for c in runs if runs.count(c) > 1})
    s["n_residus_brut"] = len(residus)
    # numerotation en ordre du fichier, par chaine
    par = collections.OrderedDict()
    for ch, seq, ic, nm in residus:
        par.setdefault(ch, []).append((int(seq) if seq.lstrip("-").isdigit() else None, ic, nm))
    s["decroissances"] = {ch: sum(1 for a, b in zip(v, v[1:]) if a[0] is not None and b[0] is not None and b[0] < a[0])
                          for ch, v in par.items()}
    s["resseq_doubles_consecutifs"] = {ch: sum(1 for a, b in zip(v, v[1:]) if a[0] == b[0]) for ch, v in par.items()}
    s["resseq_doubles_non_consecutifs"] = {ch: len(v) - len({(x[0], x[1]) for x in v}) for ch, v in par.items()}
    s["residus_par_chaine"] = {ch: len(v) for ch, v in par.items()}
    s["premier_dernier_resseq"] = {ch: [v[0][0], v[-1][0]] for ch, v in par.items()}
    return s, par


# ══════════ LA PORTE (CPU) ══════════
def porte_cpu(sequences, chemin_x0):
    """System.from_description puis map_cif_atom_order. Ni modele ni GPU. Retourne les messages LITTERAUX."""
    from bioemu2 import System, map_cif_atom_order
    out = {"from_description": None, "map_cif_atom_order": None}
    try:
        system = System.from_description(proteins=list(sequences))
        out["from_description"] = "OK"
        out["chain_ids_systeme"] = list(system.chain_ids)
    except Exception as e:
        out["from_description"] = f"{type(e).__name__}: {e}"[:700]
        return out, None
    try:
        x = map_cif_atom_order(system, str(chemin_x0))
        out["map_cif_atom_order"] = "OK"
        out["n_atomes_mappes"] = int(x.n_atoms)
    except Exception as e:
        out["map_cif_atom_order"] = f"{type(e).__name__}: {e}"[:700]
        return out, None
    # PORTE ETENDUE : les controles d'entree de `propagate_with_embeddings`
    # (bioemu2/propagation.py:401 `topology, batch = prepare_sampling_inputs(...)`, :405 topologie du x0
    # == topologie du systeme, :410 coordonnees finies), executes sur CPU AVANT tout passage du modele.
    try:
        import torch
        from bioemu2.structure.preparation import prepare_sampling_inputs
        topo, _batch = prepare_sampling_inputs(system, device=torch.device("cpu"))
        out["prepare_sampling_inputs"] = "OK"
        out["x0_topologie_egale_au_systeme"] = bool(x.topology == topo)
        out["x0_coordonnees_finies"] = bool(np.isfinite(x.xyz).all())
    except Exception as e:
        out["prepare_sampling_inputs"] = f"{type(e).__name__}: {e}"[:700]
        return out, None
    if not (out["x0_topologie_egale_au_systeme"] and out["x0_coordonnees_finies"]):
        return out, None
    return out, x


def porte_sur_brut(chemin):
    """TEMOIN NEGATIF : la porte sur la structure BRUTE, avec les sequences que donnerait
    `topology.to_fasta()` (comme le pilote 6CWG, propager.py:23). Retourne les messages litteraux."""
    t = md.load(str(chemin))
    return porte_cpu(list(t.topology.to_fasta()), chemin)[0]


# ══════════ LA CURATION ══════════
def curer(chemin, sortie_x0=None, attendu_corpus=None, ecrire_diag=None, contrat=None, xtc=None, image=0):
    """Enveloppe : `_curer` ne doit jamais faire echapper une exception ; si c'est le cas, c'est UN REFUS inscrit
    (classe ERREUR_INATTENDUE, message litteral), pas un arret du lot."""
    try:
        return _curer(chemin, sortie_x0, attendu_corpus, ecrire_diag, contrat, xtc, image)
    except SystemExit:
        raise
    except Exception as e:
        return {"entree": str(chemin), "etat": "REFUSE", "classes_refus": ["ERREUR_INATTENDUE"],
                "constats": [{"classe": "ERREUR_INATTENDUE", "gravite": "REFUS", "n": 1, "detail": "exception non prevue dans la curation",
                              "message": f"{type(e).__name__}: {e}"[:500], "exemples": None}]}


def _curer(chemin, sortie_x0=None, attendu_corpus=None, ecrire_diag=None, contrat=None, xtc=None, image=0):
    """Retourne un dict ; ne leve pas pour un refus : un refus est une DONNEE.
    xtc / image : si fournis, les COORDONNEES sont celles de l'image `image` du xtc (noms et topologie du PDB) :
    c'est l'« image 0 de la replique » du pilote, qui n'est PAS `structure.pdb` pour 36 repliques sur 84.
    attendu_corpus : {chaine_originale: {"n_residus":..., "resseq_premier":..., "resseq_dernier":...}} (manifeste)."""
    contrat = contrat or contrat_standard()
    R = {"entree": str(chemin), "sha256_entree": hashlib.sha256(Path(chemin).read_bytes()).hexdigest(),
         "constats": [], "etat": None}

    def noter(classe, gravite, n, detail, message=None, exemples=None):
        R["constats"].append({"classe": classe, "gravite": gravite, "n": n, "detail": detail,
                              "message": message, "exemples": exemples})

    scan, par_chaine_brut = scanner_brut(chemin)
    R["scan_brut"] = scan
    # ── 1. format
    if scan["n_model"] > 1:
        noter("FORMAT_MODELES", "REFUS", scan["n_model"], f"{scan['n_model']} MODEL dans le fichier (une seule structure attendue)")
    if scan["chaines_discontigues"]:
        noter("CHAINES_DISCONTIGUES", "REFUS", len(scan["chaines_discontigues"]),
              f"la chaine {scan['chaines_discontigues']} reapparait apres une autre : biotite leve « Duplicate or discontiguous chain »")
    if scan["n_altloc"]:
        noter("ALTLOC", "REFUS", scan["n_altloc"], "positions alternatives (altLoc) : atomes dupliques, objet ambigu")
    if scan["n_occ_0"]:
        noter("OCCUPATION_NULLE", "REFUS", scan["n_occ_0"], "atomes a occupation 0,00 (coordonnees non observees)")
    if scan["n_occ_lt1"] - scan["n_occ_0"] > 0:
        noter("OCCUPATION_PARTIELLE", "INFO", scan["n_occ_lt1"] - scan["n_occ_0"], "atomes a occupation partielle")
    if scan["n_nonfini"]:
        noter("COORDONNEES", "REFUS", scan["n_nonfini"], "coordonnees non finies ou illisibles")
    if scan["n_element_vide"]:
        noter("ELEMENT_VIDE", "INFO", scan["n_element_vide"], "colonne element vide : elements devines par mdtraj, ecrits explicitement dans le x0")
    if scan["n_icode"]:
        noter("CODES_INSERTION", "RENUMEROTE", scan["n_icode"], "codes d'insertion : resSeq ambigu ; renumerotation sequentielle par ordre du fichier")
    if scan["n_ter"] == 0:
        noter("TER_ABSENT", "INFO", 0, "aucun TER : les chaines sont separees par l'identifiant de chaine seul (suffisant pour biotite)")
    if any(scan["decroissances"].values()):
        noter("NUMEROTATION_DECROISSANTE", "RENUMEROTE", sum(scan["decroissances"].values()),
              f"decroissances par chaine {scan['decroissances']} ; renumerotation sequentielle si les liaisons C–N sont intactes")
    if any(scan["resseq_doubles_consecutifs"].values()) and not scan["n_icode"]:
        noter("RESSEQ_DOUBLES", "RENUMEROTE", sum(scan["resseq_doubles_consecutifs"].values()), "residus consecutifs de meme resSeq")
    if any(scan["resseq_doubles_non_consecutifs"].values()):
        noter("RESSEQ_REPETES", "RENUMEROTE", sum(scan["resseq_doubles_non_consecutifs"].values()),
              "resSeq repetes dans une chaine : biotite leve « Duplicate or discontiguous residue » ; renumerotation sequentielle")

    # ── 2. lecture mdtraj : BRUT (noms tels quels) et STANDARD (ce que mdtraj convertit implicitement)
    try:
        brut = md.load(str(chemin), standard_names=False)
        std = md.load(str(chemin))
    except Exception as e:
        noter("FORMAT_ILLISIBLE", "REFUS", 1, f"mdtraj : {type(e).__name__}: {e}"[:300])
        R["etat"] = "REFUSE"
        return R
    if brut.n_frames != 1 or std.n_frames != 1:
        noter("FORMAT_MODELES", "REFUS", brut.n_frames, f"{brut.n_frames} images")
    if xtc is not None:
        try:
            frame = md.load_frame(str(xtc), image, top=str(chemin))
        except Exception as e:
            noter("FORMAT_ILLISIBLE", "REFUS", 1, f"xtc image {image} : {type(e).__name__}: {e}"[:300])
            R["etat"] = "REFUSE"
            return R
        if frame.n_atoms != std.n_atoms:
            noter("ECART_LECTURE", "REFUS", 1, f"xtc {frame.n_atoms} atomes, PDB {std.n_atoms}")
            R["etat"] = "REFUSE"
            return R
        std = md.Trajectory(xyz=frame.xyz, topology=std.topology)     # coordonnees de l'image, noms du PDB
        R["coordonnees"] = f"{xtc} image {image}"
    else:
        R["coordonnees"] = "PDB"
    if brut.n_atoms != std.n_atoms or brut.n_residues != std.n_residues:
        noter("ECART_MDTRAJ", "REFUS", 1, "mdtraj lit un nombre d'atomes/residus different avec et sans noms standard")
        R["etat"] = "REFUSE"
        return R
    tb, ts = brut.topology, std.topology
    n_res_brut = tb.n_residues
    if n_res_brut != scan["n_residus_brut"]:
        noter("ECART_LECTURE", "REFUS", 1, f"lecture brute {scan['n_residus_brut']} residus, mdtraj {n_res_brut} : decoupage en residus different")

    # ── 3. correspondance AVANT -> APRES (mdtraj implicite), par atome
    # cles d'atomes lourds : fonction du pilote 6CWG (preparer.py:21 `cles`), importee, pas reecrite
    try:
        kb, ks = prep6.cles(tb), prep6.cles(ts)
    except SystemExit as e:
        noter("ELEMENT_INCONNU", "REFUS", 1, str(e)[:300])
        R["etat"] = "REFUSE"
        return R
    R["cles_lourdes_6cwg"] = {"avant": len(kb), "apres": len(ks), "seulement_avant": len(set(kb) - set(ks)), "seulement_apres": len(set(ks) - set(kb))}
    obs = collections.Counter()
    obsH = collections.Counter()
    for ab, as_ in zip(tb.atoms, ts.atoms):
        cle = (ab.residue.name, as_.residue.name, ab.name, as_.name)
        if as_.element is not None and as_.element.symbol == "H":
            obsH[cle] += 1
        else:
            obs[cle] += 1
    R["correspondance_mdtraj_lourds"] = [[a, b, c, d, n] for (a, b, c, d), n in sorted(obs.items())]
    R["correspondance_mdtraj_H"] = [[a, b, c, d, n] for (a, b, c, d), n in sorted(obsH.items())]

    # ── 4. chaines proteiques, dans l'ordre du fichier
    lettres = []
    for c in ts.chains:
        L = prep.chaine_id(ts, c)
        if L not in lettres:
            lettres.append(L)
    prot, non_prot = [], []
    for L in lettres:
        cand = [r for c in ts.chains if prep.chaine_id(ts, c) == L for r in c.residues if prep.dorsale(r) is not None]
        (prot if cand else non_prot).append(L)
    if non_prot:
        noms = sorted({r.name for c in ts.chains if prep.chaine_id(ts, c) in non_prot for r in c.residues})
        noter("HORS_PROTEINE", "REFUS", len(non_prot), f"chaines sans residu d'acide amine {non_prot} (residus {noms}) : ions/ligands/solvant ; les effacer changerait l'objet")
    # residus non proteiques DANS une chaine proteique
    for L in prot:
        autres = [r for c in ts.chains if prep.chaine_id(ts, c) == L for r in c.residues if prep.dorsale(r) is None]
        if autres:
            noter("HORS_PROTEINE", "REFUS", len(autres), f"chaine {L} : residus sans dorsale N-CA-C {sorted({r.name for r in autres})}", exemples=[f"{r.name}{r.resSeq}" for r in autres[:6]])

    # ── 5. connectivite, ruptures, numerotation : fonctions de preparer_cristal (importees)
    chaines = []     # [{lettre, residus (mdtraj std), d_suivant}]
    for L in prot:
        try:
            res = prep.residus_de(std, L)
        except SystemExit as e:
            noter("RESIDU_ISOLE_OU_SANS_DORSALE", "REFUS", 1, f"chaine {L}", message=str(e)[:400])
            continue
        n_tot = sum(c.n_residues for c in ts.chains if prep.chaine_id(ts, c) == L)
        if res["exclus"]:
            noter("ACIDE_AMINE_LIBRE", "REFUS", len(res["exclus"]), f"chaine {L} : acide amine non lie en peptide {[x['residu'] + str(x['resSeq']) for x in res['exclus']]}")
        rs, d = res["residus"], res["d_suivant_nm"]
        ruptures = [(k, d[k]) for k in range(len(rs) - 1) if not prep.est_lie(d[k])]
        if ruptures:
            noter("RUPTURE_CHAINE", "REFUS", len(ruptures),
                  f"chaine {L} : {len(ruptures)} rupture(s) C–N (> {prep.LIEN_CN_MAX_NM * 10:.1f} A) ; combler changerait l'objet",
                  exemples=[f"{rs[k].name}{rs[k].resSeq}->{rs[k + 1].name}{rs[k + 1].resSeq} {dd * 10:.2f} A" for k, dd in ruptures[:6]])
        try:
            g = prep.garde_numerotation(rs, d, f"chaine {L}", [])
        except SystemExit as e:
            g = None
            if not any(c["classe"] == "NUMEROTATION_DECROISSANTE" for c in R["constats"]):
                noter("NUMEROTATION_DECROISSANTE", "RENUMEROTE", 1, f"garde_numerotation : {str(e)[:200]}")
        chaines.append({"lettre": L, "residus": rs, "d_suivant_nm": d, "n_tot": n_tot, "garde": g})
        if g is not None and g["n_sauts_resSeq_sans_rupture"]:
            noter("SAUTS_RESSEQ_LIAISON_INTACTE", "RENUMEROTE", g["n_sauts_resSeq_sans_rupture"], f"chaine {L} : sauts de numerotation a liaison C–N intacte (pas des gaps)")

    # ── 6. coherence avec le manifeste du corpus (l'objet donne = l'objet du corpus)
    if attendu_corpus:
        for ch in chaines:
            a = attendu_corpus.get(ch["lettre"])
            if a is None:
                noter("COHERENCE_CORPUS", "REFUS", 1, f"chaine {ch['lettre']} absente du manifeste ({sorted(attendu_corpus)})")
            elif a["n_residus"] != len(ch["residus"]):
                noter("COHERENCE_CORPUS", "REFUS", 1, f"chaine {ch['lettre']} : {len(ch['residus'])} residus, manifeste {a['n_residus']}")
        for L in attendu_corpus:
            if L not in [c["lettre"] for c in chaines]:
                noter("COHERENCE_CORPUS", "REFUS", 1, f"chaine {L} du manifeste absente de la structure")

    # ── 7. application de la table EXPLICITE (residus puis atomes), et comparaison a mdtraj
    ecarts_mdtraj, noms_finals, nouv_res = collections.Counter(), {}, {}
    derniers = {}
    for ch in chaines:
        derniers[ch["residus"][-1].index] = ch["lettre"]
    renomme_res, renomme_at = collections.Counter(), collections.Counter()
    for r_b, r_s in zip(tb.residues, ts.residues):
        nom_res = RENOMMAGE_RESIDUS.get(r_b.name, r_b.name)
        nouv_res[r_s.index] = nom_res
        if nom_res != r_b.name:
            renomme_res[(r_b.name, nom_res)] += 1
        if nom_res != r_s.name:
            ecarts_mdtraj[("RESIDU", r_b.name, nom_res, r_s.name)] += 1
        noms_atomes = [a.name for a in r_b.atoms]
        terminal = r_b.index in derniers
        for a_b, a_s in zip(r_b.atoms, r_s.atoms):
            if a_s.element is not None and a_s.element.symbol == "H":
                continue
            nm = a_b.name
            if (nom_res, nm) in RENOMMAGE_ATOMES_PAR_RESIDU:
                nm2 = RENOMMAGE_ATOMES_PAR_RESIDU[(nom_res, nm)]
            elif terminal and nm in RENOMMAGE_ATOMES_TERMINAL_C and RENOMMAGE_ATOMES_TERMINAL_C[nm] not in noms_atomes:
                nm2 = RENOMMAGE_ATOMES_TERMINAL_C[nm]
            else:
                nm2 = nm
            noms_finals[a_b.index] = nm2
            if nm2 != nm:
                renomme_at[(nom_res, nm, nm2)] += 1
            if nm2 != a_s.name:
                ecarts_mdtraj[("ATOME", nom_res, nm, f"table={nm2} | mdtraj={a_s.name}")] += 1
    R["histidines_nom_orig"] = dict(collections.Counter(r_b.name for r_b in tb.residues if nouv_res[r_b.index] == "HIS"))
    R["renommages_residus_appliques"] = [[a, b, n] for (a, b), n in sorted(renomme_res.items())]
    R["renommages_atomes_appliques"] = [[a, b, c, n] for (a, b, c), n in sorted(renomme_at.items())]
    for k, n in sorted(ecarts_mdtraj.items()):
        pass
    if ecarts_mdtraj:
        noter("ECART_MDTRAJ", "REFUS", sum(ecarts_mdtraj.values()),
              "mdtraj convertit (ou ne convertit pas) autrement que la table explicite : renommage hors table",
              exemples=[f"{t} {a} {b} {c} x{n}" for (t, a, b, c), n in sorted(ecarts_mdtraj.items())[:12]])
    if renomme_res:
        noter("RENOMMAGE_RESIDU", "RENOMME", sum(renomme_res.values()), f"{dict((f'{a}->{b}', n) for (a, b), n in renomme_res.items())}")
    if renomme_at:
        noter("RENOMMAGE_ATOME", "RENOMME", sum(renomme_at.values()), f"{dict((f'{r}:{a}->{b}', n) for (r, a, b), n in renomme_at.items())}")

    # ── 8. conformite au contrat : residus standard, atomes attendus / en trop / manquants / dupliques
    non_std, extra, manq, dup, oxt = collections.Counter(), collections.Counter(), collections.Counter(), collections.Counter(), 0
    ex_non_std = {}
    for ch in chaines:
        for r in ch["residus"]:
            nm = nouv_res[r.index]
            atomes = [noms_finals[a.index] for a in r.atoms if a.index in noms_finals]
            if nm not in contrat:
                parent3 = None
                lettre = prep.lettre_de(nm)
                for k3, v1 in prep.TROIS_VERS_UN.items():
                    if v1 == lettre:
                        parent3 = k3
                n_par = len(contrat.get(parent3, [])) if parent3 else None
                non_std[nm] += 1
                ex_non_std.setdefault(nm, f"{nm}{r.resSeq} chaine {ch['lettre']} : {len(atomes)} atomes lourds"
                                          + (f" contre {n_par} pour {parent3} (parent CCD {lettre})" if n_par else " (aucun parent standard)"))
                continue
            comptes = collections.Counter(atomes)
            for a_nm, c in comptes.items():
                if c > 1:
                    dup[(nm, a_nm)] += 1
            attendus = set(contrat[nm])
            presents = set(atomes)
            for a_nm in sorted(presents - attendus):
                if a_nm in ATOMES_ECARTES:
                    oxt += 1
                else:
                    extra[(nm, a_nm)] += 1
            for a_nm in sorted(attendus - presents):
                manq[(nm, a_nm)] += 1
    if non_std:
        noter("RESIDU_NON_STANDARD", "REFUS", sum(non_std.values()), f"residus hors des 20 standard : {dict(non_std)} ; renommer changerait l'objet",
              exemples=list(ex_non_std.values()))
    if extra:
        noter("ATOMES_NON_CONFORMES", "REFUS", sum(extra.values()), "atomes lourds dont le nom n'est pas dans le contrat de BioEmu2 apres renommage",
              exemples=[f"{r}:{a} x{n}" for (r, a), n in sorted(extra.items())[:12]])
    if manq:
        noter("ATOMES_MANQUANTS", "REFUS", sum(manq.values()), "atomes du contrat absents du x0",
              exemples=[f"{r}:{a} x{n}" for (r, a), n in sorted(manq.items())[:12]])
    if dup:
        noter("ATOMES_DUPLIQUES", "REFUS", sum(dup.values()), "meme nom d'atome deux fois dans un residu",
              exemples=[f"{r}:{a} x{n}" for (r, a), n in sorted(dup.items())[:12]])
    if oxt:
        noter("OXT", "ECARTE", oxt, "OXT absent du contrat de BioEmu2 : ecarte du x0 (l'amont l'ignore deja, ref_struct.py:172)")
    n_H = int(sum(obsH.values()))
    if n_H:
        noter("HYDROGENES", "ECARTE", n_H, "hydrogenes ecartes du x0 : `ignore_h` de l'amont les ignore ; coordonnees des atomes lourds inchangees")

    # ── 9. atomes superposes (atomes lourds distincts a < 0,5 A) et coordonnees
    if not np.isfinite(std.xyz).all():
        noter("COORDONNEES", "REFUS", int((~np.isfinite(std.xyz).all(axis=-1)).sum()), "coordonnees non finies (nan/inf) : superposition, contacts et porte impossibles")
        R["classes_refus"] = sorted({c["classe"] for c in R["constats"] if c["gravite"] == "REFUS"})
        R["etat"] = "REFUSE"
        return R
    garde_idx, chaine_de = [], []
    for i_ch, ch in enumerate(chaines):
        for r in ch["residus"]:
            for a in r.atoms:
                if a.index in noms_finals and noms_finals[a.index] not in ATOMES_ECARTES:
                    garde_idx.append(a.index)
                    chaine_de.append(i_ch)
    chaine_de = np.array(chaine_de)
    xyz = std.xyz[0] * 10.0
    if garde_idx and len(chaines) >= 2:
        arbre5 = cKDTree(xyz[garde_idx])
        croises = [(i, j) for i, j in arbre5.query_pairs(5.0) if chaine_de[i] != chaine_de[j]]
        R["contacts_atomes_inter_chaines_5A"] = len(croises)
        if not croises:
            noter("COMPLEXE_DISLOQUE", "REFUS", 0, "aucun contact atome lourd-atome lourd < 5 A entre chaines : complexe disloque (image periodique ? dissocie ?)")
    elif len(chaines) < 2:
        noter("UNE_SEULE_CHAINE", "INFO", len(chaines), "moins de deux chaines proteiques : pas de complexe")
    if garde_idx:
        arbre = cKDTree(xyz[garde_idx])
        paires = arbre.query_pairs(SEUIL_SUPERPOSITION_A)
        if paires:
            exs = []
            for i, j in sorted(paires)[:5]:
                ai, aj = ts.atom(garde_idx[i]), ts.atom(garde_idx[j])
                exs.append(f"{ai.residue.name}{ai.residue.resSeq}:{ai.name} / {aj.residue.name}{aj.residue.resSeq}:{aj.name}")
            noter("ATOMES_SUPERPOSES", "REFUS", len(paires), f"paires d'atomes lourds a < {SEUIL_SUPERPOSITION_A} A", exemples=exs)

    # ── 10. identifiants de chaine : sequentiels A, B, ... dans l'ordre du fichier
    nouveaux = [chr(ord("A") + i) for i in range(len(chaines))]
    table_chaines = {ch["lettre"]: nouveaux[i] for i, ch in enumerate(chaines)}
    if any(k != v for k, v in table_chaines.items()):
        noter("IDENTIFIANTS_CHAINE", "RENOMME", sum(1 for k, v in table_chaines.items() if k != v),
              f"identifiants {list(table_chaines)} -> {list(table_chaines.values())} (BioEmu2 numerote A, B, ... dans l'ordre des sequences)")
    R["table_chaines"] = table_chaines
    sequences = ["".join(prep.lettre_de(nouv_res[r.index]) for r in ch["residus"]) for ch in chaines]
    R["longueurs_chaines"] = [len(s) for s in sequences]

    # ── 11. le x0 : lourds, noms explicites, chaines A.., residus 1..N par chaine
    refus_avant_porte = [c for c in R["constats"] if c["gravite"] == "REFUS"]
    x0_txt, serie, map_residus = [], 0, []
    for i, ch in enumerate(chaines):
        for o, r in enumerate(ch["residus"], start=1):
            map_residus.append({"chaine_orig": ch["lettre"], "resSeq_orig": int(r.resSeq), "nom_orig": tb.residue(r.index).name,
                                "chaine": nouveaux[i], "ordinal": o, "nom": nouv_res[r.index]})
            for a in r.atoms:
                if a.index not in noms_finals or noms_finals[a.index] in ATOMES_ECARTES:
                    continue
                serie += 1
                nm = noms_finals[a.index]
                nm_pdb = nm if len(nm) == 4 else " " + nm.ljust(3)
                el = ts.atom(a.index).element.symbol
                x, y, z = xyz[a.index]
                x0_txt.append("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (
                    serie % 100000, nm_pdb, nouv_res[r.index][:3], nouveaux[i], o, x, y, z, 1.0, 0.0, el.rjust(2)))
        x0_txt.append("TER\n")
    x0_txt.append("END\n")
    R["n_atomes_x0"] = serie
    R["map_residus"] = map_residus
    cible = None
    if sortie_x0 is not None and not refus_avant_porte:
        cible = Path(sortie_x0)
    elif ecrire_diag is not None and chaines:
        cible = Path(ecrire_diag)       # diagnostic NON LIVRABLE : sert uniquement a lire le message litteral de la porte
        R["x0_diagnostic_non_livrable"] = str(cible)
    if cible is not None and chaines:
        cible.parent.mkdir(parents=True, exist_ok=True)
        cible.write_text("".join(x0_txt))
        R["x0"] = str(cible)
        R["sha256_x0"] = hashlib.sha256(cible.read_bytes()).hexdigest()
        lettres_ok = all(set(s) <= set("ACDEFGHIKLMNPQRSTVWY") for s in sequences)
        if lettres_ok or True:
            gate, mapped = porte_cpu(sequences, cible)
            R["porte"] = gate
            if mapped is not None:
                # l'objet est preserve : les coordonnees mappees sont une permutation des coordonnees du x0
                a = np.round(np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in x0_txt if l.startswith("ATOM")]), 3)
                b = np.round(mapped.xyz[0] * 10.0, 3)
                ordre = lambda m: m[np.lexsort((m[:, 2], m[:, 1], m[:, 0]))]
                R["coordonnees_mappees_egales_au_x0"] = bool(a.shape == b.shape and np.allclose(ordre(a), ordre(b), atol=2e-3))
                if not R["coordonnees_mappees_egales_au_x0"]:
                    noter("PORTE_COORDONNEES", "REFUS", 1, "les coordonnees mappees par BioEmu2 ne sont pas celles du x0 ecrit")
            else:
                if gate["from_description"] != "OK":
                    noter("PORTE_FROM_DESCRIPTION", "REFUS", 1, "System.from_description a leve", message=gate["from_description"])
                elif gate["map_cif_atom_order"] != "OK":
                    noter("PORTE_MAP_CIF_ATOM_ORDER", "REFUS", 1, "map_cif_atom_order a leve", message=gate["map_cif_atom_order"])
                else:
                    noter("PORTE_ETENDUE", "REFUS", 1, "controles d'entree de propagate_with_embeddings (CPU)", message=json.dumps({k: v for k, v in gate.items() if k not in ("chain_ids_systeme",)})[:600])
    refus = [c for c in R["constats"] if c["gravite"] == "REFUS"]
    R["classes_refus"] = sorted({c["classe"] for c in refus})
    porte_ok = (bool(R.get("porte")) and R["porte"]["map_cif_atom_order"] == "OK" and R["porte"].get("prepare_sampling_inputs") == "OK"
                and R["porte"].get("x0_topologie_egale_au_systeme") is True and R["porte"].get("x0_coordonnees_finies") is True
                and R.get("coordonnees_mappees_egales_au_x0") is True)
    R["etat"] = "PRET" if (not refus and porte_ok and cible is not None and sortie_x0 is not None) else "REFUSE"
    return R


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("structure")
    ap.add_argument("x0")
    ap.add_argument("--rapport", default=None)
    ap.add_argument("--xtc", default=None, help="coordonnees de l'image --image de ce xtc (topologie et noms : le PDB)")
    ap.add_argument("--image", type=int, default=0)
    a = ap.parse_args()
    R = curer(a.structure, a.x0, ecrire_diag=str(Path(a.x0).with_suffix(".diagnostic.pdb")), xtc=a.xtc, image=a.image)
    if a.rapport:
        Path(a.rapport).write_text(json.dumps(R, indent=1, ensure_ascii=False, default=str) + "\n")
    print(f"{a.structure} -> {R['etat']}")
    for c in R["constats"]:
        print(f"  [{c['gravite']:<10}] {c['classe']:<30} n={c['n']} {c['detail']}" + (f" | {c['message']}" if c["message"] else "")
              + (f" | ex: {c['exemples'][:4]}" if c["exemples"] else ""))
    sys.exit(0 if R["etat"] == "PRET" else 3)


if __name__ == "__main__":
    main()
