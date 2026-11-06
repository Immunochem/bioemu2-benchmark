#!/usr/bin/env python3
"""tests_curation.py — chaque classe de REFUS de curer_entree.py est-elle VIVANTE ?

Les 28 systemes reels n'en declenchent presque aucune (0 refus sur 178 entrees) : c'est exactement
la situation d'un garde « satisfait par l'absence de cas », ou desarme. On FABRIQUE donc, a partir
d'un vrai systeme (6CWG, `structure.pdb`, 390 residus, 2 chaines A/B, hydrogenes, HSD, ILE CD), des
entrees dont on connait la reponse :
  - defauts qui DOIVENT etre refuses (altLoc, rupture, residu non standard, atome manquant ou
    duplique ou superpose, ligand, chaines entrelacees, plusieurs modeles, occupation nulle,
    complexe disloque, acide amine libre, chaine absente du manifeste, coordonnee non finie, renommage
    que mdtraj fait et que la table n'a pas declare) ;
  - defauts qui DOIVENT etre corriges SANS changer l'objet (codes d'insertion, numerotation decroissante,
    resSeq repetes, chaines B->D, protonation HSE/HIE, element vide) : PRET + classe declaree.
Un cas dont le verdict differe de l'attente est imprime KO et fait sortir 1 : c'est une donnee.
"""
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
os.environ["CUDA_VISIBLE_DEVICES"] = ""
RACINE = Path("/gpfs/phherman/bioemu2_bench")
sys.path.insert(0, str(RACINE / "prep_fix"))
import curer_entree as ce

SORTIE = RACINE / "prep_fix" / "curation" / "tests"
SORTIE.mkdir(parents=True, exist_ok=True)
SYS = RACINE / "corpus" / "systemes" / "6CWG" / "rep1" / "structure.pdb"


def residus(lines):
    """[(chaine, resSeq, icode, nom, [indices de lignes])] dans l'ordre du fichier."""
    out, prev = [], None
    for i, l in enumerate(lines):
        if not l.startswith(("ATOM  ", "HETATM")):
            continue
        k = (l[21], int(l[22:26]), l[26], l[17:20].strip())
        if k != prev:
            out.append((k[0], k[1], k[2], k[3], []))
            prev = k
        out[-1][4].append(i)
    return out


def par_chaine(lines, ch):
    return [r for r in residus(lines) if r[0] == ch]


def retirer_residu(lines, ch, o):
    r = par_chaine(lines, ch)[o]
    s = set(r[4])
    return [l for i, l in enumerate(lines) if i not in s]


def champ(l, a, b, v):
    return l[:a] + v + l[b:]


def main():
    base = SYS.read_text().splitlines(keepends=True)
    L = [json.loads(l) for l in open(RACINE / "corpus" / "manifeste_corpus.json") if l.strip()]
    x6 = [x for x in L if x["pdb"] == "6CWG"][0]
    attendu = {c["id"]: {"n_residus": c["n_residus"], "resseq_premier": c["resseq_premier"], "resseq_dernier": c["resseq_dernier"]} for c in x6["chaines"]}
    contrat = ce.contrat_standard()
    A = par_chaine(base, "A")
    B = par_chaine(base, "B")
    cas = []     # (nom, lignes, attendu_corpus, etat_attendu, classes_attendues, gravites_attendues)

    cas.append(("T00_intact", base, attendu, "PRET", [], []))
    # ── DOIT etre REFUSE
    l = list(base)
    for i in A[100][4][:3]:
        l[i] = champ(l[i], 16, 17, "A")
    cas.append(("T01_altloc", l, attendu, "REFUSE", ["ALTLOC"], []))
    cas.append(("T02_rupture_residu_retire", retirer_residu(base, "A", 50), None, "REFUSE", ["RUPTURE_CHAINE"], []))
    cas.append(("T02b_residu_retire_manifeste", retirer_residu(base, "A", 50), attendu, "REFUSE", ["RUPTURE_CHAINE", "COHERENCE_CORPUS"], []))
    l = list(base)
    for i in A[60][4]:
        l[i] = champ(l[i], 17, 20, "DDE")
    cas.append(("T03_residu_non_standard_DDE", l, attendu, "REFUSE", ["RESIDU_NON_STANDARD"], []))
    l = list(base)
    i_cb = [i for i in A[70][4] if l[i][12:16].strip() == "CB"]
    if i_cb:
        l = [x for j, x in enumerate(l) if j != i_cb[0]]
        cas.append(("T04_atome_manquant_CB", l, attendu, "REFUSE", ["ATOMES_MANQUANTS"], []))
    l = list(base)
    i_ca = [i for i in A[80][4] if l[i][12:16].strip() == "CA"][0]
    l.insert(i_ca + 1, l[i_ca])
    cas.append(("T05_atome_duplique_CA", l, attendu, "REFUSE", ["ATOMES_DUPLIQUES"], []))
    l = list(base)
    i_cb = [i for i in A[90][4] if l[i][12:16].strip() == "CB"][0]
    i_ca = [i for i in A[90][4] if l[i][12:16].strip() == "CA"][0]
    x, y, z = float(l[i_ca][30:38]), float(l[i_ca][38:46]), float(l[i_ca][46:54])
    l[i_cb] = l[i_cb][:30] + f"{x + 0.1:8.3f}{y:8.3f}{z:8.3f}" + l[i_cb][54:]
    cas.append(("T06_atomes_superposes", l, attendu, "REFUSE", ["ATOMES_SUPERPOSES"], []))
    l = list(base) + ["HETATM99999 ZN    ZN C 999      10.000  10.000  10.000  1.00  0.00          ZN\n"]
    cas.append(("T07_ion_hors_proteine", l, attendu, "REFUSE", ["HORS_PROTEINE"], []))
    # chaines entrelacees : la derniere residue de A est deplacee apres la chaine B
    l = list(base)
    derniers = set(A[-1][4])
    reste = [x for i, x in enumerate(l) if i not in derniers]
    dep = [l[i] for i in sorted(derniers)]
    idx_fin_B = max(i for i, x in enumerate(reste) if x.startswith("ATOM") and x[21] == "B")
    cas.append(("T08_chaines_discontigues", reste[: idx_fin_B + 1] + dep + reste[idx_fin_B + 1:], attendu, "REFUSE", ["CHAINES_DISCONTIGUES"], []))
    cas.append(("T09_plusieurs_modeles", ["MODEL        1\n"] + base + ["ENDMDL\n", "MODEL        2\n"] + base + ["ENDMDL\n"], attendu, "REFUSE", ["FORMAT_MODELES"], []))
    l = list(base)
    i0 = [i for i in A[10][4]][0]
    l[i0] = champ(l[i0], 54, 60, "  0.00")
    cas.append(("T10_occupation_nulle", l, attendu, "REFUSE", ["OCCUPATION_NULLE"], []))
    # complexe disloque : toute la chaine B translatee de 100 A
    l = list(base)
    for r in B:
        for i in r[4]:
            x = float(l[i][30:38])
            l[i] = l[i][:30] + f"{x + 100:8.3f}" + l[i][38:]
    cas.append(("T11_complexe_disloque", l, attendu, "REFUSE", ["COMPLEXE_DISLOQUE"], []))
    # acide amine libre (6 atomes d'une SER) insere en fin de chaine A, tres loin
    fin_A = max(i for i, x in enumerate(base) if x.startswith("ATOM") and x[21] == "A")
    ser = [("N", "N", 60.0, 60.0, 60.0), ("CA", "C", 61.4, 60.0, 60.0), ("C", "C", 62.0, 61.4, 60.0),
           ("O", "O", 61.4, 62.4, 60.0), ("CB", "C", 61.9, 59.3, 61.2), ("OG", "O", 63.3, 59.3, 61.2)]
    libre = ["ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (90000 + k, " " + a.ljust(3), "SER", "A", 999, x, y, z, 1.0, 0.0, e.rjust(2))
             for k, (a, e, x, y, z) in enumerate(ser)]
    cas.append(("T12_acide_amine_libre", base[: fin_A + 1] + libre + base[fin_A + 1:], None, "REFUSE", ["ACIDE_AMINE_LIBRE"], []))
    cas.append(("T13_chaine_absente_du_manifeste", [x for x in base if not (x.startswith("ATOM") and x[21] == "B")], attendu, "REFUSE", ["COHERENCE_CORPUS"], []))
    l = list(base)
    i0 = [i for i in A[20][4]][0]
    l[i0] = l[i0][:30] + "     nan" + l[i0][38:]
    cas.append(("T14_coordonnee_nan", l, attendu, "REFUSE", [], []))        # classe : COORDONNEES ou FORMAT_ILLISIBLE, voir resultat
    # renommage que mdtraj FAIT et que la table explicite n'a PAS declare : CYS -> CYX
    l = list(base)
    cys = [r for r in A if r[3] == "CYS"]
    for i in cys[0][4]:
        l[i] = champ(l[i], 17, 20, "CYX")
    cas.append(("T15_renommage_hors_table_CYX", l, attendu, "REFUSE", ["ECART_MDTRAJ"], []))
    # OT1 sur un residu NON terminal : mdtraj le convertit en O, la table ne le fait qu'en C-terminal
    l = list(base)
    i_o = [i for i in A[100][4] if l[i][12:16].strip() == "O"][0]
    l[i_o] = champ(l[i_o], 12, 16, " OT1")
    cas.append(("T16_OT1_residu_interne", l, attendu, "REFUSE", ["ECART_MDTRAJ"], []))
    # ── DOIT etre CORRIGE sans changer l'objet : PRET + classe declaree
    l = list(base)
    for i in A[101][4]:
        l[i] = l[i][:22] + f"{A[100][1]:>4}" + "A" + l[i][27:]
    cas.append(("T20_code_insertion", l, attendu, "PRET", [], ["CODES_INSERTION"]))
    l = list(base)
    for k, r in enumerate(B):
        for i in r[4]:
            l[i] = l[i][:22] + f"{2000 - k:>4}" + l[i][26:]
    cas.append(("T21_numerotation_decroissante", l, attendu, "PRET", [], ["NUMEROTATION_DECROISSANTE"]))
    l = list(base)
    for i in A[150][4]:
        l[i] = l[i][:22] + f"{A[10][1]:>4}" + l[i][26:]
    cas.append(("T22_resseq_repete_non_consecutif", l, attendu, "PRET", [], ["RESSEQ_REPETES"]))
    l = list(base)
    for i, x in enumerate(l):
        if x.startswith("ATOM") and x[21] == "B":
            l[i] = champ(x, 21, 22, "D")
    cas.append(("T23_chaine_B_devient_D", l, None, "PRET", [], ["IDENTIFIANTS_CHAINE"]))
    l = list(base)
    his = [r for r in A if r[3] == "HSD"]
    for i in his[0][4]:
        l[i] = champ(l[i], 17, 20, "HSE")
    cas.append(("T24_HSE_protonation", l, attendu, "PRET", [], ["RENOMMAGE_RESIDU"]))
    l = list(base)
    k = 0
    for i, x in enumerate(l):
        if x.startswith("ATOM") and x[76:78].strip() in ("C", "N") and k < 5:
            l[i] = champ(x, 76, 78, "  ")
            k += 1
    cas.append(("T25_element_vide", l, attendu, "PRET", [], ["ELEMENT_VIDE"]))

    ko, res = 0, []
    print(f"{'cas':<36} {'attendu':<7} {'obtenu':<7} classes de REFUS obtenues | constats non refus")
    for nom, lignes, att, etat_att, cl_att, gr_att in cas:
        p = SORTIE / f"{nom}.pdb"
        p.write_text("".join(lignes))
        try:
            R = ce.curer(p, sortie_x0=SORTIE / f"{nom}.x0.pdb", attendu_corpus=att, ecrire_diag=SORTIE / f"{nom}.diag.pdb", contrat=contrat)
            etat, cl = R["etat"], set(R["classes_refus"])
            autres = {c["classe"] for c in R["constats"] if c["gravite"] != "REFUS"}
            msg = [c.get("message") or c["detail"] for c in R["constats"] if c["gravite"] == "REFUS"][:2]
        except SystemExit as e:
            etat, cl, autres, msg = "SORTIE_BRUSQUE", set(), set(), [str(e)[:200]]
        except Exception as e:
            etat, cl, autres, msg = "EXCEPTION", set(), set(), [f"{type(e).__name__}: {e}"[:200]]
        ok = (etat == etat_att) and set(cl_att) <= cl and set(gr_att) <= autres
        ko += (not ok)
        print(f"{nom:<36} {etat_att:<7} {etat:<7} {sorted(cl)} | {sorted(autres & {'CODES_INSERTION','NUMEROTATION_DECROISSANTE','RESSEQ_REPETES','IDENTIFIANTS_CHAINE','RENOMMAGE_RESIDU','ELEMENT_VIDE'})}  {'OK' if ok else 'KO'}")
        if not ok or nom in ("T14_coordonnee_nan",):
            print(f"      -> {msg}")
        res.append({"cas": nom, "attendu": etat_att, "obtenu": etat, "classes_refus": sorted(cl), "autres": sorted(autres), "ok": ok, "message": msg})
    (RACINE / "prep_fix" / "curation" / "tests_curation.json").write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    print(f"\ncas : {len(cas)} ; verdict different de l'attente : {ko}")
    sys.exit(1 if ko else 0)


if __name__ == "__main__":
    main()
