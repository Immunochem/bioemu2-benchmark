#!/usr/bin/env python3
"""resume_systemes.py — un tableau par SYSTEME x par CLASSE de blocage, lu dans audit_entrees.json (aucun recalcul).

Pour chaque systeme (structure.pdb de la replique 1 = le fichier brut ; le x0 « image 0 » partage ses noms) :
tokens, table de chaines (ancien -> nouveau), renommages d'atomes lourds necessaires (ILE CD->CD1, OT1/OT2 terminaux),
residus HSD->HIS, atomes avec code d'insertion, sauts de resSeq a liaison intacte, atomes OXT et H ecartes, premiere erreur
LITTERALE de la porte sur la structure BRUTE, verdict sur les x0 cures (PRET / total), et refus (classe + message).
Ecrit curation/resume_28.tsv. Les temoins negatifs (1ZM4, 1FQ1) sont a part.
"""
import collections
import json
import re
from pathlib import Path

SORTIE = Path("/gpfs/phherman/bioemu2_bench/prep_fix/curation")


def premiere_erreur(p):
    if not p:
        return "?"
    if p.get("map_cif_atom_order") == "OK":
        return "PASSE"
    m = (p.get("from_description") if p.get("from_description") not in (None, "OK") else p.get("map_cif_atom_order")) or "?"
    m = m.split(": ", 1)[-1]
    return re.sub(r"\[np\.str_\('(\w)'\), np\.str_\('(\w)'\)\]", r"['\1', '\2']", m)[:120]


def main():
    R = json.load(open(SORTIE / "audit_entrees.json"))
    par = collections.OrderedDict()
    for r in R:
        par.setdefault((r["lot"], r["pdb"]), []).append(r)
    ent = "systeme\ttokens\tchaines_ancien_nouveau\tlongueurs\tILE_CD_vers_CD1\tOT1_OT2_terminaux\tHSD_vers_HIS\tatomes_code_insertion\tsauts_resSeq_liaison_intacte\tOXT_ecartes\tH_ecartes\tporte_brute_premiere_erreur\tx0_PRET_sur_total\tclasses_de_refus"
    lignes = [ent]
    print(ent.replace("\t", " | "))
    for (lot, pdb), rs in par.items():
        s1 = next(r for r in rs if r["kind"] == "structure" and r["replique"] == 1)
        i1 = next(r for r in rs if r["kind"] == "img0" and r["replique"] == 1)
        ren = collections.Counter()
        for res, nm, nm2, n in s1.get("renommages_atomes_appliques", []):
            ren[(res if nm == "CD" else "terminal", nm, nm2)] += n
        c = {x["classe"]: x for x in s1["constats"]}
        ile = ren.get(("ILE", "CD", "CD1"), 0)
        ot = sum(n for (k, nm, nm2), n in ren.items() if k == "terminal")
        hsd = sum(n for a, b, n in s1.get("renommages_residus_appliques", []) if a == "HSD")
        icode = c.get("CODES_INSERTION", {}).get("n", 0)
        sauts = c.get("SAUTS_RESSEQ_LIAISON_INTACTE", {}).get("n", 0)
        oxt = c.get("OXT", {}).get("n", 0)
        h = c.get("HYDROGENES", {}).get("n", 0)
        tab = s1.get("table_chaines", {})
        tc = ",".join(f"{a}->{b}" for a, b in tab.items())
        tok = sum(i1.get("longueurs_chaines", [])) if i1.get("longueurs_chaines") else "?"
        pret = sum(1 for r in rs if r["etat"] == "PRET")
        refus = collections.Counter(cl for r in rs for cl in r["classes_refus"])
        ligne = "\t".join(map(str, [pdb + ("*" if lot != "perimetre" else ""), tok, tc, "/".join(map(str, s1.get("longueurs_chaines", []))), ile, ot, hsd, icode, sauts, oxt, h,
                                    premiere_erreur(s1.get("porte_sur_brut")), f"{pret}/{len(rs)}", ",".join(f"{k}x{v}" for k, v in refus.items()) or "-"]))
        lignes.append(ligne)
        print(ligne.replace("\t", " | "))
    (SORTIE / "resume_28.tsv").write_text("\n".join(lignes) + "\n")
    msgs = collections.Counter()
    for r in R:
        for cst in r["constats"]:
            if cst["gravite"] == "REFUS":
                msgs[(r["lot"], cst["classe"], (cst["message"] or cst["detail"])[:160])] += 1
    print("\nREJETS (entrees) : classe + message litteral, par lot")
    if not msgs:
        print("  aucun")
    for (lot, cl, m), n in sorted(msgs.items()):
        print(f"  {lot:<15} x{n:<3} {cl:<26} {m}")
    cnt_p = collections.Counter(cl for r in R if r["lot"] == "perimetre" for cl in r["classes_refus"])
    print("\ncompte de rejets par classe, PERIMETRE (178 x0) :", dict(cnt_p) or "aucun")
    cnt_t = collections.Counter(cl for r in R if r["lot"] != "perimetre" for cl in r["classes_refus"])
    print("compte de rejets par classe, TEMOINS NEGATIFS (14 x0) :", dict(cnt_t))
    brute = collections.Counter(premiere_erreur(r.get("porte_sur_brut"))[:60] for r in R if r["lot"] == "perimetre" and r["kind"] != "img0")
    print("\nporte CPU sur les 94 fichiers PDB bruts du perimetre, premiere erreur :", dict(brute))
    print("x0 PRETS du perimetre :", sum(1 for r in R if r["lot"] == "perimetre" and r["etat"] == "PRET"), "/", sum(1 for r in R if r["lot"] == "perimetre"))
    his = collections.Counter()
    for r in R:
        if r["lot"] == "perimetre" and r["kind"] == "img0" and r.get("histidines_nom_orig"):
            his.update(r["histidines_nom_orig"])
    print("histidines par nom d'origine dans les 84 x0 image 0 du perimetre (champ histidines_nom_orig) :", dict(his))


if __name__ == "__main__":
    main()
