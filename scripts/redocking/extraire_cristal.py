#!/usr/bin/env python3
"""extraire_cristal.py — un mmCIF RCSB -> un PDB brut de DEUX chaines proteiques, pret pour `prep_fix/curer_entree.py`.

PREPARATION MINIMALE, SANS RIEN INVENTER : modele 1, positions alternatives (une seule par residu : '.' sinon la plus petite lettre), atomes lourds,
residus proteiques standard SEULEMENT, chaines renommees A, B dans l ordre demande, residus renumerotes 1..N par chaine (ordre label_seq_id).
Ce qui changerait l objet est REFUSE, avec sa classe, par chaine : residu non standard DANS le polymere (MSE, SEP...), saut de label_seq_id (residus non
observes internes), chaine absente. Les eaux, ions et ligands (label_seq_id '.') sont ECARTES et comptes, pas refuses : BioEmu2 ne les recoit pas.
La table (auth_asym_id, auth_seq_id, code d insertion) -> ordinal est ecrite dans le JSON compagnon : tout appariement aval est ORDINAL, jamais par resSeq.

Usage : extraire_cristal.py <entree.cif> <chaineA,chaineB> <sortie_brut.pdb> <sortie.json>   (chaines = auth_asym_id)
Code de sortie : 0 = ecrit ; 3 = REFUSE (le JSON porte les classes) ; 2 = usage.
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np
import biotite.structure.io.pdbx as pdbx

STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
UN = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def extraire(cif, chaines_auth):
    """-> (resultat, refus). resultat['chaines'] = [{auth, residus:[{nom, auth_seq, icode, label_seq, atomes:[(nom, elem, x, y, z, occ, b)]}], refus:[classes]}].
    `refus` = union des refus de toutes les chaines demandees (chaque entree porte la chaine dans son detail)."""
    refus, info = [], {}
    c = pdbx.CIFFile.read(str(cif))
    a = c.block["atom_site"]
    col = lambda k: a[k].as_array()
    num = lambda k: a[k].as_array(float)
    elem, nom_a, alt = col("type_symbol"), col("label_atom_id"), col("label_alt_id")
    nom_r, lseq, ins = col("label_comp_id"), col("label_seq_id"), col("pdbx_PDB_ins_code")
    aseq, achn, mod = col("auth_seq_id"), col("auth_asym_id"), col("pdbx_PDB_model_num")
    X, Y, Z, occ, bf = num("Cartn_x"), num("Cartn_y"), num("Cartn_z"), num("occupancy"), num("B_iso_or_equiv")
    info["n_modeles"] = int(len(set(mod.tolist())))
    ok_mod = mod == mod[0]
    info["n_atomes_modele1"] = int(ok_mod.sum())
    resultat = []
    ecartes = {"eau": 0, "non_polymere_dans_chaine": 0, "hydrogenes": 0, "occupation_nulle": 0}
    for ca in chaines_auth:
        rc = []
        m = ok_mod & (achn == ca)
        if not m.any():
            refus.append({"classe": "CHAINE_ABSENTE", "detail": f"auth_asym_id {ca} absent du modele 1"})
            continue
        poly = m & (lseq != ".")
        ecartes["non_polymere_dans_chaine"] += int((m & (lseq == ".")).sum())
        ecartes["eau"] += int((m & (nom_r == "HOH")).sum())
        par_res = {}
        for i in np.nonzero(poly)[0]:
            par_res.setdefault(int(lseq[i]), []).append(i)
        ls_tri = sorted(par_res)
        residus, non_std = [], []
        for ls in ls_tri:
            ii = par_res[ls]
            noms = {nom_r[i] for i in ii}
            if len(noms) != 1:
                rc.append("RESIDU_HETEROGENE"); refus.append({"classe": "RESIDU_HETEROGENE", "detail": f"{ca} label_seq {ls} : noms {sorted(noms)}"})
                continue
            nm = noms.pop()
            if nm not in STANDARD:
                non_std.append(f"{nm}{aseq[ii[0]]}")
                continue
            lettres = sorted({alt[i] for i in ii if alt[i] != "."})
            garde = lettres[0] if lettres else None
            atomes = []
            for i in ii:
                if alt[i] != "." and alt[i] != garde:
                    continue
                if elem[i] in ("H", "D"):
                    ecartes["hydrogenes"] += 1; continue
                if occ[i] <= 0:
                    ecartes["occupation_nulle"] += 1; continue
                atomes.append((nom_a[i], elem[i], float(X[i]), float(Y[i]), float(Z[i]), float(occ[i]), float(bf[i])))
            residus.append({"nom": nm, "auth_seq": int(aseq[ii[0]]), "icode": ("" if ins[ii[0]] in ("?", ".") else ins[ii[0]]), "label_seq": ls, "atomes": atomes})
        if non_std:
            rc.append("RESIDU_NON_STANDARD"); refus.append({"classe": "RESIDU_NON_STANDARD", "detail": f"chaine {ca} : residus de polymere hors des 20 standard {non_std[:8]} ({len(non_std)})"})
        sauts = [(a_, b_) for a_, b_ in zip(ls_tri, ls_tri[1:]) if b_ - a_ > 1]
        if sauts:
            rc.append("RESIDUS_NON_OBSERVES_INTERNES")
            refus.append({"classe": "RESIDUS_NON_OBSERVES_INTERNES", "detail": f"chaine {ca} : sauts de label_seq_id {sauts[:6]} ({sum(b - a_ - 1 for a_, b in sauts)} residus manquants)"})
        resultat.append({"auth": ca, "residus": residus, "n_residus": len(residus), "refus": rc,
                         "premier_label_seq": ls_tri[0] if ls_tri else None, "dernier_label_seq": ls_tri[-1] if ls_tri else None,
                         "n_manquants_internes": int(sum(b - a_ - 1 for a_, b in sauts)),
                         "sequence": "".join(UN[r["nom"]] for r in residus)})
    info["ecartes"] = ecartes
    return {"chaines": resultat, "info": info}, refus


def ecrire_pdb(res, chemin):
    lettres = "ABCDEFGH"
    L, n = [], 0
    for k, ch in enumerate(res["chaines"]):
        for o, r in enumerate(ch["residus"], start=1):
            for (nm, el, x, y, z, oc, b) in r["atomes"]:
                n += 1
                an = nm if len(nm) == 4 else " " + nm
                L.append(f"ATOM  {n % 100000:5d} {an:<4s} {r['nom']:>3s} {lettres[k]}{o:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{oc:6.2f}{b:6.2f}          {el:>2s}")
        L.append("TER")
    L.append("END")
    Path(chemin).write_text("\n".join(L) + "\n")


def table_ordinal(res):
    return [{"chaine_sortie": "AB"[k], "auth_asym_id": ch["auth"],
             "ordinal_vers_cristal": [[o, r["auth_seq"], r["icode"], r["nom"], r["label_seq"]] for o, r in enumerate(ch["residus"], start=1)]} for k, ch in enumerate(res["chaines"])]


def contacts_residus(ra, rb, seuil=5.0):
    """nombre de couples de residus (atomes lourds < seuil A) entre deux chaines extraites ; appariement ORDINAL (indice dans la liste)."""
    from scipy.spatial import cKDTree
    xa = [(i, a[2:5]) for i, r in enumerate(ra["residus"]) for a in r["atomes"]]
    xb = [(j, a[2:5]) for j, r in enumerate(rb["residus"]) for a in r["atomes"]]
    if not xa or not xb:
        return 0
    ta, tb = cKDTree(np.array([p for _, p in xa])), cKDTree(np.array([p for _, p in xb]))
    s = set()
    for i, l in enumerate(ta.query_ball_tree(tb, seuil)):
        for j in l:
            s.add((xa[i][0], xb[j][0]))
    return len(s)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(__doc__); sys.exit(2)
    cif, chs, out_pdb, out_json = sys.argv[1], sys.argv[2].split(","), sys.argv[3], sys.argv[4]
    res, refus = extraire(cif, chs)
    rap = {"cif": str(cif), "sha256_cif": sha256(cif), "chaines_auth_demandees": chs, "info": res["info"], "refus": refus,
           "chaines": [{"auth": c["auth"], "n_residus": c["n_residus"], "sequence": c["sequence"], "premier_label_seq": c["premier_label_seq"], "dernier_label_seq": c["dernier_label_seq"]} for c in res["chaines"]]}
    if refus:
        rap["etat"] = "REFUSE"
        Path(out_json).write_text(json.dumps(rap, indent=1, ensure_ascii=False) + "\n")
        print("REFUSE", [r["classe"] for r in refus]); sys.exit(3)
    ecrire_pdb(res, out_pdb)
    rap["etat"] = "ECRIT"; rap["sha256_pdb"] = sha256(out_pdb); rap["table_ordinal"] = table_ordinal(res)
    Path(out_json).write_text(json.dumps(rap, indent=1, ensure_ascii=False) + "\n")
    print("ECRIT", out_pdb, [(c["auth"], c["n_residus"]) for c in res["chaines"]])
