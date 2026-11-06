#!/usr/bin/env python3
"""auditer_28.py — l'etape de curation, EXERCEE exhaustivement sur les 28 systemes actifs, sur CPU.

Perimetre : manifeste `corpus/manifeste_corpus.json` (lecture seule), `statut == "actif_v2"`, moins
1ZM4 (DDE) et 1FQ1 (TPO), hors perimetre. Ces deux-la passent QUAND MEME dans le banc, a part, comme
TEMOINS NEGATIFS : la porte doit reproduire sur eux les messages litteraux des jobs GPU 11027, 11028.
Entrees par systeme : `corpus/systemes/<PDB>/rep{1,2,3}/structure.pdb` (image 0 de chaque replique) et,
pour le socle, `<PDB>_input.pdb`. Chaque entree : curation, porte CPU sur le x0 cure, et porte CPU sur la
structure BRUTE (temoin negatif). Un refus est une ENTREE du resultat, jamais un saut.
Ne charge aucun modele, ne demande aucun GPU.
"""
import collections
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
os.environ["CUDA_VISIBLE_DEVICES"] = ""
RACINE = Path("/gpfs/phherman/bioemu2_bench")
sys.path.insert(0, str(RACINE / "prep_fix"))
import curer_entree as ce

SORTIE = RACINE / "prep_fix" / "curation"
X0 = SORTIE / "x0"
DIAG = SORTIE / "diag"
HORS_PERIMETRE = ("1ZM4", "1FQ1")


def main():
    L = [json.loads(l) for l in open(RACINE / "corpus" / "manifeste_corpus.json") if l.strip()]
    actifs = [x for x in L if x.get("statut") == "actif_v2"]
    systemes = [x for x in actifs if x["pdb"] not in HORS_PERIMETRE]
    temoins = [x for x in actifs if x["pdb"] in HORS_PERIMETRE]
    print(f"manifeste : {len(L)} lignes ; actif_v2 = {len(actifs)} ; perimetre = {len(systemes)} ; temoins negatifs hors perimetre = {[t['pdb'] for t in temoins]}")
    contrat = ce.contrat_standard()
    print("contrat (BioEmu2) : %d residus standard, atomes lourds par residu : %s" % (len(contrat), {k: len(v) for k, v in sorted(contrat.items())}))
    resultats = []
    for lot, liste in (("perimetre", systemes), ("temoin_negatif", temoins)):
        for x in liste:
            attendu = {c["id"]: {"n_residus": c["n_residus"], "resseq_premier": c["resseq_premier"], "resseq_dernier": c["resseq_dernier"]} for c in x["chaines"]}
            # x0 « image 0 de la replique r » = coordonnees de l'image 0 du xtc, noms du structure.pdb (pilote 6CWG, preparer.py:34) ;
            # structure.pdb seul = coordonnees du PDB (NON egales a l'image 0 pour 36 repliques sur 84, manifeste) ;
            # input = `<PDB>_input.pdb` du socle (le « cristal » du pilote).
            entrees = []
            for r in x["repliques"]:
                d_r = Path(x["dossier"]) / f"rep{r}"
                entrees.append((f"{x['pdb']}__img0_rep{r}", d_r / "structure.pdb", r, d_r / "trajectory.xtc"))
                entrees.append((f"{x['pdb']}__structure_rep{r}", d_r / "structure.pdb", r, None))
            inp = (x.get("fichiers_staging_existants") or {}).get("input_pdb")
            if inp:
                entrees.append((f"{x['pdb']}__input", Path(inp), None, None))
            for etiquette, chemin, rep, xtc in entrees:
                if not chemin.is_file():
                    resultats.append({"pdb": x["pdb"], "etiquette": etiquette, "lot": lot, "kind": etiquette.split("__")[1].split("_rep")[0], "replique": rep, "etat": "REFUSE",
                                      "constats": [{"classe": "FICHIER_ABSENT", "gravite": "REFUS", "n": 1, "detail": str(chemin), "message": None, "exemples": None}],
                                      "classes_refus": ["FICHIER_ABSENT"]})
                    continue
                R = ce.curer(chemin, sortie_x0=X0 / f"{etiquette}.pdb", attendu_corpus=attendu,
                             ecrire_diag=DIAG / f"{etiquette}.pdb", contrat=contrat, xtc=xtc, image=0)
                R.update(pdb=x["pdb"], etiquette=etiquette, lot=lot, replique=rep, kind=etiquette.split("__")[1].split("_rep")[0])
                try:
                    R["porte_sur_brut"] = ce.porte_sur_brut(chemin)
                except Exception as e:
                    R["porte_sur_brut"] = {"from_description": f"EXCEPTION {type(e).__name__}: {e}"[:300], "map_cif_atom_order": None}
                if rep is not None:
                    R["ecart_structure_vs_image0_xtc_nm"] = x["structure_pdb_vs_image0_du_xtc"].get(f"rep{rep}", {}).get("ecart_brut_max_nm")
                if R.get("etat") == "PRET":
                    (X0 / f"{etiquette}.json").write_text(json.dumps({k: R[k] for k in ("entree", "coordonnees", "sha256_entree", "sha256_x0", "table_chaines", "longueurs_chaines", "n_atomes_x0", "map_residus", "histidines_nom_orig", "renommages_residus_appliques", "renommages_atomes_appliques", "porte")},
                                                                    indent=1, ensure_ascii=False) + "\n")
                resultats.append(R)
    (SORTIE / "audit_entrees.json").write_text(json.dumps(resultats, indent=1, ensure_ascii=False, default=str) + "\n")

    # ───────── rapport
    par_sys = collections.OrderedDict()
    for R in resultats:
        par_sys.setdefault((R["lot"], R["pdb"]), []).append(R)
    print("\n" + "=" * 110)
    print("A. VERDICT PAR SYSTEME ET PAR ENTREE (PRET = la porte CPU passe sur le x0 cure ; REFUSE = au moins une classe de refus)")
    for (lot, pdb), rs in par_sys.items():
        etat_sys = "PRET" if all(r["etat"] == "PRET" for r in rs) else ("PARTIEL" if any(r["etat"] == "PRET" for r in rs) else "REFUSE")
        print(f"{lot:<14} {pdb}  -> {etat_sys}  " + " ; ".join(f"{r['etiquette'].split('__')[1]}:{r['etat'][0]}" + (f"[{','.join(r['classes_refus'])}]" if r['classes_refus'] else "") for r in rs))
    n_sys_pret = sum(1 for (lot, pdb), rs in par_sys.items() if lot == "perimetre" and all(r["etat"] == "PRET" for r in rs))
    n_sys_md0 = sum(1 for (lot, pdb), rs in par_sys.items() if lot == "perimetre" and all(r["etat"] == "PRET" for r in rs if r["kind"] == "img0"))
    print(f"\nsystemes du perimetre : {len(systemes)} ; tous leurs x0 PRETS : {n_sys_pret} ; tous leurs x0 image 0 de replique PRETS : {n_sys_md0}")
    ent = [r for r in resultats if r["lot"] == "perimetre"]
    print(f"entrees du perimetre : {len(ent)} ; PRET {sum(1 for r in ent if r['etat'] == 'PRET')} ; REFUSE {sum(1 for r in ent if r['etat'] == 'REFUSE')}")

    print("\n" + "=" * 110)
    print("B. COMPTE DE REJETS PAR CLASSE (perimetre ; entrees puis systemes ; temoins negatifs a part)")
    for titre, sel in (("PERIMETRE", [r for r in resultats if r["lot"] == "perimetre"]), ("TEMOINS NEGATIFS (1ZM4, 1FQ1)", [r for r in resultats if r["lot"] != "perimetre"])):
        cl_e = collections.Counter(c for r in sel for c in r["classes_refus"])
        cl_s = collections.defaultdict(set)
        for r in sel:
            for c in r["classes_refus"]:
                cl_s[c].add(r["pdb"])
        print(f"  {titre} : {len(sel)} entrees, {sum(1 for r in sel if r['etat'] == 'REFUSE')} refusees")
        for c, n in sorted(cl_e.items(), key=lambda kv: -kv[1]):
            print(f"    {c:<32} entrees {n:>3} · systemes {len(cl_s[c]):>2} {sorted(cl_s[c])}")
        if not cl_e:
            print("    (aucun refus)")
        lit = collections.Counter()
        for r in sel:
            for c in r["constats"]:
                if c["gravite"] == "REFUS" and c["message"]:
                    lit[(c["classe"], c["message"][:230])] += 1
        for (cl, m), n in lit.items():
            print(f"    message litteral x{n} [{cl}] : {m}")
    print("\n  Autres constats (gravite != REFUS), par classe ; fichiers PDB uniques du perimetre (94 : 84 structure.pdb + 10 input) :")
    gr = collections.defaultdict(lambda: [0, 0, set()])
    for r in [r for r in resultats if r.get("kind") != "img0" and r["lot"] == "perimetre"]:      # fichiers PDB uniques du perimetre
        for c in r["constats"]:
            if c["gravite"] != "REFUS":
                g = gr[(c["classe"], c["gravite"])]
                g[0] += 1
                g[1] += c["n"] or 0
                g[2].add(r["pdb"])
    for (cl, gv), (ne, nn, ss) in sorted(gr.items()):
        print(f"    {cl:<32} {gv:<11} entrees {ne:>3} · occurrences {nn:>7} · systemes {len(ss)}")

    print("\n" + "=" * 110)
    print("C. LA PORTE CPU SUR LA STRUCTURE BRUTE (temoin negatif : sequences de `topology.to_fasta()`, comme le pilote)")
    cb = collections.Counter()
    ex = {}
    for r in resultats:
        if r.get("kind") == "img0":
            continue                      # meme fichier PDB que `structure` : ne pas compter deux fois la structure brute
        p = r.get("porte_sur_brut") or {}
        if p.get("map_cif_atom_order") == "OK":
            cle = "PASSE"
        elif p.get("from_description") not in (None, "OK"):
            cle = "from_description: " + (p["from_description"] or "")[:90]
        else:
            m = (p.get("map_cif_atom_order") or "")
            cle = m[:75] if "ValueError" not in m else m.split(": ", 1)[1][:75] if ": " in m else m[:75]
            import re
            cle = re.sub(r"\d+", "N", cle)
        cb[(r["lot"], cle)] += 1
        ex.setdefault((r["lot"], cle), []).append(r["etiquette"])
    for (lot, cle), n in sorted(cb.items(), key=lambda kv: (kv[0][0], -kv[1])):
        print(f"  {lot:<14} x{n:>3}  {cle}   ex: {ex[(lot, cle)][:3]}")
    print("\n  Messages litteraux complets de la porte BRUTE sur les temoins et sur les 3 cas connus :")
    for r in resultats:
        if r["pdb"] in ("1ZM4", "1FQ1", "2AJF", "5JMO") and "porte_sur_brut" in r:
            p = r["porte_sur_brut"]
            print(f"    {r['etiquette']:<22} from_description={str(p.get('from_description'))[:140]} | map_cif_atom_order={str(p.get('map_cif_atom_order'))[:260]}")

    print("\n" + "=" * 110)
    print("D. TABLE EXPLICITE DES NOMS D'ATOMES : AVANT -> APRES l'aller-retour mdtraj (toutes entrees, atomes LOURDS)")
    uniques = [r for r in resultats if r.get("kind") != "img0"]      # chaque fichier PDB compte UNE fois (img0 = meme fichier que structure)
    up = [r for r in uniques if r["lot"] == "perimetre"]
    print(f"  base : {len(up)} fichiers PDB du perimetre ({sum(1 for r in up if r['kind'] == 'structure')} structure.pdb + {sum(1 for r in up if r['kind'] == 'input')} _input.pdb) ; temoins negatifs exclus de cette table")
    chg, inch = collections.Counter(), collections.Counter()
    for r in up:
        for a, b, c, d, n in r.get("correspondance_mdtraj_lourds", []):
            (inch if (a == b and c == d) else chg)[(a, b, c, d)] += n
    print(f"  atomes lourds qui NE CHANGENT PAS : {sum(inch.values())} atomes, {len(inch)} couples (residu, nom) distincts")
    print(f"  atomes lourds qui CHANGENT        : {sum(chg.values())} atomes, {len(chg)} correspondances distinctes :")
    for (a, b, c, d), n in sorted(chg.items()):
        print(f"    residu {a}->{b} ; atome {c}->{d} : {n}")
    chgH = collections.Counter()
    for r in up:
        for a, b, c, d, n in r.get("correspondance_mdtraj_H", []):
            if (a, c) != (b, d):
                chgH[(a, b, c, d)] += n
    print(f"  hydrogenes qui changent de nom (ecartes du x0, informatif) : {sum(chgH.values())} atomes, {len(chgH)} correspondances distinctes ; 12 premieres :")
    for (a, b, c, d), n in sorted(chgH.items(), key=lambda kv: -kv[1])[:12]:
        print(f"    residu {a}->{b} ; atome {c}->{d} : {n}")
    ren_res = collections.Counter()
    for r in up:
        for a, b, n in r.get("renommages_residus_appliques", []):
            ren_res[(a, b)] += n
    print("  renommages de RESIDUS appliques (table explicite) :", dict(((f"{a}->{b}"), n) for (a, b), n in ren_res.items()))
    ren_at = collections.Counter()
    for r in up:
        for a, b, c, n in r.get("renommages_atomes_appliques", []):
            ren_at[(a, b, c)] += n
    print("  renommages d'ATOMES appliques (table explicite) :", dict(((f"{a}:{b}->{c}"), n) for (a, b, c), n in ren_at.items()))
    print("  noms restant NON CONFORMES apres passage (ATOMES_NON_CONFORMES / ATOMES_MANQUANTS / ATOMES_DUPLIQUES / ECART_MDTRAJ / RESIDU_NON_STANDARD) :")
    nc = collections.Counter()
    for r in resultats:
        for c in r["constats"]:
            if c["classe"] in ("ATOMES_NON_CONFORMES", "ATOMES_MANQUANTS", "ATOMES_DUPLIQUES", "ECART_MDTRAJ", "RESIDU_NON_STANDARD") and c["exemples"]:
                for e in c["exemples"]:
                    nc[(c["classe"], e, r["pdb"])] += 1
    if not nc:
        print("    AUCUN sur les entrees du perimetre et des temoins" )
    for (cl, e, pdb), n in sorted(nc.items()):
        print(f"    {cl:<24} {pdb}  {e}")

    print("\n" + "=" * 110)
    print("E. TABLE COMPLETE DE mdtraj (pdbNames.xml, protein) : ce que mdtraj convertit, lu dans le code installe")
    from mdtraj.formats.pdb.pdbfile import PDBTrajectoryFile
    PDBTrajectoryFile._loadNameReplacementTables()
    rr = PDBTrajectoryFile._residueNameReplacements
    ra = PDBTrajectoryFile._atomNameReplacements
    std20 = sorted(contrat)
    print("  residus : alias -> canonique, pour les canoniques parmi les 20 :", {k: v for k, v in sorted(rr.items()) if v in std20 and k != v})
    print("  atomes LOURDS convertis (alias -> canonique) par residu standard (hydrogenes exclus) :")
    for res in std20:
        d = {k: v for k, v in sorted(ra.get(res, {}).items()) if k != v and not v.startswith("H") and not k.startswith("H")}
        if d:
            print(f"    {res} : {d}")
    print("  (atomes lourds convertis communs a toutes les proteines : O, OXT)")

    print("\n" + "=" * 110)
    print("F. PROPRIETES BRUTES CHERCHEES (altLoc, codes d'insertion, occupations, TER, numerotation, chaines, elements) — total sur les entrees")
    tot = collections.Counter()
    print(f"  (fichiers PDB uniques : {len([r for r in resultats if r.get('kind') != 'img0'])}, temoins inclus)")
    for r in [r for r in resultats if r.get("kind") != "img0"]:
        s = r.get("scan_brut")
        if not s:
            continue
        for k in ("n_altloc", "n_icode", "n_occ_lt1", "n_occ_0", "n_element_vide", "n_nonfini", "n_hetatm", "n_model"):
            tot[k] += s[k]
        tot["entrees_sans_TER"] += (s["n_ter"] == 0)
        tot["chaines_discontigues"] += len(s["chaines_discontigues"])
        tot["decroissances"] += sum(s["decroissances"].values())
        tot["resseq_doubles_consecutifs"] += sum(s["resseq_doubles_consecutifs"].values())
        tot["resseq_repetes"] += sum(s["resseq_doubles_non_consecutifs"].values())
        tot["entrees_chaines_non_sequentielles"] += (list(r.get("table_chaines", {}).items()) != [(k, k) for k in r.get("table_chaines", {})] and any(k != v for k, v in r.get("table_chaines", {}).items()))
    for k, v in sorted(tot.items()):
        print(f"    {k:<38} {v}")
    ids = collections.Counter(tuple(r.get("table_chaines", {}).items()) for r in resultats if r["lot"] == "perimetre" and r.get("kind") != "img0")
    print("  tables de chaines (ancien -> nouveau) rencontrees :")
    for k, n in sorted(ids.items(), key=lambda kv: -kv[1]):
        print(f"    {dict(k)} : {n} entrees")
    ec = [r.get("ecart_structure_vs_image0_xtc_nm") for r in resultats if r.get("ecart_structure_vs_image0_xtc_nm") is not None]
    print(f"  ecart structure.pdb vs image 0 du xtc (manifeste, nm) : max {max(ec) if ec else None} sur {len(ec)} repliques")


    suite(resultats, systemes, contrat)


def suite(resultats, systemes, contrat):
    import re
    import urllib.request
    perim = [r for r in resultats if r["lot"] == "perimetre"]
    print("\n" + "=" * 110)
    print("G. MATRICE SYSTEME x CLASSE DE BLOCAGE (cellule = entrees portant la classe / entrees du systeme ; '.' = jamais)")
    REFUS = ["FICHIER_ABSENT", "FORMAT_MODELES", "FORMAT_ILLISIBLE", "CHAINES_DISCONTIGUES", "ALTLOC", "OCCUPATION_NULLE", "COORDONNEES", "ELEMENT_INCONNU",
             "ECART_LECTURE", "ECART_MDTRAJ", "HORS_PROTEINE", "RESIDU_ISOLE_OU_SANS_DORSALE", "ACIDE_AMINE_LIBRE", "RUPTURE_CHAINE", "COHERENCE_CORPUS",
             "RESIDU_NON_STANDARD", "ATOMES_NON_CONFORMES", "ATOMES_MANQUANTS", "ATOMES_DUPLIQUES", "ATOMES_SUPERPOSES", "COMPLEXE_DISLOQUE",
             "PORTE_FROM_DESCRIPTION", "PORTE_MAP_CIF_ATOM_ORDER", "PORTE_ETENDUE", "PORTE_COORDONNEES"]
    AUTRES = ["RENOMMAGE_RESIDU", "RENOMMAGE_ATOME", "IDENTIFIANTS_CHAINE", "CODES_INSERTION", "NUMEROTATION_DECROISSANTE", "RESSEQ_DOUBLES", "RESSEQ_REPETES",
              "SAUTS_RESSEQ_LIAISON_INTACTE", "HYDROGENES", "OXT", "ELEMENT_VIDE", "OCCUPATION_PARTIELLE", "TER_ABSENT", "UNE_SEULE_CHAINE"]
    cols = REFUS + AUTRES
    lignes_tsv = ["pdb\tn_entrees\tporte_brute_passe\t" + "\t".join(cols) + "\tverdict"]
    par = collections.OrderedDict()
    for r in resultats:
        par.setdefault((r["lot"], r["pdb"]), []).append(r)
    zero_refus = {c: 0 for c in REFUS}
    for (lot, pdb), rs in par.items():
        n = len(rs)
        cell = []
        for c in cols:
            k = sum(1 for r in rs if any(x["classe"] == c for x in r["constats"]))
            cell.append(f"{k}/{n}" if k else ".")
            if c in zero_refus and k:
                zero_refus[c] += 1
        passe = sum(1 for r in rs if (r.get("porte_sur_brut") or {}).get("map_cif_atom_order") == "OK")
        verdict = "PRET" if all(r["etat"] == "PRET" for r in rs) else "REFUSE/PARTIEL"
        lignes_tsv.append(f"{pdb}\t{n}\t{passe}/{n}\t" + "\t".join(cell) + f"\t{verdict}")
    (SORTIE / "matrice_systemes_classes.tsv").write_text("\n".join(lignes_tsv) + "\n")
    print("  (matrice complete : curation/matrice_systemes_classes.tsv ; colonnes REFUS d'abord, puis renommages/renumerotations/ecartes)")
    # version texte compacte : uniquement les colonnes non vides
    utiles = [c for i, c in enumerate(cols) if any(l.split("\t")[3 + i] != "." for l in lignes_tsv[1:])]
    print("  colonnes jamais declenchees (0 systeme) :", [c for c in cols if c not in utiles])
    print("  " + f"{'systeme':<8}{'brute':>7} " + " ".join(f"{c[:14]:>14}" for c in utiles))
    for l in lignes_tsv[1:]:
        v = l.split("\t")
        print("  " + f"{v[0]:<8}{v[2]:>7} " + " ".join(f"{v[3 + cols.index(c)]:>14}" for c in utiles) + f"  {v[-1]}")
    print("  classes de REFUS declenchees par systeme (hors temoins) :", {c: n for c, n in zero_refus.items() if n})

    print("\n" + "=" * 110)
    print("H. COHERENCE DES SEQUENCES : les x0 d'un meme systeme portent-ils la meme sequence (images 0 des 3 repliques, structures, input) ?")
    for (lot, pdb), rs in par.items():
        seqs = collections.defaultdict(list)
        for r in rs:
            if r.get("map_residus"):
                seq = tuple((m["chaine"], m["nom"]) for m in r["map_residus"])
                seqs[seq].append(r["etiquette"].split("__")[1])
        if len(seqs) != 1:
            print(f"  {pdb} : {len(seqs)} sequences distinctes", {i: v for i, v in enumerate(seqs.values())})
    print("  (un systeme absent de cette liste porte UNE seule sequence sur toutes ses entrees)")

    print("\n" + "=" * 110)
    print("I. HYDROGENES : noms APRES mdtraj contre le dictionnaire des composants PDB (CCD, files.rcsb.org/ligands) — informatif, ils sont ecartes du x0")
    ccd = {}
    for res in sorted(contrat):
        try:
            txt = urllib.request.urlopen(f"https://files.rcsb.org/ligands/download/{res}.cif", timeout=30).read().decode()
        except Exception as e:
            print(f"  CCD {res} indisponible : {e}")
            continue
        noms, entetes, dans = [], [], False
        for l in txt.splitlines():
            if l.startswith("_chem_comp_atom."):
                entetes.append(l.strip()); dans = True; continue
            if dans:
                if l.startswith(("#", "loop_", "_")):
                    break
                import shlex
                v = shlex.split(l)
                noms.append(v[entetes.index("_chem_comp_atom.atom_id")])
        ccd[res] = set(noms)
    nc = collections.Counter()
    ok = 0
    for r in [r for r in resultats if r.get("kind") != "img0" and r["lot"] == "perimetre"]:
        for a, b, c, d, n in r.get("correspondance_mdtraj_H", []):
            if b in ccd and d not in ccd[b]:
                nc[(b, d)] += n
            else:
                ok += n
    print(f"  hydrogenes conformes au CCD apres mdtraj : {ok} ; NON conformes : {sum(nc.values())} ({len(nc)} couples (residu, nom)) :", dict(sorted(nc.items(), key=lambda kv: -kv[1])[:20]))

    print("\n" + "=" * 110)
    print("J. TABLE EXPLICITE DES NOMS D'ATOMES -> curation/table_noms_atomes.tsv (toutes les correspondances observees, atomes lourds)")
    agr = collections.defaultdict(lambda: [0, set(), set()])
    for r in [r for r in resultats if r.get("kind") != "img0"]:
        for a, b, c, d, n in r.get("correspondance_mdtraj_lourds", []):
            agr[(a, b, c, d)][0] += n
            agr[(a, b, c, d)][1].add(r["etiquette"])
            agr[(a, b, c, d)][2].add(r["lot"])
    lg = ["residu_avant\tresidu_apres\tatome_avant\tatome_apres\tstatut\tn_atomes\tn_fichiers\tdans_le_contrat_BioEmu2_apres\tdevenir_dans_le_x0\tlot"]
    for (a, b, c, d), (n, es, lots) in sorted(agr.items()):
        statut = "inchange" if (a == b and c == d) else ("residu_renomme" if (c == d) else "atome_renomme")
        dans = (d in contrat.get(b, [])) or d in ce.ATOMES_ECARTES
        devenir = "ecarte (OXT)" if d in ce.ATOMES_ECARTES else ("ecrit" if dans else "NON CONFORME")
        lg.append(f"{a}\t{b}\t{c}\t{d}\t{statut}\t{n}\t{len(es)}\t{'oui' if dans else 'NON'}\t{devenir}\t{'+'.join(sorted(lots))}")
    (SORTIE / "table_noms_atomes.tsv").write_text("\n".join(lg) + "\n")
    nonp = [l for l in lg[1:] if l.split("\t")[7] == "NON" and "perimetre" in l.split("\t")[9]]
    nonc = [l for l in lg[1:] if l.split("\t")[7] == "NON" and "perimetre" not in l.split("\t")[9]]
    print(f"  {len(lg) - 1} correspondances ; ne tombant PAS dans le contrat apres passage : perimetre {len(nonp)} {nonp[:5]} ; temoins negatifs seuls {len(nonc)} (atomes de DDE / TPO)")
    print("\n" + "=" * 110)
    print("K. PORTE ETENDUE (CPU) sur les x0 PRETS : etapes passees")
    pr = [r for r in perim if r.get("porte")]
    for cle in ("from_description", "map_cif_atom_order", "prepare_sampling_inputs"):
        print(f"  {cle:<26} OK sur {sum(1 for r in pr if r['porte'].get(cle) == 'OK')}/{len(pr)} x0 ecrits (perimetre)")
    print(f"  x0 topologie == systeme       {sum(1 for r in pr if r['porte'].get('x0_topologie_egale_au_systeme'))}/{len(pr)} ; coordonnees finies {sum(1 for r in pr if r['porte'].get('x0_coordonnees_finies'))}/{len(pr)} ; coordonnees mappees == x0 {sum(1 for r in pr if r.get('coordonnees_mappees_egales_au_x0'))}/{len(pr)}")
    cc = [r.get("contacts_atomes_inter_chaines_5A") for r in perim if r.get("contacts_atomes_inter_chaines_5A") is not None]
    print(f"  contacts atome-atome inter-chaines < 5 A : min {min(cc)} · mediane {sorted(cc)[len(cc) // 2]} · max {max(cc)} sur {len(cc)} x0")


if __name__ == "__main__":
    main()
