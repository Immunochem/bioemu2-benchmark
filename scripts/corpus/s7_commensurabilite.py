"""J5 — corpus v2 : le critere d'inclusion est la COMMENSURABILITE, plus la completude du cristal.

Instruction pherman : « il faut donner les memes objets aux deux (la MD et BioEmu2) ».

Critere C4' (remplace `gap_interne_rcsb`) evalue sur `structure.pdb` de la replique 1 de chaque projet :
  (a) aucune RUPTURE de chaine dans la structure MD, au sens du preparateur corrige (liaison peptidique C(i)-N(i+1)) ;
  (b) accord residu par residu entre la structure MD et ce que le preparateur donnera a BioEmu2 :
      meme nombre de residus, meme sequence que `topology.to_fasta()` de la MEME structure.
Les FONCTIONS du preparateur sont importees (copie figee dans corpus/vendor/, sha256 dans SHA_PREPARATEUR.txt), jamais recopiees.
Les residus non observes d'apres RCSB deviennent une PROPRIETE declaree (internes reconstruits / absents, terminaux), plus un rejet.

Criteres inchanges : exactement 2 chaines proteiques ; <= 1030 residus ; >= 25 paires a 5,0 A ; resolution declaree.
Un rejet est une ENTREE.
"""
import collections
import io
import json
import os
import statistics
import sys

import mdtraj as md

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE) + "/vendor")
import corpus_lib as cl  # noqa: E402
import preparer_cristal as pc  # noqa: E402  (prep_fix, version lue sur le disque)
import preparer_orig as po  # noqa: E402  (production, sha f50295da...)
import s2_analyse as s2  # noqa: E402

ROOT = cl.ROOT
CACHE = cl.CACHE
FENETRES = (2.9, 1.9)   # A : defaut de la version lue ; fenetre de la version annoncee (sensibilite)


def isoler(fn, *a, **k):
    """Appelle une fonction du preparateur qui peut `sys.exit` : le refus est une DONNEE, pas un arret."""
    try:
        return fn(*a, **k), None
    except SystemExit as e:
        return None, str(e)[:300]


def evaluer_preparateur(traj, chemin, chaines_md, fasta):
    """Pour chaque chaine proteique : residus, gaps, sequence, comparaison a to_fasta. Fenetre 2,9 puis 1,9."""
    ins = pc.codes_insertion(chemin)
    out = {}
    lettres = collections.Counter(pc.chaine_id(traj.topology, c) for c in traj.topology.chains)
    for c in chaines_md:
        lettre = c["id"]
        rec = {"lettre": lettre, "n_res_mdtraj": c["n_res"], "lettre_dupliquee": lettres[lettre] > 1}
        fa = fasta.get(c["idx"])
        rec["n_fasta"] = len(fa) if fa else None
        for w in FENETRES:
            pc.LIEN_CN_MAX_NM = w / 10
            r, err = isoler(pc.residus_de, traj, lettre)
            if err:
                rec["fenetre_%.1f" % w] = {"refus": err}
                continue
            g, err = isoler(pc.garde_numerotation, r["residus"], r["d_suivant_nm"], lettre, ins.get(lettre, []))
            if err:
                rec["fenetre_%.1f" % w] = {"refus": err}
                continue
            seq = "".join(pc.lettre_de(x.name) for x in r["residus"])
            d_all = [d * 10 for d in r["d_suivant_nm"] if d is not None]
            diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(seq, fa or "")) if a != b][:5]
            rec["fenetre_%.1f" % w] = {
                "n_obs": g["n_residus_observes"], "n_gaps": g["n_gaps"],
                "gaps": [(x["apres"], x["avant"], x["d_C_N_angstrom"]) for x in g["gaps"]],
                "n_sauts_resSeq_sans_rupture": g["n_sauts_resSeq_sans_rupture"],
                "liaison_max_A": g["liaison_peptidique_mesuree_angstrom"]["max"],
                "n_jonctions_1.9_2.9": sum(1 for d in d_all if 1.9 < d <= 2.9),
                "n_exclus_libres": len(r["exclus"]), "n_modifies": len(
                    [x for x in r["residus"] if x.name not in pc.TROIS_VERS_UN]),
                "n_X": seq.count("X"), "len_seq": len(seq),
                "accord_n_res": (len(seq) == c["n_res"]), "accord_fasta_longueur": (fa is not None and len(seq) == len(fa)),
                "accord_fasta_sequence": (fa is not None and seq == fa), "diffs_vs_fasta": diffs}
        pc.LIEN_CN_MAX_NM = FENETRES[0] / 10
        # preparateur d'ORIGINE (production) : memes questions
        ro, err = isoler(po.residus_de, traj, lettre)
        if err:
            rec["origine"] = {"refus": err}
        else:
            go, err = isoler(po.garde_numerotation, ro, lettre, ins.get(lettre, []))
            seq_o = "".join(po.TROIS_VERS_UN[x.name] for x in ro)
            rec["origine"] = {"refus": err} if err else {
                "n_obs": len(ro), "n_gaps_sauts_resSeq": go["n_gaps"], "accord_n_res": len(ro) == c["n_res"],
                "accord_fasta_sequence": (fa is not None and seq_o == fa)}
        out[lettre] = rec
    return out


def reconstruction_rcsb(E, chain, rc):
    """Residus non observes RCSB, par classe, dans les TROIS numerotations (seq_id d'entite, auteur, resSeq de la MD).

    `rc` est le resultat de s2.analyser_chaine_rcsb (instance et methode d'appariement) : on ne refait pas
    l'appariement, on relit sa cible avec les memes fonctions.
    """
    pe = inst = None
    for p in E["polymer_entities"] or []:
        for i in p["polymer_entity_instances"] or []:
            if i["rcsb_id"] == rc["instance_rcsb"]:
                pe, inst = p, i
    seq_ent = pe["entity_poly"]["pdbx_seq_one_letter_code_can"]
    mp = inst["rcsb_polymer_entity_instance_container_identifiers"]["auth_to_entity_poly_seq_mapping"] or []
    sq = s2.seq_md(chain)
    if rc["methode_appariement"] == "numerotation":
        pos = s2.positions_par_numerotation(chain, sq, seq_ent, mp)[0]
    else:
        pos = s2.positions_par_alignement(sq, seq_ent)[0]
    U = s2.non_observes(inst)
    L = len(seq_ent)
    obs = [p for p in range(1, L + 1) if p not in U]
    out = {"instance_rcsb": rc["instance_rcsb"], "methode_appariement": rc["methode_appariement"], "longueur_entite": L,
           "n_non_observes_total": len(U), "evaluable": bool(obs)}
    if not obs:
        return out
    f, l = min(obs), max(obs)
    mdres = {}
    for i, p in enumerate(pos):
        if p:
            mdres[p] = (chain["resseq"][i], chain["resnames"][i], i + 1)
    inv = {i + 1: v for i, v in enumerate(mp)}
    mdpos = [p for p in pos if p]
    lo, hi = min(mdpos), max(mdpos)
    classes = {"internes": sorted(p for p in U if f < p < l), "n_terminaux": sorted(p for p in U if p < f),
               "c_terminaux": sorted(p for p in U if p > l)}
    out["premier_observe_seq_id"], out["dernier_observe_seq_id"] = f, l
    out["definition_v1_internes_dans_etendue_MD"] = len([p for p in U if lo < p < hi])
    for nom, lst in classes.items():
        runs = []
        for p in lst:
            if runs and p == runs[-1][1] + 1:
                runs[-1][1] = p
            else:
                runs.append([p, p])
        out[nom] = {
            "n": len(lst), "presents_dans_MD": sum(1 for p in lst if p in mdres), "absents_de_MD": sum(1 for p in lst if p not in mdres),
            "runs": [{"seq_id": [a, b], "auth_seq_id": [inv.get(a), inv.get(b)],
                      "resSeq_MD": [mdres[a][0] if a in mdres else None, mdres[b][0] if b in mdres else None],
                      "ordinal_MD": [mdres[a][2] if a in mdres else None, mdres[b][2] if b in mdres else None],
                      "residus_MD": "".join(s2.AA3.get(mdres[p][1], "X") if p in mdres else "-" for p in range(a, b + 1)),
                      "residus_entite": seq_ent[a - 1:b]} for a, b in runs]}
    return out


def main():
    ld = cl.charger_projets_dynabench()
    v1 = {P["accession"]: P for P in cl.jload(ROOT + "/inventaire_dynabench.json")["projets"]}
    rc = cl.jload(CACHE + "/rcsb_graphql.json")
    rcent, _ko = s2.index_rcsb(rc)
    sel1 = cl.jload(ROOT + "/selection.json")["selection"]
    sel1_acc = {s["accession"] for s in sel1}
    print("ENTREES projets=%d  v1_selection=%d  preparateur=%s  origine=%s" % (
        len(v1), len(sel1), cl.sha256_file(HERE + "/../vendor/preparer_cristal.py")[:16], cl.sha256_file(HERE + "/../vendor/preparer_orig.py")[:16]))
    res = []
    for q in ld["projets"]:
        acc = q["accession"]
        P = json.loads(json.dumps(v1[acc]))
        pdb = P["pdb"]
        chemin = "%s/structure_rep1/%s_%s.pdb" % (CACHE, pdb, acc)
        P["motifs_rejet_v1"] = list(P["motifs_rejet"])
        P["motifs_rejet"] = []
        P["notes_v2"] = []
        S = s2.analyser_structure(chemin)      # meme fonction qu'en v1 : chaines, resSeq, noms
        traj = md.load(chemin)
        fasta_list = traj.topology.to_fasta()
        fasta = {ch.index: fasta_list[k] for k, ch in enumerate(traj.topology.chains)} if len(fasta_list) == traj.topology.n_chains else {}
        prot = [c for c in S["chaines"] if c["n_res_prot"] > 0 and c["n_res_prot"] >= 0.9 * c["n_res"]]
        P["commensurabilite"] = evaluer_preparateur(traj, chemin, prot, fasta)
        # ---- RCSB : propriete declaree
        E = rcent.get(pdb)
        P["rcsb_v2"] = {"chaines": []}
        nonapp = False
        if E is None:
            P["motifs_rejet"].append("rcsb_indisponible")
        else:
            for c in prot:
                r1 = s2.analyser_chaine_rcsb(E, c)
                if not r1["evaluable"]:
                    nonapp = True
                    P["rcsb_v2"]["chaines"].append({"chaine": c["id"], "evaluable": False, "raison": r1.get("raison")})
                    continue
                d = reconstruction_rcsb(E, c, r1)
                d["chaine"] = c["id"]
                d["evaluable"] = d.get("evaluable", False)
                P["rcsb_v2"]["chaines"].append(d)
            if nonapp:
                P["motifs_rejet"].append("rcsb_non_apparie")
        # ---- criteres v2
        n_prot = len(prot)
        P["n_chaines_proteiques"] = n_prot
        if n_prot != 2:
            P["motifs_rejet"].append("chaines_proteiques_!=2")
        total = sum(c["n_res_prot"] for c in prot)
        P["n_residus_total"] = total
        if total > cl.MAX_RES:
            P["motifs_rejet"].append("taille_>1030")
        if n_prot == 2 and P["interface_paires_5A"] is not None and P["interface_paires_5A"] < cl.MIN_PAIRS:
            P["motifs_rejet"].append("interface_<25")
        if P["rcsb"] is not None and P["rcsb"]["resolution_A"] is None:
            P["motifs_rejet"].append("resolution_non_declaree")
        W = "fenetre_%.1f" % FENETRES[0]
        refus = [k for k, v in P["commensurabilite"].items() if "refus" in v.get(W, {})]
        rupt = [k for k, v in P["commensurabilite"].items() if v.get(W, {}).get("n_gaps", 0) > 0]
        incomm = []
        for k, v in P["commensurabilite"].items():
            w = v.get(W, {})
            if "refus" in w:
                continue
            if not (w.get("accord_n_res") and w.get("accord_fasta_longueur") and w.get("accord_fasta_sequence")):
                incomm.append(k)
        P["commensurabilite_resume"] = {"chaines_refusees_par_le_preparateur": refus, "chaines_a_rupture_CN": rupt,
                                        "chaines_en_desaccord_avec_to_fasta": incomm}
        if rupt:
            P["motifs_rejet"].append("rupture_chaine_MD")
        if refus or incomm:
            P["motifs_rejet"].append("non_commensurable")
        P["passe_les_criteres_v2"] = not P["motifs_rejet"]
        # evaluations alternatives (sensibilite, ne decident rien)
        W2 = "fenetre_%.1f" % FENETRES[1]

        def ok_prep(cle):
            for v in P["commensurabilite"].values():
                w = v.get(cle, {})
                if "refus" in w:
                    return False
                if w.get("n_gaps", 0) > 0 or not (w.get("accord_n_res") and w.get("accord_fasta_sequence")):
                    return False
            return bool(P["commensurabilite"])
        P["commensurable_fenetre_1.9"] = ok_prep(W2)
        o = True
        for v in P["commensurabilite"].values():
            w = v.get("origine", {})
            if "refus" in w or w.get("n_gaps_sauts_resSeq", 0) > 0 or not (w.get("accord_n_res") and w.get("accord_fasta_sequence")):
                o = False
        P["commensurable_preparateur_origine"] = o and bool(P["commensurabilite"])
        res.append(P)

    # ============ comptes ============
    n_pass = sum(1 for P in res if P["passe_les_criteres_v2"])
    ordre = ["structure_indisponible", "rcsb_indisponible", "rcsb_non_apparie", "chaines_proteiques_!=2", "taille_>1030", "interface_<25",
             "resolution_non_declaree", "rupture_chaine_MD", "non_commensurable"]
    tous = collections.Counter(m for P in res for m in P["motifs_rejet"])
    premier = collections.Counter()
    for P in res:
        for m in ordre:
            if m in P["motifs_rejet"]:
                premier[m] += 1
                break
    print("COMPTE_V2 total=%d  passent_les_5_criteres=%d  rejetes=%d  (v1: passent=%d)" % (
        len(res), n_pass, len(res) - n_pass, sum(1 for P in res if not P["motifs_rejet_v1"])))
    for m in ordre:
        print("MOTIF_V2_TOUS    %-26s %3d" % (m, tous.get(m, 0)))
    for m in ordre:
        print("MOTIF_V2_PREMIER %-26s %3d" % (m, premier.get(m, 0)))
    print("MOTIF_V2_PREMIER somme=%d (doit egaler rejetes=%d)" % (sum(premier.values()), len(res) - n_pass))
    surv = list(res)
    print("ENTONNOIR_V2 depart %d" % len(surv))
    for nom, m in (("rcsb_disponible_et_apparie", ("rcsb_indisponible", "rcsb_non_apparie")), ("2_chaines_proteiques", ("chaines_proteiques_!=2",)),
                   ("taille<=1030", ("taille_>1030",)), ("interface>=25", ("interface_<25",)), ("resolution_declaree", ("resolution_non_declaree",)),
                   ("aucune_rupture_chaine_MD", ("rupture_chaine_MD",)), ("commensurable", ("non_commensurable",))):
        surv = [P for P in surv if not any(x in P["motifs_rejet"] for x in m)]
        print("ENTONNOIR_V2 apres %-28s %3d" % (nom, len(surv)))
    # volet (a) et (b) sur les 210
    nrupt = sum(1 for P in res if P["commensurabilite_resume"]["chaines_a_rupture_CN"])
    nrefus = sum(1 for P in res if P["commensurabilite_resume"]["chaines_refusees_par_le_preparateur"])
    ndes = sum(1 for P in res if P["commensurabilite_resume"]["chaines_en_desaccord_avec_to_fasta"])
    print("VOLET_A projets_avec_rupture_CN(fenetre 2.9)=%d  chaines_refusees_par_le_preparateur(projets)=%d" % (nrupt, nrefus))
    print("VOLET_B projets_en_desaccord_avec_to_fasta=%d sur %d" % (ndes, len(res)))
    nch = nch_ok = jmax = 0
    jonc = collections.Counter()
    for P in res:
        for k, v in P["commensurabilite"].items():
            w = v.get("fenetre_2.9", {})
            if "refus" in w:
                continue
            nch += 1
            nch_ok += bool(w["accord_n_res"] and w["accord_fasta_sequence"] and w["n_gaps"] == 0)
            jonc["jonctions_1.9_2.9"] += w["n_jonctions_1.9_2.9"]
            jonc["chaines_avec_n_X"] += (w["n_X"] > 0)
            jonc["chaines_avec_modifies"] += (w["n_modifies"] > 0)
            jonc["residus_modifies_total"] += w["n_modifies"]
            jonc["residus_X_total"] += w["n_X"]
            jonc["chaines_exclus_libres"] += (w["n_exclus_libres"] > 0)
            jonc["sauts_resSeq_sans_rupture"] += w["n_sauts_resSeq_sans_rupture"]
            jmax = max(jmax, w["liaison_max_A"] or 0)
    print("CHAINES evaluees_par_le_preparateur=%d  accord_complet_et_sans_rupture=%d  liaison_C-N_max_A=%.3f  %s" % (nch, nch_ok, jmax, json.dumps(dict(jonc))))
    print("DEPENDANCE_VERSION commensurables: fenetre_2.9(version lue)=%d  fenetre_1.9(version annoncee)=%d  preparateur_d_origine(production)=%d  sur %d" % (
        sum(1 for P in res if not any(m in P["motifs_rejet"] for m in ("rupture_chaine_MD", "non_commensurable")) and P["commensurabilite"]),
        sum(1 for P in res if P["commensurable_fenetre_1.9"]), sum(1 for P in res if P["commensurable_preparateur_origine"]), len(res)))
    for P in res:
        if P["commensurable_preparateur_origine"] != (not any(m in P["motifs_rejet"] for m in ("rupture_chaine_MD", "non_commensurable")) and bool(P["commensurabilite"])):
            ch = {k: v.get("origine") for k, v in P["commensurabilite"].items()}
            print("ORIGINE_DIFFERE %s %s origine=%s version_lue=%s" % (P["pdb"], P["accession"], json.dumps(ch), json.dumps(P["commensurabilite_resume"])))
    for P in res:
        for k, v in P["commensurabilite"].items():
            w = v.get("fenetre_2.9", {})
            if "refus" in w or w.get("n_gaps") or not w.get("accord_fasta_sequence") or not w.get("accord_n_res"):
                print("DESACCORD_OU_RUPTURE %s %s chaine=%s %s" % (P["pdb"], P["accession"], k, json.dumps(w if "refus" not in w else w, ensure_ascii=False)[:400]))
    for P in res:
        for k, v in P["commensurabilite"].items():
            if v.get("fenetre_2.9", {}).get("n_jonctions_1.9_2.9"):
                print("JONCTION_ETIREE %s %s chaine=%s n=%d liaison_max=%s" % (P["pdb"], P["accession"], k, v["fenetre_2.9"]["n_jonctions_1.9_2.9"], v["fenetre_2.9"]["liaison_max_A"]))

    # ============ residus non observes RCSB : propriete du corpus ============
    def agreg(filtre, etiq):
        nc = ni = nin = npres = nabs = ncrec = nterm_n = nterm_c = nterm_n_pres = nterm_c_pres = old = nchaines_term = 0
        chaines_rec = []
        for P in res:
            if not filtre(P):
                continue
            for c in P["rcsb_v2"]["chaines"]:
                if not c.get("evaluable"):
                    continue
                nc += 1
                I, N, C_ = c["internes"], c["n_terminaux"], c["c_terminaux"]
                old += c["definition_v1_internes_dans_etendue_MD"]
                if I["n"]:
                    ni += 1
                    nin += I["n"]
                    npres += I["presents_dans_MD"]
                    nabs += I["absents_de_MD"]
                    if I["presents_dans_MD"]:
                        ncrec += 1
                        chaines_rec.append((P["pdb"], c["chaine"]))
                nterm_n += N["n"]
                nterm_c += C_["n"]
                nterm_n_pres += N["presents_dans_MD"]
                nterm_c_pres += C_["presents_dans_MD"]
                nchaines_term += bool(N["n"] or C_["n"])
        print("RECONSTRUCTION[%s] chaines_evaluables=%d  chaines_avec_residus_internes_non_observes=%d  residus_internes_non_observes=%d  "
              "dont_presents_dans_la_MD(reconstruits)=%d  dont_absents_de_la_MD=%d  chaines_avec_au_moins_un_reconstruit=%d  "
              "[ancienne_definition_etendue_MD: %d]" % (etiq, nc, ni, nin, npres, nabs, ncrec, old))
        print("TERMINAUX[%s] chaines_avec_residus_terminaux_non_observes=%d  N-terminaux=%d (presents_dans_MD=%d)  C-terminaux=%d (presents_dans_MD=%d)" % (
            etiq, nchaines_term, nterm_n, nterm_n_pres, nterm_c, nterm_c_pres))
        return chaines_rec
    agreg(lambda P: True, "210_projets")
    agreg(lambda P: P["n_chaines_proteiques"] == 2, "2_chaines")

    # ============ selection v2 : meme regle qu'en v1 ============
    CIBLE = 30
    pas = [P for P in res if P["passe_les_criteres_v2"]]
    ecartes_cat = [P for P in pas if P["categorie_code"] is None]
    elig = [P for P in pas if P["categorie_code"] is not None]
    sel, raison = [], {}

    def ajouter(P, why):
        if P["accession"] not in raison:
            sel.append(P)
            raison[P["accession"]] = why
    for P in sorted([p for p in elig if p["dans_staging12"]], key=lambda p: (p["n_residus_total"], p["pdb"])):
        ajouter(P, "deja_sur_staging")
    for P in sorted([p for p in elig if p["vhh"]], key=lambda p: (p["n_residus_total"], p["pdb"])):
        ajouter(P, "vhh")
    reste = sorted([p for p in elig if p["accession"] not in raison], key=lambda p: (p["n_residus_total"], p["pdb"], p["accession"]))
    for P in reste:
        if len(sel) >= CIBLE:
            break
        ajouter(P, "plus_petits")
    reserve = [p for p in reste if p["accession"] not in raison][:10]
    for P in res:
        P["selectionne_v2"] = P["accession"] in raison
        P["raison_selection_v2"] = raison.get(P["accession"])
        P["statut_vs_v1"] = ("reste" if (P["selectionne_v2"] and P["accession"] in sel1_acc) else "entre" if P["selectionne_v2"]
                             else "sort" if P["accession"] in sel1_acc else None)
    entrent = [P for P in res if P["statut_vs_v1"] == "entre"]
    sortent = [P for P in res if P["statut_vs_v1"] == "sort"]
    restent = [P for P in res if P["statut_vs_v1"] == "reste"]
    print("SELECTION_V2 n=%d  (deja_staging=%d, vhh=%d, plus_petits=%d)  cible=%d  passent=%d  eligibles_avec_categorie=%d  ecartes_categorie_inconnue=%d" % (
        len(sel), sum(1 for v in raison.values() if v == "deja_sur_staging"), sum(1 for v in raison.values() if v == "vhh"),
        sum(1 for v in raison.values() if v == "plus_petits"), CIBLE, len(pas), len(elig), len(ecartes_cat)))
    print("DELTA_V1_V2 entrent=%d  sortent=%d  restent=%d  (v1=%d, v2=%d)" % (len(entrent), len(sortent), len(restent), len(sel1), len(sel)))
    for P in sorted(entrent, key=lambda p: p["n_residus_total"]):
        print("ENTRE %s %s raison=%s total=%d cat=%s cls=%s vhh=%s staging12=%s" % (P["pdb"], P["accession"], P["raison_selection_v2"], P["n_residus_total"], P["categorie_code"], P["classe_difficulte"], P["vhh"], P["dans_staging12"]))
    for P in sorted(sortent, key=lambda p: p["n_residus_total"]):
        print("SORT %s %s total=%d raison_v1=%s cat=%s cls=%s (passe_v2=%s)" % (P["pdb"], P["accession"], P["n_residus_total"], P["raison_selection"], P["categorie_code"], P["classe_difficulte"], P["passe_les_criteres_v2"]))
    for P in sorted(restent, key=lambda p: p["n_residus_total"]):
        print("RESTE %s %s raison_v2=%s raison_v1=%s total=%d" % (P["pdb"], P["accession"], P["raison_selection_v2"], P["raison_selection"], P["n_residus_total"]))
    for P in reserve:
        print("RESERVE_V2 %s %s total=%d cat=%s cls=%s" % (P["pdb"], P["accession"], P["n_residus_total"], P["categorie_code"], P["classe_difficulte"]))
    for P in ecartes_cat:
        print("ECARTE_CATEGORIE_INCONNUE_V2 %s %s total=%d" % (P["pdb"], P["accession"], P["n_residus_total"]))
    for pdb in "6CWG 3K75 1ZM4 1BKD 2AJF 2C0L 2OT3".split():
        P = next(p for p in res if p["pdb"] == pdb)
        print("STAGING12_V2 %s passe_v2=%s motifs=%s total=%d" % (pdb, P["passe_les_criteres_v2"], ",".join(P["motifs_rejet"]) or "-", P["n_residus_total"]))
    print("COMPOSITION_V2 categorie=%s difficulte=%s" % (
        json.dumps(dict(sorted(collections.Counter(str(P["categorie_code"]) for P in sel).items()))),
        json.dumps(dict(sorted(collections.Counter(str(P["classe_difficulte"]) for P in sel).items())))))
    tail = [P["n_residus_total"] for P in sel]
    print("TAILLE_V2 n=%d min=%d mediane=%.1f max=%d somme=%d  VHH=%d" % (len(tail), min(tail), statistics.median(tail), max(tail), sum(tail), sum(1 for P in sel if P["vhh"])))
    # ---- 6CWG : les trois numerotations d'un meme residu
    P = next(p for p in res if p["pdb"] == "6CWG")
    for c in P["rcsb_v2"]["chaines"]:
        if c.get("evaluable"):
            print("NUMEROTATION_6CWG chaine=%s internes=%s n_terminaux=%s c_terminaux=%s" % (c["chaine"], json.dumps(c["internes"], ensure_ascii=False), json.dumps(c["n_terminaux"], ensure_ascii=False), json.dumps(c["c_terminaux"], ensure_ascii=False)))

    # ============ ecritures ============
    cl.jdump({"meta": {"date": cl.now_iso(), "criteres": "v2 : 2 chaines proteiques ; <=1030 ; >=25 paires 5,0 A ; resolution declaree ; commensurabilite (aucune rupture C-N + accord to_fasta)",
                       "preparateur_sha256": cl.sha256_file(HERE + "/../vendor/preparer_cristal.py"), "preparateur_origine_sha256": cl.sha256_file(HERE + "/../vendor/preparer_orig.py"),
                       "n_projets": len(res), "n_passent": n_pass}, "projets": res}, ROOT + "/inventaire_dynabench_v2.json")
    cl.jdump({"cible": CIBLE, "selection": [{"accession": P["accession"], "pdb": P["pdb"], "raison": raison[P["accession"]], "n_residus_total": P["n_residus_total"],
                                              "statut_vs_v1": P["statut_vs_v1"]} for P in sel],
              "sortent_de_la_selection_v1": [{"accession": P["accession"], "pdb": P["pdb"]} for P in sortent],
              "reserve": [{"accession": P["accession"], "pdb": P["pdb"], "n_residus_total": P["n_residus_total"]} for P in reserve],
              "ecartes_categorie_inconnue": [{"accession": P["accession"], "pdb": P["pdb"]} for P in ecartes_cat]}, ROOT + "/selection_v2.json")
    cl.jdump({"tous_motifs": dict(tous), "premier_motif": dict(premier), "ordre": ordre,
              "rejets": [{"accession": P["accession"], "pdb": P["pdb"], "motifs": P["motifs_rejet"]} for P in res if P["motifs_rejet"]]}, ROOT + "/rejets_par_motif_v2.json")
    cols = ["pdb", "accession", "categorie_code", "classe_difficulte", "n_chaines_proteiques", "n_residus_total", "interface_paires_5A", "passe_v2", "motifs_v2", "motifs_v1",
            "selectionne_v2", "raison_v2", "statut_vs_v1", "vhh"]
    with open(ROOT + "/inventaire_dynabench_v2.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for P in res:
            f.write("\t".join(str(x) for x in [P["pdb"], P["accession"], P["categorie_code"], P["classe_difficulte"], P["n_chaines_proteiques"], P["n_residus_total"],
                                               P["interface_paires_5A"], P["passe_les_criteres_v2"], ";".join(P["motifs_rejet"]), ";".join(P["motifs_rejet_v1"]),
                                               P["selectionne_v2"], P["raison_selection_v2"], P["statut_vs_v1"], P["vhh"]]) + "\n")
    print("FINI s7 %s" % cl.now_iso())


if __name__ == "__main__":
    main()
