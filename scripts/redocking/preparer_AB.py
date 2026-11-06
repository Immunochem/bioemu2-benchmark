#!/usr/bin/env python3
"""preparer_AB.py — CPU. Prepare les x0 CRISTALLINS du redocking BioEmu2, chacun passe par la curation (prep_fix/curer_entree.py, IMPORTEE).

TROIS STRATES (la date de DEPOT decide, pas une etiquette) :
  A_corpus  : les 28 systemes du corpus (actif_v2 hors 1ZM4, 1FQ1). Deux chaines d interface choisies parmi les copies de l unite asymetrique (sans refus d extraction d abord, puis le plus
              de contacts). Date de depot relevee a l API RCSB. Beaucoup sont refuses par la preparation (voir manifeste) : c est un RESULTAT, pas un incident.
  A_apparie : complexes RCSB deposes AVANT 2021-01-01 (marge de 9 mois sous la coupure presumee), tires par LA MEME procedure que B (memes filtres RCSB, meme extraction, meme curation,
              memes refus) : seule la date differe. C est la strate contaminee comparable.
  B         : complexes RCSB deposes APRES 2021-09-30 (strate tenue a l ecart).
PROCEDURE RCSB commune a A_apparie et B : 2 chaines proteiques / 2 entites / 2 instances deposees, X-ray <= 3,0 A, 150-835 residus deposes ; >= 30 residus par chaine, >= 95 % modelises,
sequences sur les 20 standard ; UN seul complexe par groupe de sequences a 30 % d identite (aucune famille deux fois dans une strate) ; tirage STRATIFIE par periode de depot, graine fixe ;
telechargement UN cristal a la fois (pause 1 s) ; extraction ; >= 25 couples de residus en contact (< 5 A) ; curation. Tout rejet est une ENTREE du manifeste, avec sa classe.
"""
import collections, hashlib, importlib.util, itertools, json, os, random, sys, time, urllib.request
from pathlib import Path
sys.dont_write_bytecode = True
os.environ["CUDA_VISIBLE_DEVICES"] = ""
B = Path("/gpfs/phherman/bioemu2_bench"); R = B / "redocking"
sys.path.insert(0, str(R / "scripts"))
import extraire_cristal as E

SEUIL_CONTACTS = 25
GRAINE = 20261105
HORS = {"1ZM4", "1FQ1"}
STRATES_RCSB = {
    "B": {"operateur": "greater", "date": "2021-09-30T23:59:59Z", "cible": 30, "max_traites": 160, "max_meta": None,
          "bornes": (("2022-07-01", "P1_2021-10_2022-06"), ("2023-04-01", "P2_2022-07_2023-03"), ("2024-07-01", "P3_2023-04_2024-06"), ("2025-07-01", "P4_2024-07_2025-06")), "dernier": "P5_2025-07_2026"},
    "A_apparie": {"operateur": "less_or_equal", "date": "2020-12-31T23:59:59Z", "cible": 30, "max_traites": 160, "max_meta": 1800,
                  "bornes": (("2009-01-01", "Q1_avant_2009"), ("2012-01-01", "Q2_2009_2011"), ("2015-01-01", "Q3_2012_2014"), ("2018-01-01", "Q4_2015_2017")), "dernier": "Q5_2018_2020"},
}


def charger_curer():
    sp = importlib.util.spec_from_file_location("curer_entree", B / "prep_fix/curer_entree.py")
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def post_json(url, obj):
    req = urllib.request.Request(url, data=json.dumps(obj).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read())


def dates_rcsb(ids):
    q = "query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id rcsb_accession_info{deposit_date initial_release_date} } }"
    r = post_json("https://data.rcsb.org/graphql", {"query": q, "variables": {"ids": list(ids)}})
    return {e["rcsb_id"]: e["rcsb_accession_info"] for e in (r["data"]["entries"] or []) if e}


GQL = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
 rcsb_accession_info{deposit_date initial_release_date}
 rcsb_entry_info{resolution_combined deposited_polymer_monomer_count deposited_modeled_polymer_monomer_count}
 struct{title}
 polymer_entities{ rcsb_id entity_poly{pdbx_seq_one_letter_code_can rcsb_entity_polymer_type pdbx_strand_id}
   rcsb_polymer_entity_group_membership{aggregation_method group_id similarity_cutoff}
   rcsb_entity_source_organism{ncbi_scientific_name}
   rcsb_polymer_entity{pdbx_description} } } }"""


def candidats(nom, cfg):
    """candidats RCSB de la strate (fichier cache strate_<nom>/candidats.json ; B est deja produit par choisir_B.py)."""
    d = R / f"strate_{nom}"
    d.mkdir(parents=True, exist_ok=True)
    f = d / "candidats.json"
    if f.is_file():
        return json.loads(f.read_text())
    T = lambda a, o, v: {"type": "terminal", "service": "text", "parameters": {"attribute": a, "operator": o, "value": v}}
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        T("rcsb_accession_info.deposit_date", cfg["operateur"], cfg["date"]),
        T("rcsb_entry_info.polymer_entity_count_protein", "equals", 2), T("rcsb_entry_info.polymer_entity_count", "equals", 2),
        T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 2), T("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        T("rcsb_entry_info.resolution_combined", "less_or_equal", 3.0),
        T("rcsb_entry_info.deposited_polymer_monomer_count", "range", {"from": 150, "include_lower": True, "to": 835, "include_upper": True})]},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 1000}, "sort": [{"sort_by": "rcsb_accession_info.deposit_date", "direction": "asc"}]}}
    ids, start, total = [], 0, None
    while True:
        q["request_options"]["paginate"]["start"] = start
        r = post_json("https://search.rcsb.org/rcsbsearch/v2/query", q)
        total = r["total_count"]; ids += [x["identifier"] for x in r["result_set"]]
        start += 1000
        if start >= total or start >= 10000:
            break
        time.sleep(1)
    n_ids_total = total
    if cfg["max_meta"] and len(ids) > cfg["max_meta"]:
        random.Random(GRAINE).shuffle(ids)
        ids = ids[:cfg["max_meta"]]
    cand, rej = {}, []
    for i in range(0, len(ids), 100):
        r = post_json("https://data.rcsb.org/graphql", {"query": GQL, "variables": {"ids": ids[i:i + 100]}})
        if "errors" in r and not r.get("data"):
            rej.append({"lot": ids[i], "erreur": str(r["errors"])[:300]}); continue
        for e in r["data"]["entries"] or []:
            if e:
                cand[e["rcsb_id"]] = e
        time.sleep(0.5)
    out = {"coupure": cfg["date"], "operateur": cfg["operateur"], "total_count_recherche": n_ids_total, "n_ids_metadonnees_demandees": len(ids),
           "n_metadonnees": len(cand), "rejets_lots": rej, "entrees": cand}
    f.write_text(json.dumps(out))
    return out


def curer_systeme(C, d, res2):
    d.mkdir(parents=True, exist_ok=True)
    E.ecrire_pdb(res2, d / "brut.pdb")
    rap = {"cif": res2.get("cif"), "chaines": [{"auth": c["auth"], "n_residus": c["n_residus"], "sequence": c["sequence"]} for c in res2["chaines"]],
           "info": res2["info"], "sha256_brut": E.sha256(d / "brut.pdb"), "table_ordinal": E.table_ordinal(res2)}
    (d / "extraction.json").write_text(json.dumps(rap, indent=1, ensure_ascii=False) + "\n")
    cur = C.curer(str(d / "brut.pdb"), sortie_x0=str(d / "x0.pdb"), ecrire_diag=str(d / "x0.diagnostic.pdb"))
    (d / "curation.json").write_text(json.dumps(cur, indent=1, ensure_ascii=False, default=str) + "\n")
    return cur["etat"], cur.get("classes_refus", [])


# ───────────────────────── A_corpus ─────────────────────────
def strate_A_corpus(C):
    man = {x["pdb"]: x for x in (json.loads(l) for l in open(B / "corpus/manifeste_corpus.json") if l.strip())
           if x.get("statut") == "actif_v2" and x["pdb"] not in HORS}
    assert len(man) == 28, len(man)
    cro = json.load(open(B / "generation/cristal/croisement_cristal_psf.json"))
    dates = dates_rcsb(sorted(man))
    out = {}
    for pdb in sorted(man):
        cif = B / f"generation/cristal/{pdb}.cif"
        cv = {k: v for k, v in cro[pdb]["chaines_cristal_vers_md"].items() if v and v.get("chaine_md_origine")}
        md_ch = [c["id"] for c in man[pdb]["chaines"]]
        cand = {m: [k for k, v in cv.items() if v["chaine_md_origine"] == m] for m in md_ch}
        tous = sorted({c for l in cand.values() for c in l})
        ent = {"pdb": pdb, "md_chaines": md_ch, "date": dates.get(pdb), "sha256_cif": E.sha256(cif), "rcsb_methode": man[pdb]["rcsb"]["methode"], "rcsb_resolution_A": man[pdb]["rcsb"]["resolution_A"]}
        if not tous or any(not cand[m] for m in md_ch):
            ent.update(etat="REFUSE", classes_refus=["CHAINES_NON_APPARIEES_AU_CRISTAL"], choix=None, essais=[]); out[pdb] = ent
            print(f"A_corpus {pdb} REFUSE CHAINES_NON_APPARIEES_AU_CRISTAL", flush=True); continue
        res, _ = E.extraire(cif, tous)
        par = {c["auth"]: c for c in res["chaines"]}
        essais = []
        for ca, cb in itertools.product(cand[md_ch[0]], cand[md_ch[1]]):
            if ca == cb or ca not in par or cb not in par:
                continue
            essais.append({"chaines": [ca, cb], "refus_extraction": sorted(set(par[ca]["refus"]) | set(par[cb]["refus"])), "contacts_5A": E.contacts_residus(par[ca], par[cb]),
                           "n_residus": [par[ca]["n_residus"], par[cb]["n_residus"]], "meme_lettre": [cv[ca]["meme_lettre"], cv[cb]["meme_lettre"]]})
        sans = [e for e in essais if not e["refus_extraction"]]
        pool = sans if sans else essais
        best = max(pool, key=lambda e: (e["contacts_5A"], sum(e["meme_lettre"]))) if pool else None
        ent.update(essais=essais, choix=best)
        if best is None:
            ent.update(etat="REFUSE", classes_refus=["AUCUNE_PAIRE_DE_CHAINES"])
        elif best["refus_extraction"]:
            ent.update(etat="REFUSE", classes_refus=best["refus_extraction"])
        else:
            r2, _ = E.extraire(cif, best["chaines"])
            r2["cif"] = str(cif)
            etat, cl = curer_systeme(C, R / "strate_A_corpus" / pdb, r2)
            ent.update(etat=etat, classes_refus=cl, sequences=[c["sequence"] for c in r2["chaines"]], tokens=sum(c["n_residus"] for c in r2["chaines"]))
        out[pdb] = ent
        print(f"A_corpus {pdb} {ent['etat']} {ent['classes_refus']} choix={best['chaines'] if best else None} contacts={best['contacts_5A'] if best else None} "
              f"depot={ent['date']['deposit_date'][:10] if ent['date'] else None}", flush=True)
    (R / "strate_A_corpus").mkdir(exist_ok=True)
    (R / "strate_A_corpus/manifeste.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    c = collections.Counter(v["etat"] for v in out.values())
    cl = collections.Counter(k for v in out.values() for k in v["classes_refus"])
    apres = [p for p, v in out.items() if v["date"] and v["date"]["deposit_date"] > "2021-09-30T23:59:59Z"]
    print(f"BILAN_A_corpus : {dict(c)} · classes de refus {dict(cl)} · deposes APRES 2021-09-30 : {apres} · REJETS = {c.get('REFUSE', 0)}", flush=True)
    return out


# ───────────────────────── A_apparie et B (meme procedure) ─────────────────────────
def periode(dep, cfg):
    for lim, nom in cfg["bornes"]:
        if dep < lim:
            return nom
    return cfg["dernier"]


def strate_rcsb(C, nom):
    cfg = STRATES_RCSB[nom]
    cj = candidats(nom, cfg)
    E_ = cj["entrees"]
    std = set("ACDEFGHIKLMNPQRSTVWY")
    ecartes_pre, eligibles = collections.Counter(), []
    for pid, e in E_.items():
        pe = e["polymer_entities"]
        if len(pe) != 2:
            ecartes_pre["pas_2_entites"] += 1; continue
        seqs = [x["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "") for x in pe]
        strands = [x["entity_poly"]["pdbx_strand_id"] for x in pe]
        if any("," in s for s in strands):
            ecartes_pre["entite_a_plusieurs_chaines"] += 1; continue
        if any(not set(s) <= std for s in seqs):
            ecartes_pre["sequence_non_standard"] += 1; continue
        if min(len(s) for s in seqs) < 30:
            ecartes_pre["chaine_de_moins_de_30_residus"] += 1; continue
        ei = e["rcsb_entry_info"]
        if ei["deposited_modeled_polymer_monomer_count"] < 0.95 * ei["deposited_polymer_monomer_count"]:
            ecartes_pre["moins_de_95pc_modelise"] += 1; continue
        g30 = []
        for x in pe:
            g = [m["group_id"] for m in (x["rcsb_polymer_entity_group_membership"] or []) if m["aggregation_method"] == "sequence_identity" and m["similarity_cutoff"] == 30.0]
            g30.append(g[0] if g else f"sans_groupe_{x['rcsb_id']}")
        dep = e["rcsb_accession_info"]["deposit_date"][:10]
        eligibles.append({"pdb": pid, "depot": dep, "publication": e["rcsb_accession_info"]["initial_release_date"][:10], "periode": periode(dep, cfg),
                          "resolution_A": ei["resolution_combined"][0], "titre": e["struct"]["title"], "chaines_auth": strands, "sequences_entite": seqs, "groupes_30pc": g30,
                          "description": [x["rcsb_polymer_entity"]["pdbx_description"] for x in pe],
                          "organismes": [[o["ncbi_scientific_name"] for o in (x["rcsb_entity_source_organism"] or [])][:2] for x in pe]})
    rnd = random.Random(GRAINE)
    par_per = collections.defaultdict(list)
    for x in eligibles:
        par_per[x["periode"]].append(x)
    for l in par_per.values():
        rnd.shuffle(l)
    ordre, groupes_pris, saut_famille = [], set(), 0
    pers = sorted(par_per)
    while len(ordre) < cfg["max_traites"] and any(par_per[p] for p in pers):
        for p in pers:
            while par_per[p]:
                x = par_per[p].pop()
                if any(g in groupes_pris for g in x["groupes_30pc"]):
                    saut_famille += 1; continue
                groupes_pris.update(x["groupes_30pc"]); ordre.append(x); break
    print(f"{nom} : {len(E_)} candidats avec metadonnees (recherche RCSB : {cj.get('total_count_recherche', cj.get('n_ids_recherche'))}) ; ecartes avant telechargement {dict(ecartes_pre)} ; "
          f"eligibles {len(eligibles)} ; famille deja prise {saut_famille} ; ordre de traitement {len(ordre)} (par periode {dict(collections.Counter(x['periode'] for x in ordre))})", flush=True)
    (R / f"strate_{nom}/cif").mkdir(parents=True, exist_ok=True)
    out, n_pret = {}, 0
    for x in ordre:
        if n_pret >= cfg["cible"]:
            break
        pid = x["pdb"]
        ent = dict(x)
        cif = R / f"strate_{nom}/cif/{pid}.cif"
        try:
            if not cif.is_file():
                with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.cif", timeout=90) as r:
                    cif.write_bytes(r.read())
                time.sleep(1.0)
            ent["sha256_cif"] = E.sha256(cif)
            res, refus = E.extraire(cif, x["chaines_auth"])
        except Exception as ex:
            ent.update(etat="REFUSE", classes_refus=["TELECHARGEMENT_OU_LECTURE"], detail=f"{type(ex).__name__}: {str(ex)[:200]}")
            out[pid] = ent; print(f"{nom} {pid} REFUSE TELECHARGEMENT_OU_LECTURE {ent['detail']}", flush=True); continue
        if len(res["chaines"]) != 2 or refus:
            cl = sorted({r["classe"] for r in refus})
            ent.update(etat="REFUSE", classes_refus=cl, detail=[r["detail"] for r in refus][:4])
            out[pid] = ent; print(f"{nom} {pid} REFUSE {cl}", flush=True); continue
        nct = E.contacts_residus(res["chaines"][0], res["chaines"][1])
        tok = sum(c["n_residus"] for c in res["chaines"])
        ent.update(contacts_5A=nct, tokens=tok, sequences=[c["sequence"] for c in res["chaines"]])
        if nct < SEUIL_CONTACTS:
            ent.update(etat="REFUSE", classes_refus=["INTERFACE_TROP_PETITE"], detail=f"{nct} couples de residus < {SEUIL_CONTACTS}")
            out[pid] = ent; print(f"{nom} {pid} REFUSE INTERFACE_TROP_PETITE {nct}", flush=True); continue
        if tok > 835:
            ent.update(etat="REFUSE", classes_refus=["AU_DELA_DU_PLAFOND_835"], detail=f"{tok} residus")
            out[pid] = ent; print(f"{nom} {pid} REFUSE >835 {tok}", flush=True); continue
        res["cif"] = str(cif)
        etat, cl = curer_systeme(C, R / f"strate_{nom}" / pid, res)
        ent.update(etat=etat, classes_refus=cl)
        out[pid] = ent
        n_pret += etat == "PRET"
        print(f"{nom} {pid} {etat} {cl} contacts={nct} tokens={tok} depot={x['depot']} {x['periode']} pret={n_pret}", flush=True)
    (R / f"strate_{nom}/manifeste.json").write_text(json.dumps({"graine": GRAINE, "operateur_date": cfg["operateur"], "date_coupure": cfg["date"], "ecartes_avant_telechargement": dict(ecartes_pre),
                                                                  "n_eligibles": len(eligibles), "n_traites": len(out), "systemes": out}, indent=1, ensure_ascii=False) + "\n")
    c = collections.Counter(v["etat"] for v in out.values())
    cl = collections.Counter(k for v in out.values() for k in v["classes_refus"])
    print(f"BILAN_{nom} : traites {len(out)} · {dict(c)} · classes de refus {dict(cl)} · REJETS = {c.get('REFUSE', 0)} · cible {cfg['cible']} {'ATTEINTE' if n_pret >= cfg['cible'] else 'NON ATTEINTE'}", flush=True)
    return out


if __name__ == "__main__":
    C = charger_curer()
    quoi = sys.argv[1:] or ["A_corpus", "B", "A_apparie"]
    bilan = {}
    if "A_corpus" in quoi:
        bilan["A_corpus"] = strate_A_corpus(C)
    for nom in ("B", "A_apparie"):
        if nom in quoi:
            bilan[nom] = strate_rcsb(C, nom)
    print("FIN : " + " · ".join(f"{k} pret {sum(v['etat'] == 'PRET' for v in d.values())} sur {len(d)} traites" for k, d in bilan.items()), flush=True)
    sys.exit(0)
