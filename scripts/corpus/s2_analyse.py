"""J2 — inventaire des 210 projets DynaBench, criteres d'inclusion, selection.

Entrees (aucune requete reseau) :
  L0a_work/all_projects.json            cache de l'API (liste des projets) — lu, jamais reinterroge
  L0a_work/dynabench_index.tsv          jointure DB5.5 deja faite en L0a (classe, I-RMSD, complexe)
  L0a_work/db55/db55_classes.json       table BM5.5 parsee (categorie fonctionnelle, noms des partenaires)
  cache/fetch_structures.jsonl          resultat de J1 (structure.pdb replique 1)
  cache/rcsb_graphql.json               reponse RCSB de J1

Criteres (chiffres du brief) :
  C1  exactement 2 chaines proteiques                    (structure.pdb, mdtraj)
  C2  <= 1030 residus au total                           (structure.pdb, mdtraj)
  C3  >= 25 paires de residus a 5,0 A, atomes lourds     (structure.pdb, scipy cKDTree)
  C4  0 residu NON OBSERVE interne a la chaine MD        (API RCSB uniquement : UNOBSERVED_RESIDUE_XYZ)
  C5  resolution declaree                                (API RCSB)

Un rejet est une ENTREE. Le compte par motif est imprime trois fois : tous motifs, premier motif, entonnoir.
"""
import collections
import csv
import json
import os
import re
import sys

import mdtraj as md
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import (CACHE, CONTACT_A, L0A, MAX_RES, MIN_PAIRS, ROOT, STAGING12,  # noqa: E402
                        charger_projets_dynabench, jdump, jload, now_iso)

CAT_LABEL = {  # codes MESURES dans la feuille ; les libelles sont une HYPOTHESE (legende absente du cache)
    "AA": "Antibody-Antigen (Fab/IgG)", "AS": "Antibody-Antigen, single-domain (nanobody / VNAR)",
    "EI": "Enzyme-Inhibitor", "ER": "Enzyme-Regulator/Activator", "ES": "Enzyme-Substrate",
    "OG": "Other, G-protein", "OR": "Other, Receptor", "OX": "Other",
}
CAMELIDS = {"Vicugna pacos", "Lama glama", "Camelus dromedarius", "Camelus bactrianus",
            "Vicugna vicugna", "Lama guanicoe", "Camelidae", "Lama pacos"}
KW_VHH = re.compile(r"nanobod|VHH|camelid|single[- ]domain|sybody|\bcAb\b|llama|alpaca|dromedar|\bNb[ .\-_]?\w*", re.I)
KW_VNAR = re.compile(r"shark|IgNAR|VNAR|(new|novel) antigen receptor", re.I)
FR4_FIN = re.compile(r"VTVS[SA]")   # fin de FR4 d'un domaine V (VHH/VH), cherchee dans les 30 derniers residus


def lire_index():
    rows = {}
    with open(L0A.replace("/L0a_work", "") + "/dynabench_index.tsv", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[r["accession"]] = r
    return rows


def appariement_db55(idx_row, db55):
    """Cherche l'entree de la table BM5.5 qui correspond au complexe donne par l'index L0a."""
    cx = (idx_row.get("db55_complexe") or "").strip()
    pdb = idx_row["pdb"]
    for e in db55.get(pdb, []):
        if e["complex"].strip() == cx:
            return e
    return None


def analyser_structure(path):
    t = md.load(path)
    top = t.topology
    xyz = t.xyz[0] * 10.0  # nm -> A
    elem = [(a.element.symbol if a.element is not None else None) for a in top.atoms]
    heavy = np.array([e not in ("H", "D", None) for e in elem])
    res_of = np.array([a.residue.index for a in top.atoms])
    chains = []
    for ch in top.chains:
        res = list(ch.residues)
        atoms = [a.index for a in ch.atoms]
        prot = [r for r in res if r.is_protein]
        ca = [a.index for a in ch.atoms if a.name == "CA" and a.residue.is_protein]
        d = np.linalg.norm(np.diff(xyz[ca], axis=0), axis=1) if len(ca) > 1 else np.array([])
        chains.append({
            "idx": ch.index, "id": getattr(ch, "chain_id", None),
            "n_res": len(res), "n_res_prot": len(prot),
            "resseq": [r.resSeq for r in res], "resnames": [r.name for r in res],
            "resname_first": res[0].name, "resname_last": res[-1].name,
            "n_atomes": len(atoms), "n_lourds": int(heavy[atoms].sum()),
            "n_ca": len(ca), "ca_ruptures_md": int((d > 4.5).sum()), "ca_dmax": float(d.max()) if d.size else None,
            "_atomes_lourds": [i for i in atoms if heavy[i]],
        })
    prot_chains = [c for c in chains if c["n_res_prot"] > 0 and c["n_res_prot"] >= 0.9 * c["n_res"]]
    nonprot_res = sum(c["n_res"] - c["n_res_prot"] for c in chains)
    # paires de residus a 5 A entre chaines proteiques, toutes combinaisons
    paires = {}
    for i in range(len(prot_chains)):
        for j in range(i + 1, len(prot_chains)):
            a, b = prot_chains[i], prot_chains[j]
            ta, tb = cKDTree(xyz[a["_atomes_lourds"]]), cKDTree(xyz[b["_atomes_lourds"]])
            lst = ta.query_ball_tree(tb, CONTACT_A + 1e-6)
            s = set()
            for ia, jb in enumerate(lst):
                for k in jb:
                    s.add((res_of[a["_atomes_lourds"][ia]], res_of[b["_atomes_lourds"][k]]))
            paires["%s-%s" % (a["id"], b["id"])] = len(s)
    for c in chains:
        c.pop("_atomes_lourds")
    return {"n_atomes": top.n_atoms, "n_residus": top.n_residues, "n_chaines": len(chains),
            "n_elements_inconnus": int(sum(1 for e in elem if e is None)),
            "chaines": chains, "n_chaines_proteiques": len(prot_chains),
            "n_residus_proteine": sum(c["n_res_prot"] for c in prot_chains),
            "n_residus_non_proteine": nonprot_res, "paires_5A": paires}


def index_rcsb(rc):
    ent, ko = {}, []
    for ch in rc["chunks"]:
        if ch["reponse"] and ch["reponse"].get("data") and ch["reponse"]["data"].get("entries"):
            for e in ch["reponse"]["data"]["entries"]:
                if e:
                    ent[e["rcsb_id"]] = e
        else:
            ko.append({"ids": ch["ids_demandes"], "erreur": ch["erreur"], "http": ch["http"]})
    return ent, ko


AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
       "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K"}


def seq_md(chain):
    return "".join(AA3.get(n, "X") for n in chain["resnames"])


def _aligneur():
    from Bio import Align
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 2.0, -1.0
    al.open_gap_score, al.extend_gap_score = -8.0, -0.5
    al.end_deletion_score = 0.0     # cible = entite RCSB : ses extremites peuvent deborder la chaine MD (teste : score 40 contre -75)
    return al


_AL = None


def positions_par_numerotation(chain, sq_md, seq_ent, mp):
    pos_of = {}
    for i, v in enumerate(mp):
        if v not in (None, "", "?") and v not in pos_of:
            pos_of[v] = i + 1
    pos, ident = [], 0
    for rs, aa in zip(chain["resseq"], sq_md):
        p = pos_of.get(str(rs))
        pos.append(p)
        if p and seq_ent[p - 1] == aa:
            ident += 1
    mapped = [p for p in pos if p]
    mono = all(b > a for a, b in zip(mapped, mapped[1:]))
    return pos, ident / len(pos), len(mapped) / len(pos), mono


def positions_par_alignement(sq_md, seq_ent):
    global _AL
    if _AL is None:
        _AL = _aligneur()
    aln = _AL.align(seq_ent, sq_md)[0]
    pos = [None] * len(sq_md)
    ident = 0
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        t0, t1, q0, q1 = int(t0), int(t1), int(q0), int(q1)
        for k in range(t1 - t0):
            pos[q0 + k] = t0 + k + 1
            if seq_ent[t0 + k] == sq_md[q0 + k]:
                ident += 1
    mapped = sum(1 for p in pos if p)
    return pos, ident / len(sq_md), mapped / len(sq_md)


def non_observes(inst):
    u = set()
    for f in inst["rcsb_polymer_instance_feature"] or []:
        if f["type"] == "UNOBSERVED_RESIDUE_XYZ":
            for fp in f["feature_positions"] or []:
                u.update(range(fp["beg_seq_id"], fp["end_seq_id"] + 1))
    return u


def analyser_chaine_rcsb(entry, chain):
    """Apparie la chaine MD a une instance RCSB (par NUMEROTATION auteur verifiee par la sequence, sinon par
    ALIGNEMENT), puis lit les residus NON OBSERVES internes dans l'API RCSB. La verite des gaps vient de RCSB ;
    l'appariement n'est qu'un moyen de savoir de quelle chaine on parle (les lettres de chaine de DynaBench ne
    sont PAS celles de RCSB : 3K75, 1BKD, 2AJF, 2I25, 5HGG).
    """
    out = {"chaine": chain["id"], "evaluable": False}
    sq = seq_md(chain)
    cands = []
    for pe in entry["polymer_entities"] or []:
        if (pe["entity_poly"] or {}).get("rcsb_entity_polymer_type") != "Protein":
            continue
        seq_ent = pe["entity_poly"]["pdbx_seq_one_letter_code_can"] or ""
        for inst in pe["polymer_entity_instances"] or []:
            ids = inst["rcsb_polymer_entity_instance_container_identifiers"]
            pos, ident, mapped, mono = positions_par_numerotation(chain, sq, seq_ent, ids["auth_to_entity_poly_seq_mapping"] or [])
            ok = ident >= 0.9 and mapped >= 0.9 and mono
            cands.append({"pe": pe, "inst": inst, "ids": ids, "seq": seq_ent, "methode": "numerotation" if ok else None,
                          "pos": pos if ok else None, "ident": ident, "mapped": mapped, "ok": ok,
                          "ident_num": ident, "mapped_num": mapped, "mono": mono})
    if not any(c["ok"] for c in cands):
        done = {}
        for c in cands:
            key = c["ids"]["entity_id"]
            if key not in done:
                pos, ident, mapped = positions_par_alignement(sq, c["seq"])
                done[key] = (pos, ident, mapped)
            pos, ident, mapped = done[key]
            if ident >= 0.9 and mapped >= 0.9:
                c.update({"ok": True, "methode": "alignement", "pos": pos, "ident": ident, "mapped": mapped})
    ok = [c for c in cands if c["ok"]]
    if not ok:
        best = max(cands, key=lambda c: c["ident_num"]) if cands else None
        out["raison"] = ("aucune instance RCSB ne correspond (meilleure identite par numerotation %.2f, residus mappes %.2f)" % (
            best["ident_num"], best["mapped_num"])) if best else "aucune entite proteique dans l'entree RCSB"
        return out
    for c in ok:
        u = non_observes(c["inst"])
        c["u"] = u
        c["conflits"] = sum(1 for p in c["pos"] if p and p in u)
    ok.sort(key=lambda c: (0 if c["methode"] == "numerotation" else 1, -round(c["ident"], 2), c["conflits"],
                           0 if c["ids"]["auth_asym_id"] == chain["id"] else 1))
    c = ok[0]
    pe, ids, u = c["pe"], c["ids"], c["u"]
    md_pos = [p for p in c["pos"] if p]
    org = sorted({n for n in ((o or {}).get("ncbi_scientific_name") for o in (pe.get("rcsb_entity_source_organism") or []) if o) if n})
    lo, hi = min(md_pos), max(md_pos)
    interne = sorted(p for p in u if lo < p < hi)
    mdset = set(md_pos)
    absents_span = [p for p in range(lo, hi + 1) if p not in mdset]
    runs = []
    for p in interne:
        if runs and p == runs[-1][1] + 1:
            runs[-1][1] = p
        else:
            runs.append([p, p])
    inv = {i + 1: v for i, v in enumerate(ids["auth_to_entity_poly_seq_mapping"] or [])}
    out.update({"evaluable": True, "instance_rcsb": c["inst"]["rcsb_id"], "chaine_auteur_rcsb": ids["auth_asym_id"],
                "meme_lettre": ids["auth_asym_id"] == chain["id"], "methode_appariement": c["methode"],
                "identite": round(c["ident"], 4), "fraction_mappee": round(c["mapped"], 4),
                "n_candidats_acceptes": len(ok), "conflits_md_non_observes": c["conflits"],
                "entite": ids["entity_id"], "description": (pe["rcsb_polymer_entity"] or {}).get("pdbx_description"),
                "organismes": org, "longueur_entite": len(c["seq"]), "sequence_entite": c["seq"],
                "span_entite": [lo, hi], "n_non_observes_total": len(u), "n_non_observes_internes": len(interne),
                "runs_internes_entite": runs, "runs_internes_auteur": [[inv.get(a), inv.get(b)] for a, b in runs],
                "internes_presents_dans_MD": sum(1 for p in interne if p in mdset),
                "residus_absents_du_span_MD": len(absents_span)})
    return out


def vhh_chaine(info, cat):
    seq = info.get("sequence_entite") or ""
    desc = info.get("description") or ""
    org = set(info.get("organismes") or [])
    sig = {"fr4": bool(FR4_FIN.search(seq[-30:])), "longueur_v": 0 < len(seq) <= 160,
           "camelide": bool(org & CAMELIDS), "mot_cle": bool(KW_VHH.search(desc)),
           "vnar": bool(KW_VNAR.search(desc)) or any("Ginglymostoma" in o or "Orectolobus" in o for o in org),
           "cat_AS": cat == "AS"}
    sig["V_like"] = sig["fr4"] and sig["longueur_v"]
    sig["vhh"] = bool(sig["V_like"] and not sig["vnar"] and (sig["camelide"] or sig["mot_cle"] or sig["cat_AS"]))
    conv = bool(org & {"Homo sapiens", "Mus musculus", "Rattus norvegicus", "Oryctolagus cuniculus"})
    sig["vhh_probable_non_confirme"] = bool(sig["V_like"] and not sig["vnar"] and not sig["vhh"] and not conv)
    return sig


def main():
    ld = charger_projets_dynabench()
    projs = ld["projets"]
    idx = lire_index()
    db55 = jload(L0A + "/db55/db55_classes.json")
    fetch = {}
    for line in open(CACHE + "/fetch_structures.jsonl", encoding="utf-8"):
        r = json.loads(line)
        fetch[r["accession"]] = r
    rc = jload(CACHE + "/rcsb_graphql.json")
    rcent, rcko = index_rcsb(rc)
    print("ENTREES projets=%d cache_sha256=%s  fetch=%d  rcsb_entrees=%d  rcsb_paquets_en_echec=%d" % (
        len(projs), ld["sha256"], len(fetch), len(rcent), len(rcko)))

    inv = []
    for q in projs:
        acc, m = q["accession"], q["metadata"]
        pdb = m["PDBIDS"][0]
        ir = idx[acc]
        e55 = appariement_db55(ir, db55)
        P = {"accession": acc, "pdb": pdb, "nom": m.get("NAME"), "mds": q.get("mds"),
             "n_repliques": len(q.get("mds") or []), "protres_api": m.get("PROTRES"), "mdAtoms_api": m.get("mdAtoms"),
             "dans_staging12": pdb in STAGING12, "motifs_rejet": [], "notes": []}
        if e55:
            P.update({"classe_difficulte": e55["cls"], "irmsd_db55": e55["irmsd"], "complexe_db55": e55["complex"].strip(),
                      "categorie_code": e55["cat"], "categorie_label_hypothese": CAT_LABEL.get(e55["cat"]),
                      "partenaire1": e55["p1"], "partenaire2": e55["p2"]})
        else:
            P.update({"classe_difficulte": None, "irmsd_db55": None, "complexe_db55": None, "categorie_code": None,
                      "categorie_label_hypothese": None, "partenaire1": None, "partenaire2": None})
            P["notes"].append("absent de la table BM5.5 sous ce code (index L0a: %s)" % ir["db55_classe"])
        # ---- structure
        fr = fetch.get(acc)
        if not fr or not fr.get("conforme"):
            P["structure"] = None
            P["motifs_rejet"].append("structure_indisponible")
        else:
            try:
                P["structure"] = analyser_structure(fr["dest"])
                P["structure_sha256"] = fr["sha256"]
                P["structure_octets"] = fr["octets"]
            except Exception as ex:  # noqa: BLE001
                P["structure"] = None
                P["motifs_rejet"].append("structure_illisible")
                P["notes"].append("mdtraj: %s: %s" % (type(ex).__name__, ex))
        S = P["structure"]
        if S:
            P["chaines"] = [{"id": c["id"], "n_res": c["n_res"], "n_res_prot": c["n_res_prot"],
                             "resseq_premier": c["resseq"][0], "resseq_dernier": c["resseq"][-1],
                             "ca_ruptures_md": c["ca_ruptures_md"],
                             "ca_dmax_A": round(c["ca_dmax"], 3) if c["ca_dmax"] is not None else None,
                             "sauts_numerotation": sum(1 for a, b in zip(c["resseq"], c["resseq"][1:]) if b - a != 1)}
                            for c in S["chaines"]]
            P["n_residus_total"] = S["n_residus_proteine"]
            P["coherence_PROTRES"] = (S["n_residus_proteine"] == m.get("PROTRES"))
            P["coherence_mdAtoms"] = (S["n_atomes"] == m.get("mdAtoms"))
            if not P["coherence_PROTRES"]:
                P["notes"].append("PROTRES API=%s mais mdtraj=%d residus proteiques" % (m.get("PROTRES"), S["n_residus_proteine"]))
            P["paires_5A"] = S["paires_5A"]
            if S["n_chaines_proteiques"] == 2:
                P["interface_paires_5A"] = list(S["paires_5A"].values())[0]
            else:
                P["interface_paires_5A"] = None
        # ---- RCSB
        E = rcent.get(pdb)
        P["rcsb"] = None
        if E is None:
            P["motifs_rejet"].append("rcsb_indisponible")
        else:
            res = (E["rcsb_entry_info"] or {}).get("resolution_combined") or []
            P["rcsb"] = {"titre": (E["struct"] or {}).get("title"),
                         "methode": ((E["exptl"] or [{}])[0] or {}).get("method"),
                         "resolution_A": res[0] if res else None,
                         "ls_d_res_high": (((E.get("refine") or [{}]) or [{}])[0] or {}).get("ls_d_res_high")}
            if S:
                P["rcsb"]["chaines"] = [analyser_chaine_rcsb(E, c) for c in S["chaines"]
                                        if c["n_res_prot"] > 0]
                P["gaps_internes_total"] = sum(c.get("n_non_observes_internes", 0) for c in P["rcsb"]["chaines"])
                P["gap_evaluable"] = all(c["evaluable"] for c in P["rcsb"]["chaines"])
        # ---- VHH
        cat = P["categorie_code"]
        P["vhh_chaines"] = []
        if P["rcsb"] and P["rcsb"].get("chaines"):
            for c in P["rcsb"]["chaines"]:
                if c["evaluable"]:
                    sg = vhh_chaine(c, cat)
                    sg["chaine"] = c["chaine"]
                    sg["description"] = c.get("description")
                    sg["organismes"] = c.get("organismes")
                    P["vhh_chaines"].append(sg)
        P["vhh"] = any(s["vhh"] for s in P["vhh_chaines"])
        P["vhh_probable_non_confirme"] = (not P["vhh"]) and any(s["vhh_probable_non_confirme"] for s in P["vhh_chaines"])
        P["vnar"] = any(s["vnar"] for s in P["vhh_chaines"]) or bool(KW_VNAR.search("%s %s" % (P["partenaire1"], P["partenaire2"])))
        inv.append(P)

    # ---------------- criteres ----------------
    for P in inv:
        S = P["structure"]
        mot = P["motifs_rejet"]
        if S:
            if S["n_chaines_proteiques"] != 2:
                mot.append("chaines_proteiques_!=2")
            if P["n_residus_total"] > MAX_RES:
                mot.append("taille_>1030")
            if S["n_chaines_proteiques"] == 2 and P["interface_paires_5A"] < MIN_PAIRS:
                mot.append("interface_<25")
        if P["rcsb"]:
            if P["rcsb"]["resolution_A"] is None:
                mot.append("resolution_non_declaree")
            if S:
                if not P.get("gap_evaluable"):
                    mot.append("gap_non_evaluable")
                elif P["gaps_internes_total"] > 0:
                    mot.append("gap_interne_rcsb")
        P["passe_les_criteres"] = len(mot) == 0

    # ---------------- comptes ----------------
    ordre = ["structure_indisponible", "structure_illisible", "rcsb_indisponible", "chaines_proteiques_!=2",
             "taille_>1030", "interface_<25", "resolution_non_declaree", "gap_non_evaluable", "gap_interne_rcsb"]
    tous = collections.Counter(m for P in inv for m in P["motifs_rejet"])
    premier = collections.Counter()
    for P in inv:
        for m in ordre:
            if m in P["motifs_rejet"]:
                premier[m] += 1
                break
    n_pass = sum(1 for P in inv if P["passe_les_criteres"])
    print("COMPTE total=%d  passent_les_5_criteres=%d  rejetes=%d" % (len(inv), n_pass, len(inv) - n_pass))
    for m in ordre:
        print("MOTIF_TOUS   %-26s %3d" % (m, tous.get(m, 0)))
    for m in ordre:
        print("MOTIF_PREMIER %-26s %3d" % (m, premier.get(m, 0)))
    print("MOTIF_PREMIER somme=%d (doit egaler rejetes=%d)" % (sum(premier.values()), len(inv) - n_pass))
    # entonnoir sequentiel
    surv = [P for P in inv]
    print("ENTONNOIR depart %d" % len(surv))
    for nom, test in [
        ("structure_lisible", lambda P: P["structure"] is not None),
        ("2_chaines_proteiques", lambda P: P["structure"]["n_chaines_proteiques"] == 2),
        ("taille<=1030", lambda P: P["n_residus_total"] <= MAX_RES),
        ("interface>=25", lambda P: P["interface_paires_5A"] >= MIN_PAIRS),
        ("resolution_declaree", lambda P: P["rcsb"] is not None and P["rcsb"]["resolution_A"] is not None),
        ("0_gap_interne_RCSB", lambda P: P.get("gap_evaluable") and P["gaps_internes_total"] == 0)]:
        surv = [P for P in surv if test(P)]
        print("ENTONNOIR apres %-22s %3d" % (nom, len(surv)))
    # sensibilite (separee, ne modifie rien)
    for g in (0, 1, 2, 3, 5, 10):
        n = sum(1 for P in inv if P["structure"] and P["structure"]["n_chaines_proteiques"] == 2
                and P["n_residus_total"] <= MAX_RES and P["interface_paires_5A"] >= MIN_PAIRS
                and P["rcsb"] and P["rcsb"]["resolution_A"] is not None and P.get("gap_evaluable")
                and P["gaps_internes_total"] <= g)
        print("SENSIBILITE tolerance_gap<=%2d residu(s) -> %d systemes" % (g, n))
    for s in (15, 20, 25):
        n = sum(1 for P in inv if P["structure"] and P["structure"]["n_chaines_proteiques"] == 2
                and P["n_residus_total"] <= MAX_RES and P["interface_paires_5A"] >= s
                and P["rcsb"] and P["rcsb"]["resolution_A"] is not None and P.get("gap_evaluable")
                and P["gaps_internes_total"] == 0)
        print("SENSIBILITE interface>=%2d paires (gap strict) -> %d systemes" % (s, n))

    # ---------------- selection ----------------
    CIBLE = 30
    pas = [P for P in inv if P["passe_les_criteres"]]
    ecartes_sel = [P for P in pas if P["categorie_code"] is None]
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
    for P in inv:
        P["selectionne"] = P["accession"] in raison
        P["raison_selection"] = raison.get(P["accession"])
        P["ecarte_de_la_selection"] = ("categorie_inconnue" if P in ecartes_sel else None)
    print("SELECTION n=%d  (deja_staging=%d, vhh=%d, plus_petits=%d)  cible=%d  passent=%d  eligibles_avec_categorie=%d  ecartes_categorie_inconnue=%d" % (
        len(sel), sum(1 for v in raison.values() if v == "deja_sur_staging"), sum(1 for v in raison.values() if v == "vhh"),
        sum(1 for v in raison.values() if v == "plus_petits"), CIBLE, len(pas), len(elig), len(ecartes_sel)))

    # ---------------- sorties ----------------
    for P in inv:
        if P["structure"]:
            P["structure"] = {k: v for k, v in P["structure"].items() if k != "chaines"}
            # les chaines detaillees sont dans P["chaines"]; on retire la liste complete des resSeq (volumineuse)
    jdump({"meta": {"date": now_iso(), "cache_projets": ld["fichier"], "cache_projets_sha256": ld["sha256"],
                    "n_projets": len(inv), "n_passent": n_pass, "criteres": {"chaines_proteiques": 2, "max_residus": MAX_RES,
                    "min_paires_5A": MIN_PAIRS, "contact_A": CONTACT_A, "gap_interne": "RCSB UNOBSERVED_RESIDUE_XYZ",
                    "resolution": "declaree"}, "rcsb_paquets_en_echec": rcko},
           "projets": inv}, ROOT + "/inventaire_dynabench.json")
    jdump({"cible": CIBLE, "selection": [{"accession": P["accession"], "pdb": P["pdb"], "raison": raison[P["accession"]],
                                           "n_residus_total": P["n_residus_total"]} for P in sel],
           "reserve": [{"accession": P["accession"], "pdb": P["pdb"], "n_residus_total": P["n_residus_total"]} for P in reserve],
           "ecartes_categorie_inconnue": [{"accession": P["accession"], "pdb": P["pdb"], "n_residus_total": P["n_residus_total"]} for P in ecartes_sel]},
          ROOT + "/selection.json")
    jdump({"tous_motifs": dict(tous), "premier_motif": dict(premier), "ordre": ordre,
           "rejets": [{"accession": P["accession"], "pdb": P["pdb"], "motifs": P["motifs_rejet"]} for P in inv if P["motifs_rejet"]]},
          ROOT + "/rejets_par_motif.json")
    cols = ["pdb", "accession", "categorie_code", "classe_difficulte", "irmsd_db55", "chaines", "n_residus_total", "interface_paires_5A",
            "gaps_internes_total", "resolution_A", "methode", "vhh", "vnar", "dans_staging12", "passe", "motifs", "selectionne", "raison"]
    with open(ROOT + "/inventaire_dynabench.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for P in inv:
            ch = ",".join("%s:%d" % (c["id"], c["n_res"]) for c in (P.get("chaines") or []))
            f.write("\t".join(str(x) for x in [
                P["pdb"], P["accession"], P["categorie_code"], P["classe_difficulte"], P["irmsd_db55"], ch, P.get("n_residus_total"),
                P.get("interface_paires_5A"), P.get("gaps_internes_total"), (P["rcsb"] or {}).get("resolution_A"),
                (P["rcsb"] or {}).get("methode"), P["vhh"], P["vnar"], P["dans_staging12"], P["passe_les_criteres"],
                ";".join(P["motifs_rejet"]), P["selectionne"], P["raison_selection"]]) + "\n")

    # ---------------- journal lisible (une ligne par projet, citable) ----------------
    for P in inv:
        ch = ",".join("%s:%d" % (c["id"], c["n_res"]) for c in (P.get("chaines") or []))
        print("INV %s %s cat=%s cls=%s chaines=[%s] total=%s protres_api=%s paires5A=%s gaps=%s res=%s meth=%s vhh=%s motifs=%s sel=%s" % (
            P["pdb"], P["accession"], P["categorie_code"], P["classe_difficulte"], ch, P.get("n_residus_total"), P["protres_api"],
            P.get("interface_paires_5A") if P.get("interface_paires_5A") is not None else (P.get("paires_5A")),
            P.get("gaps_internes_total"), (P["rcsb"] or {}).get("resolution_A"), (P["rcsb"] or {}).get("methode"), P["vhh"],
            ",".join(P["motifs_rejet"]) or "-", P["raison_selection"] or "-"))
    print("--- appariement chaines MD -> RCSB")
    cm = collections.Counter()
    for P in inv:
        for c in ((P["rcsb"] or {}).get("chaines") or []):
            cm[(c.get("methode_appariement") or "AUCUN", bool(c.get("meme_lettre")))] += 1
            if not c["evaluable"]:
                print("NONAPP %s %s chaine=%s raison=%s" % (P["pdb"], P["accession"], c["chaine"], c.get("raison")))
    for k, v in sorted(cm.items(), key=lambda kv: str(kv[0])):
        print("APPARIEMENT methode=%s meme_lettre=%s : %d chaines" % (k[0], k[1], v))
    print("--- coherence PROTRES/mdAtoms")
    nb = [P for P in inv if P.get("coherence_PROTRES") is False or P.get("coherence_mdAtoms") is False]
    print("INCOHERENCES PROTRES/mdAtoms: %d" % len(nb))
    for P in nb:
        print("INCOH %s %s PROTRES_api=%s mdtraj=%s mdAtoms_api=%s mdtraj_atomes=%s" % (
            P["pdb"], P["accession"], P["protres_api"], P.get("n_residus_total"), P["mdAtoms_api"], P["structure"]["n_atomes"] if P["structure"] else None))
    print("--- VHH")
    for P in inv:
        if P["vhh"] or P["vhh_probable_non_confirme"] or P["vnar"] or P["categorie_code"] == "AS":
            print("VHH %s %s cat=%s vhh=%s probable_non_confirme=%s vnar=%s chaines=%s motifs=%s" % (
                P["pdb"], P["accession"], P["categorie_code"], P["vhh"], P["vhh_probable_non_confirme"], P["vnar"],
                json.dumps([(s["chaine"], s["description"], s["organismes"], s["V_like"], s["camelide"], s["mot_cle"]) for s in P["vhh_chaines"]], ensure_ascii=False),
                ",".join(P["motifs_rejet"]) or "-"))
    print("--- SELECTION")
    for P in sel:
        ch = ",".join("%s:%d" % (c["id"], c["n_res"]) for c in P["chaines"])
        print("SEL %s %s %s cat=%s cls=%s total=%d chaines=[%s] paires5A=%d vhh=%s" % (
            P["raison_selection"], P["pdb"], P["accession"], P["categorie_code"], P["classe_difficulte"], P["n_residus_total"], ch,
            P["interface_paires_5A"], P["vhh"]))
    for P in reserve:
        print("RESERVE %s %s cat=%s cls=%s total=%d paires5A=%d" % (P["pdb"], P["accession"], P["categorie_code"], P["classe_difficulte"], P["n_residus_total"], P["interface_paires_5A"]))
    for P in ecartes_sel:
        print("ECARTE_CATEGORIE_INCONNUE %s %s total=%d" % (P["pdb"], P["accession"], P["n_residus_total"]))
    print("FINI s2 %s" % now_iso())


if __name__ == "__main__":
    main()
