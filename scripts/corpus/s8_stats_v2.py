"""Statistiques v2 : on RELIT inventaire_dynabench_v2.json (rien n'est recalcule), sauf la comparaison
`_input.pdb` / `structure.pdb` des 12 de staging, qui appelle les MEMES fonctions que s7 (importees).
Chaque ligne a une etiquette pour pouvoir etre citee `fichier:ligne`.
"""
import collections
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
import preparer_cristal as pc  # noqa: E402
import s2_analyse as s2  # noqa: E402
import s7_commensurabilite as s7  # noqa: E402

ROOT = cl.ROOT
inv = cl.jload(ROOT + "/inventaire_dynabench_v2.json")["projets"]
sel = cl.jload(ROOT + "/selection_v2.json")
vsel = {s["accession"]: s for s in sel["selection"]}
par_pdb = {P["pdb"]: P for P in inv}

print("N projets=%d  passent_v2=%d  selectionnes_v2=%d" % (len(inv), sum(1 for P in inv if P["passe_les_criteres_v2"]), len(vsel)))
print("PREPARATEUR sha256_version_lue=%s  sha256_production=%s" % (cl.sha256_file(ROOT + "/vendor/preparer_cristal.py"), cl.sha256_file(ROOT + "/vendor/preparer_orig.py")))

# ---------- origine (production) contre version lue, sur les structures MD ----------
cause = collections.Counter()
liste_cause = collections.defaultdict(list)
for P in inv:
    if P["commensurable_preparateur_origine"] == (not any(m in P["motifs_rejet"] for m in ("rupture_chaine_MD", "non_commensurable")) and bool(P["commensurabilite"])):
        continue
    c_mod = c_saut = c_refus = False
    for v in P["commensurabilite"].values():
        o = v.get("origine", {})
        if "refus" in o:
            c_refus = True
            continue
        if not o.get("accord_n_res") or not o.get("accord_fasta_sequence"):
            c_mod = True
        if o.get("n_gaps_sauts_resSeq", 0) > 0:
            c_saut = True
    k = "+".join(x for x, f in (("residus_modifies_perdus", c_mod), ("sauts_resSeq_comptes_comme_gaps", c_saut), ("refus", c_refus)) if f)
    cause[k] += 1
    liste_cause[k].append(P["pdb"])
print("ORIGINE_VS_LUE projets_ou_la_production_et_la_version_lue_different=%d  %s" % (sum(cause.values()), json.dumps(dict(cause))))
for k, v in liste_cause.items():
    if k != "sauts_resSeq_comptes_comme_gaps":
        print("ORIGINE_VS_LUE_DETAIL %s : %s" % (k, ",".join(sorted(v))))
for pdb in ("1FQ1", "1ZM4"):
    P = par_pdb[pdb]
    print("ORIGINE_SUR_STRUCTURE_PDB %s origine=%s" % (pdb, json.dumps({k: v.get("origine") for k, v in P["commensurabilite"].items()})))
n_sel_diff = [P["pdb"] for P in inv if P["accession"] in vsel and P["commensurable_preparateur_origine"] is False]
print("ORIGINE_VS_LUE parmi_les_30_retenus_la_production_declarerait_non_commensurables=%d : %s" % (len(n_sel_diff), ",".join(sorted(n_sel_diff))))
mods = [(P["pdb"], k, v["fenetre_2.9"]["n_modifies"]) for P in inv for k, v in P["commensurabilite"].items() if v.get("fenetre_2.9", {}).get("n_modifies")]
print("MODIFIES chaines=%d residus=%d : %s" % (len(mods), sum(m[2] for m in mods), ",".join("%s/%s:%d" % m for m in mods)))

# ---------- `_input.pdb` des 12 contre `structure.pdb` (rep1) ----------
print("--- INPUT.pdb (evomd-data) contre structure.pdb : meme preparateur, deux fenetres")
for pdb in cl.STAGING12:
    P = par_pdb[pdb]
    f_in = "%s/%s/%s_input.pdb" % (cl.STAGING_INPUTS, pdb, pdb)
    f_st = "%s/structure_rep1/%s_%s.pdb" % (cl.CACHE, pdb, P["accession"])
    ti, ts = md.load(f_in), md.load(f_st)
    fi, fs = ti.topology.to_fasta(), ts.topology.to_fasta()
    Si = s2.analyser_structure(f_in)
    prot = [c for c in Si["chaines"] if c["n_res_prot"] > 0 and c["n_res_prot"] >= 0.9 * c["n_res"]]
    fasta_i = {ch.index: fi[k] for k, ch in enumerate(ti.topology.chains)}
    ev = s7.evaluer_preparateur(ti, f_in, prot, fasta_i)
    for k, v in ev.items():
        w29, w19 = v["fenetre_2.9"], v["fenetre_1.9"]
        idx_s = [c for c in ts.topology.chains if pc.chaine_id(ts.topology, c) == k][0].index
        seq_i = fasta_i[[c for c in ti.topology.chains if pc.chaine_id(ti.topology, c) == k][0].index]
        seq_s = fs[idx_s]
        print("INPUT_VS_STRUCTURE %s chaine=%s n_res_input=%d n_res_structure=%d sequence_to_fasta_identique=%s | input: gaps(2.9)=%s gaps(1.9)=%s jonctions_1.9-2.9=%s liaison_max_A=%s | structure.pdb: gaps(2.9)=%s" % (
            pdb, k, v["n_res_mdtraj"], len(list([r for r in ts.topology.chains if pc.chaine_id(ts.topology, r) == k][0].residues)), seq_i == seq_s,
            w29.get("n_gaps"), w19.get("n_gaps"), w29.get("n_jonctions_1.9_2.9"), w29.get("liaison_max_A"),
            P["commensurabilite"][k]["fenetre_2.9"]["n_gaps"]))
        if w19.get("n_gaps") or w29.get("n_gaps"):
            print("INPUT_GAPS %s chaine=%s fenetre_2.9=%s fenetre_1.9=%s" % (pdb, k, json.dumps(w29.get("gaps")), json.dumps(w19.get("gaps"))))

# ---------- reconstruction : le corpus v2 ----------
def agreg(filtre, etiq):
    nc = ni = nin = npres = nabs = ncrec = tn = tc = tnp = tcp = nct = 0
    for P in inv:
        if not filtre(P):
            continue
        for c in P["rcsb_v2"]["chaines"]:
            if not c.get("evaluable"):
                continue
            nc += 1
            I, N, C_ = c["internes"], c["n_terminaux"], c["c_terminaux"]
            if I["n"]:
                ni += 1
                nin += I["n"]
                npres += I["presents_dans_MD"]
                nabs += I["absents_de_MD"]
                ncrec += bool(I["presents_dans_MD"])
            tn += N["n"]
            tc += C_["n"]
            tnp += N["presents_dans_MD"]
            tcp += C_["presents_dans_MD"]
            nct += bool(N["n"] or C_["n"])
    print("RECONSTRUCTION[%s] chaines=%d  avec_residus_internes_non_observes=%d  residus_internes_non_observes=%d  reconstruits(presents_dans_la_MD)=%d  absents_de_la_MD=%d  chaines_avec_au_moins_un_reconstruit=%d" % (etiq, nc, ni, nin, npres, nabs, ncrec))
    print("TERMINAUX[%s] chaines_avec_residus_terminaux_non_observes=%d  N-terminaux=%d (dans_la_MD=%d)  C-terminaux=%d (dans_la_MD=%d)" % (etiq, nct, tn, tnp, tc, tcp))


agreg(lambda P: True, "210")
agreg(lambda P: P["accession"] in vsel, "corpus_v2_30")
agreg(lambda P: P["accession"] in vsel and P["statut_vs_v1"] == "reste", "corpus_v2_restent_19")
agreg(lambda P: P["accession"] in vsel and P["statut_vs_v1"] == "entre", "corpus_v2_entrants_11")
for s in sel["selection"]:
    P = par_pdb[s["pdb"]]
    det = []
    for c in P["rcsb_v2"]["chaines"]:
        if c.get("evaluable") and (c["internes"]["n"] or c["n_terminaux"]["n"] or c["c_terminaux"]["n"]):
            det.append("%s:internes=%d(reconstruits=%d,absents=%d)/N-term=%d/C-term=%d" % (c["chaine"], c["internes"]["n"], c["internes"]["presents_dans_MD"], c["internes"]["absents_de_MD"], c["n_terminaux"]["n"], c["c_terminaux"]["n"]))
    print("CORPUS_V2 %-5s %s statut=%-6s raison=%-16s total=%4d  %s" % (P["pdb"], P["accession"], P["statut_vs_v1"], s["raison"], P["n_residus_total"], " ; ".join(det) or "aucun residu non observe"))
n_corpus_rec = sum(1 for s in sel["selection"] for c in par_pdb[s["pdb"]]["rcsb_v2"]["chaines"] if c.get("evaluable") and c["internes"]["presents_dans_MD"])
print("CORPUS_V2 chaines_avec_residus_reconstruits=%d  systemes_concernes=%d sur %d" % (
    n_corpus_rec, sum(1 for s in sel["selection"] if any(c.get("evaluable") and c["internes"]["presents_dans_MD"] for c in par_pdb[s["pdb"]]["rcsb_v2"]["chaines"])), len(sel["selection"])))
print("CORPUS_V2 systemes_avec_residus_internes_absents_de_la_MD=%s" % ",".join(
    s["pdb"] for s in sel["selection"] if any(c.get("evaluable") and c["internes"]["absents_de_MD"] for c in par_pdb[s["pdb"]]["rcsb_v2"]["chaines"])))
print("FINI s8")
