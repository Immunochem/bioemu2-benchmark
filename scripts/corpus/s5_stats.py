"""Statistiques de l'inventaire, lues dans inventaire_dynabench.json (aucun recalcul : on RELIT ce que s2 a ecrit).
Chaque ligne commence par une etiquette pour pouvoir etre citee `fichier:ligne`.
"""
import collections
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import CACHE, ROOT, STAGING12, jload  # noqa: E402

inv = jload(ROOT + "/inventaire_dynabench.json")["projets"]
sel = jload(ROOT + "/selection.json")
vsel = {P["accession"] for P in inv if P["selectionne"]}

print("N projets=%d  passent=%d  selectionnes=%d" % (len(inv), sum(1 for P in inv if P["passe_les_criteres"]), len(vsel)))
print("HIST_CHAINES_PROTEIQUES " + json.dumps(dict(sorted(collections.Counter(P["structure"]["n_chaines_proteiques"] for P in inv).items()))))
print("HIST_N_REPLIQUES " + json.dumps(dict(sorted(collections.Counter(P["n_repliques"] for P in inv).items()))))
for nom, filtre in (("210", lambda P: True), ("passent", lambda P: P["passe_les_criteres"]), ("selection", lambda P: P["accession"] in vsel)):
    sub = [P for P in inv if filtre(P)]
    print("CATEGORIE[%s] %s" % (nom, json.dumps(dict(sorted(collections.Counter(str(P["categorie_code"]) for P in sub).items())))))
    print("DIFFICULTE[%s] %s" % (nom, json.dumps(dict(sorted(collections.Counter(str(P["classe_difficulte"]) for P in sub).items())))))
    tail = [P["n_residus_total"] for P in sub if P.get("n_residus_total")]
    if tail:
        print("TAILLE[%s] n=%d min=%d mediane=%.1f max=%d somme=%d" % (nom, len(tail), min(tail), statistics.median(tail), max(tail), sum(tail)))
# ---- rupture CA-CA dans les structures MD : combien de chaines, et gaps RCSB correspondants
n_ch = n_rupt = n_gap = n_gap_sans_rupt = 0
dmax = 0.0
pres = intern = absents = 0
for P in inv:
    if not P["rcsb"] or not P["rcsb"].get("chaines"):
        continue
    for c in P["rcsb"]["chaines"]:
        if not c["evaluable"]:
            continue
        ch = next(x for x in P["chaines"] if x["id"] == c["chaine"])
        n_ch += 1
        dmax = max(dmax, ch["ca_dmax_A"] or 0)
        if ch["ca_ruptures_md"] > 0:
            n_rupt += 1
        if c["n_non_observes_internes"] > 0:
            n_gap += 1
            intern += c["n_non_observes_internes"]
            pres += c["internes_presents_dans_MD"]
            if ch["ca_ruptures_md"] == 0:
                n_gap_sans_rupt += 1
        absents += c["residus_absents_du_span_MD"]
print("CA_RUPTURES chaines_evaluables=%d  avec_rupture_CA_gt_4.5A_dans_structure.pdb=%d  dmax_CA_CA_A=%.3f" % (n_ch, n_rupt, dmax))
print("GAPS_RCSB chaines_avec_residus_non_observes_internes=%d  dont_sans_aucune_rupture_dans_la_MD=%d" % (n_gap, n_gap_sans_rupt))
print("GAPS_RCSB residus_internes_non_observes=%d  dont_presents_dans_la_MD(modelises)=%d  residus_du_span_absents_de_la_MD=%d" % (intern, pres, absents))
# ---- appariement
cm = collections.Counter()
for P in inv:
    for c in ((P["rcsb"] or {}).get("chaines") or []):
        cm[(c.get("methode_appariement") or "AUCUN", "meme_lettre" if c.get("meme_lettre") else "lettre_differente")] += 1
print("APPARIEMENT " + json.dumps({"%s/%s" % k: v for k, v in sorted(cm.items(), key=lambda kv: str(kv[0]))}))
# ---- les 12 de staging
for pdb in STAGING12:
    P = next(p for p in inv if p["pdb"] == pdb)
    print("STAGING12 %s %s passe=%s motifs=%s gaps_internes=%s paires5A=%s total=%s" % (
        pdb, P["accession"], P["passe_les_criteres"], ",".join(P["motifs_rejet"]) or "-", P.get("gaps_internes_total"),
        P.get("interface_paires_5A"), P.get("n_residus_total")))
# ---- tous les VHH / single-domain
for P in sorted(inv, key=lambda p: p["pdb"]):
    if P["vhh"] or P["vnar"] or P["categorie_code"] == "AS":
        ab = [x for x in P["vhh_chaines"] if x["vhh"] or x["vnar"]]
        print("ANTICORPS_DOMAINE_UNIQUE %s %s cat=%s vhh=%s vnar=%s chaines_prot=%s total=%s passe=%s motifs=%s selectionne=%s nom_ab=%s" % (
            P["pdb"], P["accession"], P["categorie_code"], P["vhh"], P["vnar"], P["structure"]["n_chaines_proteiques"], P.get("n_residus_total"),
            P["passe_les_criteres"], ",".join(P["motifs_rejet"]) or "-", P["accession"] in vsel,
            json.dumps([(x["description"], x["organismes"]) for x in ab], ensure_ascii=False)))
# ---- rejetes UNIQUEMENT par le gap interne : candidats a une tolerance (separe, ne change aucun critere)
solo = sorted([P for P in inv if P["motifs_rejet"] == ["gap_interne_rcsb"]], key=lambda p: (p["gaps_internes_total"], p["n_residus_total"]))
print("REJETES_SEULEMENT_PAR_GAP n=%d ; gap<=1: %d ; gap<=3: %d ; gap<=5: %d ; gap<=10: %d" % (
    len(solo), sum(1 for p in solo if p["gaps_internes_total"] <= 1), sum(1 for p in solo if p["gaps_internes_total"] <= 3),
    sum(1 for p in solo if p["gaps_internes_total"] <= 5), sum(1 for p in solo if p["gaps_internes_total"] <= 10)))
for P in solo:
    modelises = sum(c["internes_presents_dans_MD"] for c in P["rcsb"]["chaines"] if c["evaluable"])
    print("PRESQUE %s %s gap_residus=%d dont_modelises_dans_MD=%d total=%d cat=%s cls=%s vhh=%s staging12=%s" % (
        P["pdb"], P["accession"], P["gaps_internes_total"], modelises, P["n_residus_total"], P["categorie_code"], P["classe_difficulte"],
        P["vhh"], P["dans_staging12"]))
# ---- marges de l'interface : tout systeme a moins de 40 % au-dessus du seuil, parmi ceux qui passent
for P in sorted([p for p in inv if p["passe_les_criteres"] and p["interface_paires_5A"] < 35], key=lambda p: p["interface_paires_5A"]):
    print("MARGE_INTERFACE %s %s paires5A=%d (seuil 25) selectionne=%s" % (P["pdb"], P["accession"], P["interface_paires_5A"], P["accession"] in vsel))
for P in inv:
    if "interface_<25" in P["motifs_rejet"]:
        print("REJET_INTERFACE %s %s paires5A=%s total=%s" % (P["pdb"], P["accession"], P["interface_paires_5A"], P.get("n_residus_total")))
for P in inv:
    if "taille_>1030" in P["motifs_rejet"]:
        print("REJET_TAILLE %s %s total=%s chaines=%s" % (P["pdb"], P["accession"], P["n_residus_total"], P["chaines"] and [(c["id"], c["n_res"]) for c in P["chaines"]]))
for P in inv:
    if "rcsb_indisponible" in P["motifs_rejet"]:
        print("REJET_RCSB_INDISPONIBLE %s %s motifs=%s" % (P["pdb"], P["accession"], P["motifs_rejet"]))
# ---- resolution des selectionnes
res = [P["rcsb"]["resolution_A"] for P in inv if P["accession"] in vsel]
print("RESOLUTION[selection] n=%d min=%.2f mediane=%.2f max=%.2f" % (len(res), min(res), statistics.median(res), max(res)))
# ---- fiche de chaque systeme retenu (une ligne = une citation)
ordre = {a["accession"]: i for i, a in enumerate(sel["selection"])}
for P in sorted([p for p in inv if p["accession"] in vsel], key=lambda p: ordre[p["accession"]]):
    ch = ";".join("%s:%d" % (c["id"], c["n_res"]) for c in P["chaines"])
    rc = P["rcsb"]["chaines"]
    print("FICHE %-6s %s raison=%s chaines=[%s] total=%d cat=%s(%s) classe=%s irmsd=%s vhh=%s vnar=%s resolution=%s %s paires5A_rep1_structure=%d gaps_internes_rcsb=%d absents_du_span_MD=%d" % (
        P["pdb"], P["accession"], P["raison_selection"], ch, P["n_residus_total"], P["categorie_code"], P["categorie_label_hypothese"],
        P["classe_difficulte"], P["irmsd_db55"], P["vhh"], P["vnar"], P["rcsb"]["resolution_A"], P["rcsb"]["methode"],
        P["interface_paires_5A"], P["gaps_internes_total"], sum(c["residus_absents_du_span_MD"] for c in rc)))
    for c in rc:
        print("FICHE_CHAINE %-6s chaine_MD=%s -> %s (%s, identite=%.3f, meme_lettre=%s) %s | %s | non_observes_total=%d internes=%d absents_du_span_MD=%d" % (
            P["pdb"], c["chaine"], c["instance_rcsb"], c["methode_appariement"], c["identite"], c["meme_lettre"], c["description"],
            "/".join(c["organismes"]), c["n_non_observes_total"], c["n_non_observes_internes"], c["residus_absents_du_span_MD"]))
for P in [p for p in inv if p["accession"] in vsel]:
    for c in P["rcsb"]["chaines"]:
        if c["residus_absents_du_span_MD"] or c.get("fraction_mappee", 1) < 1:
            print("ALERTE_SELECTION %s chaine=%s residus_absents_du_span_MD=%d fraction_mappee=%s identite=%s" % (P["pdb"], c["chaine"], c["residus_absents_du_span_MD"], c.get("fraction_mappee"), c.get("identite")))
# ======== apres telechargement : ce que disent les 90 repliques ========
vr = [json.loads(l) for l in open(CACHE + "/verif_replicas.jsonl", encoding="utf-8")]
print("REPLIQUES verifiees=%d  conformes=%d  anomalies=%d" % (len(vr), sum(1 for r in vr if r["conforme"]), sum(1 for r in vr if not r["conforme"])))
print("REPLIQUES images=%s atomes_eq_PDB=%s chaines_prot=%s non_prot_total=%d eau_total=%d ruptures_CA_total=%d" % (
    sorted({r["n_images"] for r in vr}), all(r["n_atomes"] == r["n_atomes_pdb"] for r in vr), sorted({r["n_chaines_proteiques"] for r in vr}),
    sum(r["n_residus_non_proteine"] for r in vr), sum(r["n_eau"] for r in vr),
    sum(x["images_avec_rupture_gt_4.5A"] for r in vr for x in r["ca_ruptures_par_chaine"])))
print("REPLIQUES time_axe premier=%s dernier=%s pas=%s" % (sorted({r["time_premier"] for r in vr}), sorted({r["time_dernier"] for r in vr}), sorted({r["time_pas"] for r in vr})))
i0 = [r["paires_5A"]["image0"] for r in vr]
dip = [(r["pdb"], r["replique"], r["paires_5A"]["min"]) for r in vr if r["paires_5A"]["min"] < 25]
print("INTERFACE_TRAJ image0: min=%d mediane=%.1f ; repliques_dont_une_image_echantillonnee_<25: %d/%d ; systemes_concernes=%d" % (
    min(i0), statistics.median(i0), len(dip), len(vr), len({d[0] for d in dip})))
for d in dip:
    print("INTERFACE_TRAJ_SOUS_25 %s rep%d min_sur_20_images=%d" % d)
vs = jload(CACHE + "/verif_systemes.json")["systemes"]
print("ENTRE_REPLIQUES systemes=%d  toutes_distinctes=%s  ecart_image0_min_nm=%.3f  ecart_image500_min_nm=%.3f" % (
    len(vs), all(v["repliques_distinctes"] for v in vs.values()),
    min(x["image0_nm"] for v in vs.values() for x in v["ecart_entre_repliques_nm"].values()),
    min(x["image500_nm"] for v in vs.values() for x in v["ecart_entre_repliques_nm"].values())))
for pdb, v in vs.items():
    if v.get("staging"):
        print("DCD_STAGING %s dcd_est_replique_1=%s dcd_vs_xtc_nm=%s dcd_octets=%d" % (pdb, v["staging"]["dcd_est_replique_1"], json.dumps(v["staging"]["dcd_vs_xtc_nm_images_0_500_999"]), v["staging"]["dcd_octets"]))
s0 = [json.loads(l) for l in open(CACHE + "/verif_structure_vs_image0.jsonl", encoding="utf-8")]
ident = [r for r in s0 if r["ecart_brut_max_nm"] < 0.002]
diff = [r for r in s0 if r["ecart_brut_max_nm"] >= 0.002]
print("STRUCTUREPDB_VS_IMAGE0 repliques=%d  identiques_a_l_image0=%d  differentes=%d  (rmsd_apres_superposition: min=%.3f mediane=%.3f max=%.3f A ; translation_max=%.2f A)" % (
    len(s0), len(ident), len(diff), min(r["rmsd_apres_superposition_image0_A"] for r in diff), statistics.median(r["rmsd_apres_superposition_image0_A"] for r in diff),
    max(r["rmsd_apres_superposition_image0_A"] for r in diff), max(r["translation_centroide_A"] for r in diff)))
print("STRUCTUREPDB_VS_IMAGE0 parmi_les_differentes: structure.pdb egale a une image 1..30 du xtc (rmsd<0.05A)=%d ; sans image proche (>=0.05A)=%d" % (
    sum(1 for r in diff if r["rmsd_image_la_plus_proche_A"] < 0.05), sum(1 for r in diff if r["rmsd_image_la_plus_proche_A"] >= 0.05)))
dl = [json.loads(l) for l in open(CACHE + "/download_log.jsonl", encoding="utf-8")]
last = {}
for r in dl:
    last[r["chemin"]] = r
print("TELECHARGEMENT fichiers=%d  verifies=%d  octets=%d (%.3f Go)  double_telechargement_identique=%d  passes_3=%d" % (
    len(last), sum(1 for r in last.values() if r["verifie"]), sum(r["octets"] for r in last.values()), sum(r["octets"] for r in last.values()) / 1e9,
    sum(1 for r in last.values() if r["egal"]), sum(1 for r in last.values() if r.get("sha256_3"))))
ext = collections.Counter(os.path.splitext(r["chemin"])[1] for r in last.values())
print("TELECHARGEMENT par_extension=%s" % json.dumps(dict(ext)))
xt = [r["debit_Mo_s_1"] for r in last.values() if r["fichier"] == "trajectory.xtc" and r.get("debit_Mo_s_1")]
print("TELECHARGEMENT debit_xtc_Mo_s min=%.1f mediane=%.1f max=%.1f" % (min(xt), statistics.median(xt), max(xt)))
psf = collections.defaultdict(set)
for r in last.values():
    if r["fichier"] == "topology.psf":
        psf[r["pdb"]].add(r["sha256"])
print("PSF systemes=%d  psf_identique_entre_repliques=%d" % (len(psf), sum(1 for v in psf.values() if len(v) == 1)))
# ---- passants non retenus (pour reequilibrer la difficulte si le guide le demande)
nr = sorted([P for P in inv if P["passe_les_criteres"] and P["accession"] not in vsel], key=lambda p: (str(p["classe_difficulte"]), p["n_residus_total"]))
print("PASSANTS_NON_RETENUS n=%d  par_classe=%s" % (len(nr), json.dumps(dict(collections.Counter(str(p["classe_difficulte"]) for p in nr)))))
for P in nr:
    if P["classe_difficulte"] in ("Difficult", "Medium Difficulty"):
        print("PASSANT_NON_RETENU %s %s classe=%s cat=%s total=%d paires5A=%d" % (P["pdb"], P["accession"], P["classe_difficulte"], P["categorie_code"], P["n_residus_total"], P["interface_paires_5A"]))
print("FINI s5")
