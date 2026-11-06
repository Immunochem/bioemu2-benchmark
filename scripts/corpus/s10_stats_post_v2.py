"""Statistiques v2 APRES telechargement du delta : on relit les fichiers `_v2` ecrits par s9 (aucun recalcul)."""
import collections
import json
import os
import re
import statistics
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import CACHE, L0A, ROOT, jload  # noqa: E402

sel = jload(ROOT + "/selection_v2.json")["selection"]
statut = {s["pdb"]: s["statut_vs_v1"] for s in sel}
vr = [json.loads(l) for l in open(CACHE + "/verif_replicas_v2.jsonl", encoding="utf-8")]
print("REPLIQUES_V2 verifiees=%d  conformes=%d  anomalies=%d  systemes=%d" % (len(vr), sum(1 for r in vr if r["conforme"]), sum(1 for r in vr if not r["conforme"]), len({r["pdb"] for r in vr})))
print("REPLIQUES_V2 images=%s atomes_eq_PDB=%s chaines_prot=%s non_prot_total=%d eau_total=%d ruptures_CA_total=%d" % (
    sorted({r["n_images"] for r in vr}), all(r["n_atomes"] == r["n_atomes_pdb"] for r in vr), sorted({r["n_chaines_proteiques"] for r in vr}),
    sum(r["n_residus_non_proteine"] for r in vr), sum(r["n_eau"] for r in vr), sum(x["images_avec_rupture_gt_4.5A"] for r in vr for x in r["ca_ruptures_par_chaine"])))
print("REPLIQUES_V2 time_axe premier=%s dernier=%s pas=%s  dmax_CA_CA_A=%.3f" % (sorted({r["time_premier"] for r in vr}), sorted({r["time_dernier"] for r in vr}), sorted({r["time_pas"] for r in vr}),
                                                                           max(x["dmax_A"] for r in vr for x in r["ca_ruptures_par_chaine"])))
i0 = [r["paires_5A"]["image0"] for r in vr]
dip = [(r["pdb"], r["replique"], r["paires_5A"]["min"]) for r in vr if r["paires_5A"]["min"] < 25]
print("INTERFACE_V2 image0: min=%d mediane=%.1f ; repliques_avec_une_image_echantillonnee_<25: %d/%d ; systemes_concernes=%d" % (min(i0), statistics.median(i0), len(dip), len(vr), len({d[0] for d in dip})))
for d in sorted(dip):
    print("INTERFACE_V2_SOUS_25 %s rep%d min_sur_20_images=%d" % d)
par = collections.defaultdict(dict)
for r in vr:
    par[r["pdb"]][r["replique"]] = r["paires_5A"]
for s in sel:
    p = par[s["pdb"]]
    print("INTERFACE_PLAGE %-5s %-6s min=%3d  max=%3d  | %s" % (s["pdb"], statut[s["pdb"]], min(x["min"] for x in p.values()), max(x["max"] for x in p.values()),
                                                             "  ".join("rep%d[%d/%d/%d]" % (k, v["min"], v["mediane"], v["max"]) for k, v in sorted(p.items()))))
vs = jload(CACHE + "/verif_systemes_v2.json")["systemes"]
print("ENTRE_REPLIQUES_V2 systemes=%d  toutes_distinctes=%s  ecart_image0_min_nm=%.3f  ecart_image500_min_nm=%.3f" % (
    len(vs), all(v["repliques_distinctes"] for v in vs.values()), min(x["image0_nm"] for v in vs.values() for x in v["ecart_entre_repliques_nm"].values()),
    min(x["image500_nm"] for v in vs.values() for x in v["ecart_entre_repliques_nm"].values())))
n_dcd = 0
for pdb, v in vs.items():
    if v.get("staging"):
        n_dcd += 1
        print("DCD_STAGING_V2 %s dcd_est_replique_1=%s dcd_vs_xtc_nm=%s dcd_octets=%d" % (pdb, v["staging"]["dcd_est_replique_1"], json.dumps(v["staging"]["dcd_vs_xtc_nm_images_0_500_999"]), v["staging"]["dcd_octets"]))
print("DCD_STAGING_V2 systemes_testes=%d  dcd_est_replique_1_pour_tous=%s" % (n_dcd, all(v["staging"]["dcd_est_replique_1"] for v in vs.values() if v.get("staging"))))
s0 = [json.loads(l) for l in open(CACHE + "/verif_structure_vs_image0_v2.jsonl", encoding="utf-8")]
ident = [r for r in s0 if r["ecart_brut_max_nm"] < 0.002]
diff = [r for r in s0 if r["ecart_brut_max_nm"] >= 0.002]
print("STRUCTUREPDB_VS_IMAGE0_V2 repliques=%d  identiques_a_l_image0=%d  differentes=%d  (rmsd_apres_superposition: min=%.3f mediane=%.3f max=%.3f A ; translation_max=%.2f A)" % (
    len(s0), len(ident), len(diff), min(r["rmsd_apres_superposition_image0_A"] for r in diff), statistics.median(r["rmsd_apres_superposition_image0_A"] for r in diff),
    max(r["rmsd_apres_superposition_image0_A"] for r in diff), max(r["translation_centroide_A"] for r in diff)))
print("STRUCTUREPDB_VS_IMAGE0_V2 parmi_les_differentes: egale_a_une_image_proche(rmsd<0.05A)=%d ; sans_image_proche=%d" % (
    sum(1 for r in diff if r["rmsd_image_la_plus_proche_A"] < 0.05), sum(1 for r in diff if r["rmsd_image_la_plus_proche_A"] >= 0.05)))
print("STRUCTUREPDB_VS_IMAGE0_V2 par_statut: " + json.dumps({st: [sum(1 for r in diff if statut[r["pdb"]] == st), sum(1 for r in s0 if statut[r["pdb"]] == st)] for st in ("entre", "reste")}))
dl = [json.loads(l) for l in open(CACHE + "/download_log.jsonl", encoding="utf-8")]
last = {}
for r in dl:
    last[r["chemin"]] = r
v2 = {s["pdb"] for s in sel}
mine = {c: r for c, r in last.items() if r["pdb"] in v2}
entr = {c: r for c, r in mine.items() if statut[r["pdb"]] == "entre"}
rest = {c: r for c, r in mine.items() if statut[r["pdb"]] == "reste"}
print("TELECHARGEMENT_V2 fichiers=%d verifies=%d  double_telechargement_identique=%d  passes_3=%d  octets=%d (%.3f Go)" % (
    len(mine), sum(1 for r in mine.values() if r["verifie"]), sum(1 for r in mine.values() if r["egal"]), sum(1 for r in mine.values() if r.get("sha256_3")),
    sum(r["octets"] for r in mine.values()), sum(r["octets"] for r in mine.values()) / 1e9))
print("TELECHARGEMENT_V2 entrants: fichiers=%d telecharges_ce_passage(reprise=False)=%d octets=%d (%.3f Go) ; restent: fichiers=%d repris_sans_retelechargement(reprise=True)=%d" % (
    len(entr), sum(1 for r in entr.values() if not r.get("reprise")), sum(r["octets"] for r in entr.values()), sum(r["octets"] for r in entr.values()) / 1e9,
    len(rest), sum(1 for r in rest.values() if r.get("reprise"))))
print("TELECHARGEMENT_V2 par_extension=%s" % json.dumps(dict(collections.Counter(os.path.splitext(c)[1] for c in mine))))
psf = collections.defaultdict(set)
for r in mine.values():
    if r["fichier"] == "topology.psf":
        psf[r["pdb"]].add(r["sha256"])
print("PSF_V2 systemes=%d  psf_identique_entre_repliques=%d" % (len(psf), sum(1 for v in psf.values() if len(v) == 1)))
# recoupement avec L0a (telechargements independants, plus tot dans la journee)
sums = {}
for l in open("/gpfs/phherman/bioemu2_bench/dynabench/SHA256SUMS", encoding="utf-8"):
    h, _, p = l.strip().partition("  ")
    m = re.search(r"([0-9A-Z]{4}/rep\d/[a-z_.]+)$", p)
    if m:
        sums[m.group(1)] = h
bon = mauvais = 0
for c, r in last.items():
    if r["pdb"] in ("6CWG", "1B6C", "2OT3") and r["fichier"] in ("trajectory.xtc", "structure.pdb", "topology.psf"):
        k = "%s/rep%d/%s" % (r["pdb"], r["replique"], r["fichier"])
        if k in sums:
            if sums[k] == r["sha256"]:
                bon += 1
            else:
                mauvais += 1
                print("RECOUPEMENT_L0A DIFFERENT %s" % k)
print("RECOUPEMENT_L0A fichiers_comparables(1B6C,6CWG,2OT3 x 3 repliques x 3 fichiers)=%d  identiques=%d  differents=%d" % (bon + mauvais, bon, mauvais))
print("FINI s10")
