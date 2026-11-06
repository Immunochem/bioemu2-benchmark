"""J4 — manifeste JSONL append-only : une ligne par systeme retenu.

Rien n'est jamais reecrit : une ligne dont le contenu STABLE (hors dates) existe deja est sautee ; une ligne dont le
contenu a change est AJOUTEE avec `remplace_digest`. Les sha256 sont RECALCULES sur le disque au moment de l'ecriture
(la promesse de reutilisabilite se verifie, elle ne se recopie pas) et compares a ceux du telechargement.
Chaque ligne porte la promesse complete : machine . chemin absolu . compte + extension . date de verification.
"""
import collections
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import CACHE, ROOT, jload, now_iso, sha256_file  # noqa: E402

MANIFESTE = ROOT + "/manifeste_corpus.json"
MACHINE = "staging-server"


def digest(rec):
    stable = json.loads(json.dumps({k: v for k, v in rec.items()
                                    if k not in ("ecrit_le", "verifie_le", "digest", "remplace_digest", "promesse")}))
    for f in stable.get("fichiers", []):
        f.pop("verifie_sur_disque_le", None)
    return hashlib.sha256(json.dumps(stable, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def main():
    sel = jload(ROOT + "/selection.json")["selection"]
    inv = {P["accession"]: P for P in jload(ROOT + "/inventaire_dynabench.json")["projets"]}
    verif_sys = jload(CACHE + "/verif_systemes.json")["systemes"]
    verif_rep = collections.defaultdict(dict)
    for line in open(CACHE + "/verif_replicas.jsonl", encoding="utf-8"):
        r = json.loads(line)
        verif_rep[r["pdb"]][r["replique"]] = r
    s0 = collections.defaultdict(dict)
    if os.path.exists(CACHE + "/verif_structure_vs_image0.jsonl"):
        for line in open(CACHE + "/verif_structure_vs_image0.jsonl", encoding="utf-8"):
            r = json.loads(line)
            s0[r["pdb"]][r["replique"]] = {k: r[k] for k in ("ecart_brut_max_nm", "translation_centroide_A", "rmsd_apres_superposition_image0_A",
                                                             "image_la_plus_proche", "rmsd_image_la_plus_proche_A")}
    dl = collections.defaultdict(dict)
    for line in open(CACHE + "/download_log.jsonl", encoding="utf-8"):
        r = json.loads(line)
        dl[r["pdb"]][r["chemin"]] = r            # la derniere ligne par chemin fait foi
    existants = {}
    if os.path.exists(MANIFESTE):
        for line in open(MANIFESTE, encoding="utf-8"):
            r = json.loads(line)
            existants.setdefault(r["pdb"] + "/" + r["accession"], []).append(r)
    deja = {r["digest"] for rs in existants.values() for r in rs}
    n_ecrit = n_saute = 0
    date = now_iso()
    with open(MANIFESTE, "a", encoding="utf-8") as fo:
        for s in sel:
            pdb, acc = s["pdb"], s["accession"]
            P = inv[acc]
            base = "%s/systemes/%s" % (ROOT, pdb)
            fichiers = []
            ext = collections.Counter()
            tous_egaux = True
            for dp, _dn, fn in sorted(os.walk(base)):
                for f in sorted(fn):
                    p = os.path.join(dp, f)
                    h = sha256_file(p)
                    rec_dl = dl[pdb].get(p)
                    egal_dl = (rec_dl is not None and rec_dl.get("sha256") == h) if rec_dl else None
                    if rec_dl is not None and not egal_dl:
                        tous_egaux = False
                    rep = os.path.basename(dp)
                    fichiers.append({"chemin": p, "replique": rep, "extension": os.path.splitext(f)[1], "octets": os.path.getsize(p),
                                     "sha256": h, "declare_serveur_octets": rec_dl.get("declare_serveur") if rec_dl else None,
                                     "double_telechargement_identique": rec_dl.get("egal") if rec_dl else None,
                                     "sha256_egal_au_telechargement": egal_dl,
                                     "provenance": ("telecharge " + rec_dl["url"]) if rec_dl else "copie du cache API L0a (filenotes)",
                                     "verifie_sur_disque_le": date})
                    ext[os.path.splitext(f)[1]] += 1
            vr = verif_rep.get(pdb, {})
            vs = verif_sys.get(pdb, {})
            anomalies = {("rep%d" % k): v["anomalies"] for k, v in vr.items() if v.get("anomalies")}
            verifie = bool(fichiers) and tous_egaux and all(
                (f["double_telechargement_identique"] is True) for f in fichiers if f["double_telechargement_identique"] is not None) \
                and all(v.get("conforme") for v in vr.values()) and bool(vr)
            chaines = []
            for c in P["chaines"]:
                rc = next((x for x in (P["rcsb"] or {}).get("chaines", []) if x["chaine"] == c["id"]), {})
                chaines.append({"id": c["id"], "n_residus": c["n_res"], "resseq_premier": c["resseq_premier"], "resseq_dernier": c["resseq_dernier"],
                                "description_rcsb": rc.get("description"), "organismes": rc.get("organismes"),
                                "instance_rcsb": rc.get("instance_rcsb"), "appariement": rc.get("methode_appariement"),
                                "residus_non_observes_internes_rcsb": rc.get("n_non_observes_internes")})
                if rc.get("residus_absents_du_span_MD"):
                    chaines[-1]["residus_absents_du_span_MD"] = rc["residus_absents_du_span_MD"]
                    chaines[-1]["note_absents"] = ("positions d'entite RCSB comprises entre le premier et le dernier residu de la chaine MD "
                                                   "et ABSENTES de structure.pdb (observees dans le cristal : non_observes_total=%s) ; "
                                                   "cause non verifiee" % rc.get("n_non_observes_total"))
            rec = {
                "schema": "corpus-bioemu2-dynabench/1", "machine": MACHINE, "pdb": pdb, "accession": acc,
                "raison_selection": s["raison"],
                "chaines": chaines, "n_residus_total": P["n_residus_total"],
                "categorie_code": P["categorie_code"], "categorie_label_hypothese": P["categorie_label_hypothese"],
                "classe_difficulte": P["classe_difficulte"], "irmsd_db55_A": P["irmsd_db55"], "complexe_db55": P["complexe_db55"],
                "partenaires_db55": [P["partenaire1"], P["partenaire2"]],
                "vhh": P["vhh"], "vnar": P["vnar"],
                "vhh_evidence": [{"chaine": x["chaine"], "description": x["description"], "organismes": x["organismes"],
                                  "domaine_V_fr4": x["V_like"], "camelide": x["camelide"], "mot_cle": x["mot_cle"], "cat_AS": x["cat_AS"]}
                                 for x in P["vhh_chaines"] if x["vhh"]],
                "interface_paires_5A": {"structure_pdb_image0": P["interface_paires_5A"],
                                        "trajectoire_par_replique": {("rep%d" % k): v.get("paires_5A") for k, v in vr.items()}},
                "rcsb": {"methode": P["rcsb"]["methode"], "resolution_A": P["rcsb"]["resolution_A"], "titre": P["rcsb"]["titre"],
                         "residus_non_observes_internes_total": P["gaps_internes_total"]},
                "repliques": sorted(vr.keys()),
                "dossier": base, "fichiers": fichiers, "compte_fichiers": {"total": len(fichiers), **{k.lstrip("."): v for k, v in sorted(ext.items())}},
                "fichiers_staging_existants": vs.get("staging"),
                "verification_contenu": {("rep%d" % k): {kk: v.get(kk) for kk in (
                    "n_images", "n_atomes", "n_chaines_proteiques", "n_residus_proteine", "n_residus_non_proteine", "n_eau",
                    "time_premier", "time_dernier", "time_pas", "ca_ruptures_par_chaine", "distance_centres_chaines_A",
                    "ecart_max_xtc0_vs_structurepdb_nm", "conforme")} for k, v in vr.items()},
                "structure_pdb_vs_image0_du_xtc": {("rep%d" % k): v for k, v in sorted(s0.get(pdb, {}).items())},
                "ecart_entre_repliques_nm": vs.get("ecart_entre_repliques_nm"), "repliques_distinctes": vs.get("repliques_distinctes"),
                "anomalies": anomalies, "verifie": verifie,
                "axe_time_xtc": "numero d'image 0..999, PAS des ns ; 0,1 ns par image (FRAMESTEP API) ; 100 ns par replique",
                "licence": "indeterminee (LICENSE: null dans l'API) : usage non commercial, ne pas republier les trajectoires",
                "verifie_le": date,
            }
            rec["promesse"] = ("machine=%s . dossier=%s . %s . verifie le %s" % (
                MACHINE, base, " + ".join("%d .%s" % (v, k) for k, v in rec["compte_fichiers"].items() if k != "total")
                + " (%d fichiers)" % rec["compte_fichiers"]["total"], date))
            rec["digest"] = digest(rec)
            prev = [r["digest"] for r in existants.get(pdb + "/" + acc, [])]
            if rec["digest"] in deja:
                n_saute += 1
                print("MANIFESTE saute %s : contenu identique deja ecrit" % pdb, flush=True)
                continue
            if prev:
                rec["remplace_digest"] = prev[-1]
            rec["ecrit_le"] = date
            fo.write(json.dumps(rec, ensure_ascii=False) + "\n")
            fo.flush()
            n_ecrit += 1
            print("MANIFESTE ecrit %s total=%d verifie=%s anomalies=%s compte=%s" % (
                pdb, len(fichiers), verifie, anomalies or "-", json.dumps(rec["compte_fichiers"])), flush=True)
    print("BILAN manifeste lignes_ecrites=%d sautees=%d  fichier=%s" % (n_ecrit, n_saute, MANIFESTE), flush=True)
    print("FINI s4 %s" % now_iso(), flush=True)


if __name__ == "__main__":
    main()
