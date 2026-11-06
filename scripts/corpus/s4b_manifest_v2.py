"""J7 — manifeste v2 : lignes AJOUTEES (jamais reecrites) au meme manifeste JSONL.

  - chaque systeme retenu en v2 recoit une nouvelle ligne (statut `actif_v2`) portant les champs ajoutes :
    commensurabilite (preparateur + sha256), residus non observes RCSB (internes reconstruits / absents, terminaux,
    DANS LES TROIS NUMEROTATIONS), plage d'occupation de l'interface sur les 20 images echantillonnees ;
  - chaque systeme v1 qui sort de la selection recoit une ligne `retire_v2` (ses fichiers restent sur le disque) ;
  - une ligne dont le contenu STABLE existe deja est sautee.
Les sha256 sont recalcules sur le disque. Pour lire le manifeste : derniere ligne de chaque `pdb`.
"""
import collections
import json
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from corpus_lib import CACHE, ROOT, jload, now_iso, sha256_file  # noqa: E402
import s4_manifest as s4  # noqa: E402  (digest, MANIFESTE, MACHINE : memes definitions qu'en v1)

NUMEROTATIONS = {
    "resSeq_structure_pdb": "numero de residu de `structure.pdb` (numerotation d'AUTEUR du cristal, SANS code d'insertion : mdtraj l'ecrase) ; "
                            "c'est le `resSeq` que lit le preparateur",
    "seq_id_rcsb": "`label_seq_id` RCSB = position 1..L dans l'ENTITE (sequence deposee) ; c'est ce que l'API RCSB "
                   "`UNOBSERVED_RESIDUE_XYZ` renvoie dans `feature_positions`",
    "auth_seq_id_rcsb": "numero d'auteur RCSB (`auth_to_entity_poly_seq_mapping`) ; dans la plupart des chaines = resSeq_structure_pdb ; "
                        "peut differer de seq_id_rcsb d'un decalage constant (ex. 6CWG chaine A : seq_id = auteur + 1 car l'auteur numerote depuis 0)",
    "ordinal_MD": "position 1..N du residu dans la chaine de `structure.pdb` ; c'est la cle de correspondance avec un echantillon BioEmu2 "
                  "(le preparateur apparie sur l'ordinal, jamais sur resSeq)",
    "regle": "un index n'est pas transferable d'une numerotation a l'autre : chaque run de residus non observes est donne dans les quatre",
}


def plage_interface(vr):
    reps = {}
    mins, maxs, med = [], [], []
    nimg = 0
    for k, v in sorted(vr.items()):
        p = v.get("paires_5A")
        if not p:
            continue
        reps["rep%d" % k] = {"image0": p["image0"], "min": p["min"], "mediane": p["mediane"], "max": p["max"],
                             "images_testees": p["images_testees"], "stride": p["stride"]}
        mins.append(p["min"])
        maxs.append(p["max"])
        med.append(p["mediane"])
        nimg += p["images_testees"]
    return {"par_replique": reps,
            "plage_globale": {"min": min(mins), "max": max(maxs), "mediane_des_medianes": sorted(med)[len(med) // 2], "images_echantillonnees": nimg,
                              "seuil_critere": 25,
                              "repliques_dont_le_min_est_sous_le_seuil": [r for r, x in reps.items() if x["min"] < 25]},
            "avertissement": "20 images par replique (une toutes les 50) : le minimum reel peut etre plus bas ; le critere statique (>= 25 sur structure.pdb rep1) est pose sur un objet qui bouge"}


def main():
    sel = jload(ROOT + "/selection_v2.json")
    selected = {s["accession"]: s for s in sel["selection"]}
    inv = {P["accession"]: P for P in jload(ROOT + "/inventaire_dynabench_v2.json")["projets"]}
    meta_inv = jload(ROOT + "/inventaire_dynabench_v2.json")["meta"]
    verif_sys = jload(CACHE + "/verif_systemes_v2.json")["systemes"]
    verif_rep = collections.defaultdict(dict)
    for line in open(CACHE + "/verif_replicas_v2.jsonl", encoding="utf-8"):
        r = json.loads(line)
        verif_rep[r["pdb"]][r["replique"]] = r
    s0 = collections.defaultdict(dict)
    for line in open(CACHE + "/verif_structure_vs_image0_v2.jsonl", encoding="utf-8"):
        r = json.loads(line)
        s0[r["pdb"]][r["replique"]] = {k: r[k] for k in ("ecart_brut_max_nm", "translation_centroide_A", "rmsd_apres_superposition_image0_A",
                                                         "image_la_plus_proche", "rmsd_image_la_plus_proche_A")}
    dl = collections.defaultdict(dict)
    for line in open(CACHE + "/download_log.jsonl", encoding="utf-8"):
        r = json.loads(line)
        dl[r["pdb"]][r["chemin"]] = r
    existants = collections.defaultdict(list)
    for line in open(s4.MANIFESTE, encoding="utf-8"):
        r = json.loads(line)
        existants[r["pdb"] + "/" + r["accession"]].append(r)
    deja = {r["digest"] for rs in existants.values() for r in rs}
    date = now_iso()
    n_ecrit = n_saute = 0
    with open(s4.MANIFESTE, "a", encoding="utf-8") as fo:
        def ecrire(rec, prev_key):
            nonlocal n_ecrit, n_saute
            rec["digest"] = s4.digest(rec)
            if rec["digest"] in deja:
                n_saute += 1
                print("MANIFESTE_V2 saute %s : contenu identique deja ecrit" % rec["pdb"], flush=True)
                return
            prev = [r["digest"] for r in existants.get(prev_key, [])]
            if prev:
                rec["remplace_digest"] = prev[-1]
            rec["ecrit_le"] = date
            fo.write(json.dumps(rec, ensure_ascii=False) + "\n")
            fo.flush()
            deja.add(rec["digest"])
            n_ecrit += 1
            print("MANIFESTE_V2 ecrit %s statut=%s total=%d verifie=%s" % (rec["pdb"], rec["statut"], rec["compte_fichiers"]["total"], rec["verifie"]), flush=True)

        # ---------------- systemes actifs v2 ----------------
        for s in sel["selection"]:
            pdb, acc = s["pdb"], s["accession"]
            P = inv[acc]
            base = "%s/systemes/%s" % (ROOT, pdb)
            fichiers, ext, tous_egaux = [], collections.Counter(), True
            for dp, _dn, fn in sorted(os.walk(base)):
                for f in sorted(fn):
                    p = os.path.join(dp, f)
                    h = sha256_file(p)
                    rec_dl = dl[pdb].get(p)
                    egal_dl = (rec_dl.get("sha256") == h) if rec_dl else None
                    if rec_dl is not None and not egal_dl:
                        tous_egaux = False
                    fichiers.append({"chemin": p, "replique": os.path.basename(dp), "extension": os.path.splitext(f)[1], "octets": os.path.getsize(p),
                                     "sha256": h, "declare_serveur_octets": rec_dl.get("declare_serveur") if rec_dl else None,
                                     "double_telechargement_identique": rec_dl.get("egal") if rec_dl else None,
                                     "sha256_egal_au_telechargement": egal_dl,
                                     "provenance": ("telecharge " + rec_dl["url"]) if rec_dl else "copie du cache API L0a (filenotes)",
                                     "verifie_sur_disque_le": date})
                    ext[os.path.splitext(f)[1]] += 1
            vr, vs = verif_rep.get(pdb, {}), verif_sys.get(pdb, {})
            anomalies = {("rep%d" % k): v["anomalies"] for k, v in vr.items() if v.get("anomalies")}
            verifie = bool(fichiers) and tous_egaux and bool(vr) and all(v.get("conforme") for v in vr.values()) and all(
                f["double_telechargement_identique"] is True for f in fichiers if f["double_telechargement_identique"] is not None)
            chaines = []
            for c in P["chaines"]:
                rc = next((x for x in P["rcsb_v2"]["chaines"] if x["chaine"] == c["id"]), {})
                cm = P["commensurabilite"].get(c["id"], {}).get("fenetre_2.9", {})
                chaines.append({"id": c["id"], "n_residus": c["n_res"], "resseq_premier": c["resseq_premier"], "resseq_dernier": c["resseq_dernier"],
                                "instance_rcsb": rc.get("instance_rcsb"), "appariement": rc.get("methode_appariement"),
                                "description_rcsb": next((x["description"] for x in P["rcsb"]["chaines"] if x["chaine"] == c["id"]), None),
                                "organismes": next((x["organismes"] for x in P["rcsb"]["chaines"] if x["chaine"] == c["id"]), None),
                                "commensurabilite": {
                                    "n_residus_preparateur": cm.get("n_obs"), "n_residus_to_fasta": P["commensurabilite"][c["id"]].get("n_fasta"),
                                    "sequence_identique_a_to_fasta": cm.get("accord_fasta_sequence"), "n_gaps_C_N": cm.get("n_gaps"),
                                    "liaison_peptidique_max_A": cm.get("liaison_max_A"), "n_residus_modifies_lies": cm.get("n_modifies"),
                                    "n_residus_X": cm.get("n_X"), "n_sauts_resSeq_sans_rupture": cm.get("n_sauts_resSeq_sans_rupture")},
                                "residus_non_observes_rcsb": ({
                                    "internes_encadres_par_des_residus_observes": rc["internes"],
                                    "n_terminaux": rc["n_terminaux"], "c_terminaux": rc["c_terminaux"]} if rc.get("evaluable") else None)})
            rec = {
                "schema": "corpus-bioemu2-dynabench/2", "criteres_version": "v2-commensurabilite", "statut": "actif_v2", "statut_vs_v1": s["statut_vs_v1"],
                "machine": s4.MACHINE, "pdb": pdb, "accession": acc, "raison_selection": s["raison"],
                "chaines": chaines, "n_residus_total": P["n_residus_total"],
                "categorie_code": P["categorie_code"], "categorie_label_hypothese": P["categorie_label_hypothese"],
                "classe_difficulte": P["classe_difficulte"], "irmsd_db55_A": P["irmsd_db55"], "complexe_db55": P["complexe_db55"],
                "partenaires_db55": [P["partenaire1"], P["partenaire2"]], "vhh": P["vhh"], "vnar": P["vnar"],
                "vhh_evidence": [{"chaine": x["chaine"], "description": x["description"], "organismes": x["organismes"], "domaine_V_fr4": x["V_like"],
                                  "camelide": x["camelide"], "mot_cle": x["mot_cle"], "cat_AS": x["cat_AS"]} for x in P["vhh_chaines"] if x["vhh"]],
                "commensurabilite_critere": {
                    "definition": "(a) aucune rupture de chaine C(i)-N(i+1) dans structure.pdb (preparateur, fenetre par defaut) ET (b) meme nombre de residus "
                                  "et meme sequence que topology.to_fasta() de la meme structure",
                    "preparateur_sha256": meta_inv["preparateur_sha256"], "preparateur_origine_sha256": meta_inv["preparateur_origine_sha256"],
                    "version_annoncee_non_retrouvee": "f8e45af307ccb55d... (le fichier lu sur le disque a un autre sha256)",
                    "fenetre_C_N_angstrom": [1.0, 2.9], "resultat_fenetre_1.9": P["commensurable_fenetre_1.9"],
                    "resultat_preparateur_origine_production": P["commensurable_preparateur_origine"],
                    "resume": P["commensurabilite_resume"]},
                "residus_non_observes_rcsb_declaration": "les residus INTERNES non observes qui sont PRESENTS dans la MD ont ete RECONSTRUITS avant la MD ; "
                                                         "l'auteur de cette reconstruction (DB5.5, DynaBench, preparation NAMD) n'est PAS verifie ; "
                                                         "les residus internes ABSENTS de la MD l'ont ete retires sans rupture de chaine",
                "numerotations": NUMEROTATIONS,
                "interface_paires_5A": {"structure_pdb_rep1": P["interface_paires_5A"], "occupation_trajectoire": plage_interface(vr)},
                "rcsb": {"methode": P["rcsb"]["methode"], "resolution_A": P["rcsb"]["resolution_A"], "titre": P["rcsb"]["titre"]},
                "repliques": sorted(vr.keys()), "dossier": base, "fichiers": fichiers,
                "compte_fichiers": {"total": len(fichiers), **{k.lstrip("."): v for k, v in sorted(ext.items())}},
                "fichiers_staging_existants": vs.get("staging"),
                "verification_contenu": {("rep%d" % k): {kk: v.get(kk) for kk in (
                    "n_images", "n_atomes", "n_chaines_proteiques", "n_residus_proteine", "n_residus_non_proteine", "n_eau", "time_premier", "time_dernier",
                    "time_pas", "ca_ruptures_par_chaine", "distance_centres_chaines_A", "ecart_max_xtc0_vs_structurepdb_nm", "conforme")} for k, v in vr.items()},
                "structure_pdb_vs_image0_du_xtc": {("rep%d" % k): v for k, v in sorted(s0.get(pdb, {}).items())},
                "ecart_entre_repliques_nm": vs.get("ecart_entre_repliques_nm"), "repliques_distinctes": vs.get("repliques_distinctes"),
                "anomalies": anomalies, "verifie": verifie,
                "axe_time_xtc": "numero d'image 0..999, PAS des ns ; 0,1 ns par image (FRAMESTEP API) ; 100 ns par replique",
                "licence": "indeterminee (LICENSE: null dans l'API) : usage non commercial, ne pas republier les trajectoires",
                "verifie_le": date,
            }
            rec["promesse"] = ("machine=%s . dossier=%s . %s . verifie le %s" % (
                s4.MACHINE, base, " + ".join("%d .%s" % (v, k) for k, v in rec["compte_fichiers"].items() if k != "total")
                + " (%d fichiers)" % rec["compte_fichiers"]["total"], date))
            ecrire(rec, pdb + "/" + acc)

        # ---------------- systemes v1 qui sortent ----------------
        for P in inv.values():
            if P["statut_vs_v1"] != "sort":
                continue
            key = P["pdb"] + "/" + P["accession"]
            derniere = existants[key][-1]
            rec = json.loads(json.dumps(derniere))
            for k in ("digest", "remplace_digest", "ecrit_le"):
                rec.pop(k, None)
            rec.update({"schema": "corpus-bioemu2-dynabench/2", "criteres_version": "v2-commensurabilite", "statut": "retire_v2", "statut_vs_v1": "sort",
                        "motif_retrait": "sort de la selection v2 (30 places) : les 12 de staging et les VHH passent desormais le critere de commensurabilite et ont "
                                         "priorite ; ce systeme PASSE toujours les criteres v2 (passe_les_criteres_v2=%s) mais est au-dela de la 30e place par la regle "
                                         "« plus petits d'abord ». Ses fichiers restent sur le disque, verifies, mais ne font plus partie du corpus." % P["passe_les_criteres_v2"],
                        "verifie_le": date})
            for f in rec["fichiers"]:
                h = sha256_file(f["chemin"])
                f["sha256_egal_au_manifeste_precedent"] = (h == f["sha256"])
                f["sha256"] = h
                f["verifie_sur_disque_le"] = date
            rec["promesse"] = ("RETIRE du corpus v2 . machine=%s . dossier=%s . fichiers conserves (%d) . verifie le %s" % (
                s4.MACHINE, rec["dossier"], rec["compte_fichiers"]["total"], date))
            ecrire(rec, key)
    print("BILAN manifeste_v2 lignes_ecrites=%d sautees=%d fichier=%s" % (n_ecrit, n_saute, s4.MANIFESTE), flush=True)
    print("FINI s4b %s" % now_iso(), flush=True)


if __name__ == "__main__":
    main()
