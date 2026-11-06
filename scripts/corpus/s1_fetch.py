"""J1 — un seul passage reseau pour l'inventaire.

  (a) structure.pdb de la replique 1 de chacun des 210 projets (point d'acces FICHIERS ;
      la LISTE des projets n'est PAS reinterrogee : on lit le cache L0a) ;
  (b) UNE requete GraphQL RCSB par paquet de 30 entrees : titre, methode, resolution,
      entites (description, organisme, sequence), instances (chaine auteur,
      correspondance auteur -> position d'entite, residus NON observes).

Tout est mis en cache sous corpus/cache/. Idempotent : un passage deja complet n'est pas refait.
Un echec est une ligne du resultat, jamais un saut.
"""
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import (API, CACHE, ROOT, declared_length, http_get, jdump,  # noqa: E402
                        charger_projets_dynabench, now_iso)

os.makedirs(CACHE + "/structure_rep1", exist_ok=True)

GQL_Q = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info { resolution_combined polymer_entity_count_protein polymer_entity_count_nucleic_acid
                      deposited_polymer_entity_instance_count experimental_method }
    refine { ls_d_res_high }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      entity_poly { rcsb_entity_polymer_type pdbx_seq_one_letter_code_can }
      rcsb_polymer_entity_container_identifiers { entity_id asym_ids auth_asym_ids }
      rcsb_entity_source_organism { ncbi_scientific_name }
      polymer_entity_instances {
        rcsb_id
        rcsb_polymer_entity_instance_container_identifiers { entity_id asym_id auth_asym_id auth_to_entity_poly_seq_mapping }
        rcsb_polymer_instance_feature { type name feature_positions { beg_seq_id end_seq_id } }
      }
    }
  }
}"""


def gql(ids, tries=3):
    body = json.dumps({"query": GQL_Q, "variables": {"ids": ids}}).encode()
    last = None
    for t in range(1, tries + 1):
        t0 = time.time()
        try:
            req = urllib.request.Request("https://data.rcsb.org/graphql", data=body,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                txt = r.read()
                return {"http": r.status, "secondes": round(time.time() - t0, 3), "essai": t,
                        "octets": len(txt), "reponse": json.loads(txt), "erreur": None}
        except Exception as e:  # noqa: BLE001
            last = {"http": getattr(e, "code", None), "secondes": round(time.time() - t0, 3),
                    "essai": t, "octets": 0, "reponse": None, "erreur": "%s: %s" % (type(e).__name__, e)}
            time.sleep(3 * t)
    return last


def main():
    ld = charger_projets_dynabench()
    print("cache projets: %s  sha256=%s  total=%d  dynabench=%d" % (
        ld["fichier"], ld["sha256"], ld["n_total_cache"], ld["n_dynabench"]), flush=True)
    if ld["n_dynabench"] != 210:
        print("ARRET: %d projets DynaBench au lieu de 210" % ld["n_dynabench"], flush=True)
        sys.exit(3)
    projs = ld["projets"]

    # ---------------- (a) structure.pdb replique 1 ----------------
    out_jsonl = CACHE + "/fetch_structures.jsonl"
    deja = {}
    if os.path.exists(out_jsonl):
        for line in open(out_jsonl, encoding="utf-8"):
            r = json.loads(line)
            deja[r["accession"]] = r
    n_ok = n_ko = 0
    t_tot0 = time.time()
    octets_tot = 0
    with open(out_jsonl, "a", encoding="utf-8") as fo:
        for q in projs:
            acc = q["accession"]
            pdb = q["metadata"]["PDBIDS"][0]
            if acc in deja and deja[acc].get("conforme"):
                n_ok += 1
                continue
            dest = "%s/structure_rep1/%s_%s.pdb" % (CACHE, pdb, acc)
            url = "%s/projects/%s.1/files/structure.pdb" % (API, acc)
            r = http_get(url, dest=dest, timeout=300)
            decl = declared_length(acc, 1, "structure.pdb")
            r.update({"accession": acc, "pdb": pdb, "mdName_1": q["mds"][0] if q.get("mds") else None,
                      "declare_serveur": decl, "dest": dest, "date": now_iso()})
            r["conforme"] = (r["status"] == 200 and decl is not None and r["octets"] == decl)
            fo.write(json.dumps(r, ensure_ascii=False) + "\n")
            fo.flush()
            octets_tot += r["octets"] or 0
            if r["conforme"]:
                n_ok += 1
            else:
                n_ko += 1
                print("ECHEC structure %s %s status=%s octets=%s declare=%s erreur=%s" % (
                    pdb, acc, r["status"], r["octets"], decl, r["erreur"]), flush=True)
            time.sleep(0.15)
    dt = time.time() - t_tot0
    print("structures: conformes=%d  non_conformes=%d  octets_telecharges=%d  secondes=%.1f  debit=%.2f Mo/s" % (
        n_ok, n_ko, octets_tot, dt, (octets_tot / 1e6 / dt) if dt > 0 and octets_tot else 0.0), flush=True)

    # ---------------- (b) RCSB, un seul passage ----------------
    rc_path = CACHE + "/rcsb_graphql.json"
    if os.path.exists(rc_path):
        print("RCSB: cache deja present (%s), on ne reinterroge pas" % rc_path, flush=True)
    else:
        ids = sorted({q["metadata"]["PDBIDS"][0] for q in projs})
        chunks = []
        for i in range(0, len(ids), 30):
            part = ids[i:i + 30]
            r = gql(part)
            got = []
            if r["reponse"] and r["reponse"].get("data") and r["reponse"]["data"].get("entries"):
                got = [e["rcsb_id"] for e in r["reponse"]["data"]["entries"] if e]
            r.update({"ids_demandes": part, "ids_recus": got,
                      "ids_absents": sorted(set(part) - set(got))})
            chunks.append(r)
            print("RCSB paquet %d: demandes=%d recus=%d absents=%s http=%s erreur=%s graphql_errors=%s" % (
                i // 30 + 1, len(part), len(got), r["ids_absents"], r["http"], r["erreur"],
                json.dumps(r["reponse"].get("errors"))[:200] if r["reponse"] else None), flush=True)
            time.sleep(0.5)
        jdump({"date": now_iso(), "source": "https://data.rcsb.org/graphql", "requete": GQL_Q,
               "n_ids": len(ids), "chunks": chunks}, rc_path)
        print("RCSB: cache ecrit %s" % rc_path, flush=True)
    print("FINI s1 %s" % now_iso(), flush=True)


if __name__ == "__main__":
    main()
