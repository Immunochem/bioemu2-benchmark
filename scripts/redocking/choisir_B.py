#!/usr/bin/env python3
"""choisir_B.py — strate B (tenue a l ecart) : candidats RCSB deposes APRES 2021-09-30, 2 chaines proteiques / 2 entites / 2 instances,
X-ray <= 3,0 A, 150-835 residus deposes. Interroge l API publique (pas de calcul), ecrit candidats.json (TOUS les candidats + metadonnees).
Le filtre est la DATE DE DEPOT (et non de publication) : deposee apres la coupure => publiee apres, la contrainte la plus stricte."""
import json, sys, time, urllib.request, collections
OUT = "/gpfs/phherman/bioemu2_bench/redocking/strate_B/candidats.json"
CUT = "2021-09-30T23:59:59Z"
def post(url, obj):
    req = urllib.request.Request(url, data=json.dumps(obj).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read())
T = lambda a, o, v: {"type": "terminal", "service": "text", "parameters": {"attribute": a, "operator": o, "value": v}}
q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
    T("rcsb_accession_info.deposit_date", "greater", CUT),
    T("rcsb_entry_info.polymer_entity_count_protein", "equals", 2),
    T("rcsb_entry_info.polymer_entity_count", "equals", 2),
    T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 2),
    T("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
    T("rcsb_entry_info.resolution_combined", "less_or_equal", 3.0),
    T("rcsb_entry_info.deposited_polymer_monomer_count", "range", {"from": 150, "include_lower": True, "to": 835, "include_upper": True})]},
    "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 1000}, "sort": [{"sort_by": "rcsb_accession_info.deposit_date", "direction": "asc"}]}}
ids, total, start = [], None, 0
while True:
    q["request_options"]["paginate"]["start"] = start
    r = post("https://search.rcsb.org/rcsbsearch/v2/query", q)
    total = r["total_count"]; ids += [x["identifier"] for x in r["result_set"]]
    start += 1000
    if start >= total: break
    time.sleep(1)
print("recherche :", len(ids), "identifiants sur total_count", total, flush=True)
GQL = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
 rcsb_accession_info{deposit_date initial_release_date}
 rcsb_entry_info{resolution_combined deposited_polymer_monomer_count deposited_modeled_polymer_monomer_count}
 struct{title}
 polymer_entities{ rcsb_id entity_poly{pdbx_seq_one_letter_code_can rcsb_entity_polymer_type pdbx_strand_id}
   rcsb_polymer_entity_group_membership{aggregation_method group_id similarity_cutoff}
   rcsb_entity_source_organism{ncbi_scientific_name}
   rcsb_polymer_entity{pdbx_description} } } }"""
cand, rej = {}, []
for i in range(0, len(ids), 100):
    lot = ids[i:i + 100]
    r = post("https://data.rcsb.org/graphql", {"query": GQL, "variables": {"ids": lot}})
    if "errors" in r and not r.get("data"):
        rej.append({"lot": lot[0], "erreur": str(r["errors"])[:300]}); continue
    for e in r["data"]["entries"] or []:
        if e is None: continue
        cand[e["rcsb_id"]] = e
    time.sleep(0.5)
manque = [x for x in ids if x not in cand]
json.dump({"coupure": CUT, "n_ids_recherche": len(ids), "n_metadonnees": len(cand), "ids_sans_metadonnees": manque, "rejets_lots": rej, "entrees": cand}, open(OUT, "w"))
par_an = collections.Counter(v["rcsb_accession_info"]["deposit_date"][:4] for v in cand.values())
print("metadonnees :", len(cand), "ids sans metadonnees :", len(manque), "lots rejetes :", len(rej))
print("par annee de depot :", dict(sorted(par_an.items())))
