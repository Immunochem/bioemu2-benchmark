#!/usr/bin/env python3
"""homologie_B.py — une strate « tenue a l ecart » par la DATE ne l est pas necessairement par la SEQUENCE. Pour chaque systeme PRET de B : pour chacune des 2 chaines, les entrees PDB DEPOSEES
AVANT 2021-10-01 qui portent un homologue de cette chaine (recherche de sequence RCSB, mmseqs2, e-value <= 1, identite >= seuil), pour 3 seuils (0,9 / 0,5 / 0,3). Une entree ANTERIEURE qui porte un homologue
des DEUX chaines est une « paire homologue anterieure » : le modele a pu voir cette paire. Rend les comptes et les 5 premiers identifiants ; ecrit strate_B/homologie.json. API publique, une requete a la fois (pause 1 s)."""
import json, sys, time, urllib.request
from pathlib import Path
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
CUT = "2021-09-30T23:59:59Z"
SEUILS = (0.9, 0.5, 0.3)


def post(obj):
    req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query", data=json.dumps(obj).encode(), headers={"Content-Type": "application/json"})
    for essai in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                txt = r.read()
                return json.loads(txt) if txt.strip() else {"total_count": 0, "result_set": []}
        except Exception as e:
            err = e
            time.sleep(5 * (essai + 1))
    raise err


def entrees_anterieures(seq, ident):
    nodes = [{"type": "terminal", "service": "sequence", "parameters": {"evalue_cutoff": 1, "identity_cutoff": ident, "sequence_type": "protein", "value": seq}},
             {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_accession_info.deposit_date", "operator": "less_or_equal", "value": CUT}}]
    out, start, total = set(), 0, None
    while True:
        q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes}, "return_type": "entry", "request_options": {"paginate": {"start": start, "rows": 1000}}}
        r = post(q)
        total = r.get("total_count", 0)
        out |= {x["identifier"] for x in r.get("result_set", [])}
        start += 1000
        time.sleep(1.0)
        if start >= total or start >= 10000:
            break
    return out, total


m = json.loads((R / "strate_B/manifeste.json").read_text())
res, erreurs = {}, {}
for pid, e in sorted(m["systemes"].items()):
    if e["etat"] != "PRET":
        continue
    ent = {"depot": e["depot"], "tokens": e["tokens"], "titre": e["titre"][:90]}
    try:
        for ident in SEUILS:
            sa, ta = entrees_anterieures(e["sequences"][0], ident)
            sb, tb = entrees_anterieures(e["sequences"][1], ident)
            both = sorted(sa & sb)
            ent[f"id{int(ident * 100)}"] = {"entrees_chaine_A": len(sa), "entrees_chaine_B": len(sb), "total_count_A": ta, "total_count_B": tb, "entrees_avec_les_deux": len(both), "exemples": both[:5],
                                           "tronque_A": ta > len(sa), "tronque_B": tb > len(sb)}
        res[pid] = ent
        print(pid, e["depot"], " · ".join(f"id>={int(i * 100)}%: A {ent[f'id{int(i * 100)}']['entrees_chaine_A']} B {ent[f'id{int(i * 100)}']['entrees_chaine_B']} PAIRES {ent[f'id{int(i * 100)}']['entrees_avec_les_deux']}" for i in SEUILS), flush=True)
    except Exception as ex:
        erreurs[pid] = f"{type(ex).__name__}: {str(ex)[:300]}"
        print("ERREUR", pid, erreurs[pid], flush=True)
(R / "strate_B/homologie.json").write_text(json.dumps({"coupure_depot": CUT, "seuils_identite": SEUILS, "systemes": res, "erreurs": erreurs}, indent=1, ensure_ascii=False) + "\n")
for s in SEUILS:
    k = f"id{int(s * 100)}"
    n = sum(1 for v in res.values() if v[k]["entrees_avec_les_deux"] > 0)
    print(f"BILAN_HOMOLOGIE identite >= {int(s * 100)} % : {n} systemes B sur {len(res)} ont une paire homologue DEPOSEE avant la coupure ; erreurs {len(erreurs)} (REJETS = {len(erreurs)})", flush=True)
sys.exit(0 if not erreurs else 1)
