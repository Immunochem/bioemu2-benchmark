#!/usr/bin/env python3
"""choisir_A.py — strate A : pour chacun des 28 systemes du corpus, choisir les DEUX chaines du cristal RCSB qui portent l interface (contacts atomes lourds < 5 A
entre residus, maximum parmi les copies de l unite asymetrique), les extraire (extraire_cristal.py) et ecrire le bilan, REFUS compris. Aucun refus n est silencieux."""
import json, sys, itertools, collections
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/redocking/scripts")
import numpy as np
from scipy.spatial import cKDTree
import extraire_cristal as E
B = Path("/gpfs/phherman/bioemu2_bench"); R = B / "redocking"
HORS = {"1ZM4", "1FQ1"}
man = {x["pdb"]: x for x in (json.loads(l) for l in open(B / "corpus/manifeste_corpus.json") if l.strip()) if x.get("statut") == "actif_v2" and x["pdb"] not in HORS}
assert len(man) == 28, len(man)
cro = json.load(open(B / "generation/cristal/croisement_cristal_psf.json"))
def contacts(ra, rb, seuil=5.0):
    xa = [(i, a[2:5]) for i, r in enumerate(ra["residus"]) for a in r["atomes"]]
    xb = [(j, a[2:5]) for j, r in enumerate(rb["residus"]) for a in r["atomes"]]
    if not xa or not xb: return 0
    ta, tb = cKDTree(np.array([p for _, p in xa])), cKDTree(np.array([p for _, p in xb]))
    s = set()
    for i, l in enumerate(ta.query_ball_tree(tb, seuil)):
        for j in l: s.add((xa[i][0], xb[j][0]))
    return len(s)
bilan = {}
for pdb in sorted(man):
    cif = B / f"generation/cristal/{pdb}.cif"
    cv = cro[pdb]["chaines_cristal_vers_md"]
    md_chaines = [c["id"] for c in man[pdb]["chaines"]]
    cand = {m: [k for k, v in cv.items() if v["chaine_md_origine"] == m] for m in md_chaines}
    essais = []
    for ca, cb in itertools.product(cand[md_chaines[0]], cand[md_chaines[1]]):
        res, refus = E.extraire(cif, [ca, cb])
        if len(res["chaines"]) < 2:
            essais.append({"chaines": [ca, cb], "refus": refus, "contacts": None}); continue
        essais.append({"chaines": [ca, cb], "refus": [r["classe"] for r in refus], "detail_refus": [r["detail"] for r in refus],
                       "contacts": contacts(res["chaines"][0], res["chaines"][1]), "n_residus": [c["n_residus"] for c in res["chaines"]],
                       "meme_lettre": [cv[ca]["meme_lettre"], cv[cb]["meme_lettre"]]})
    ok = [e for e in essais if e["contacts"] is not None]
    best = max(ok, key=lambda e: (e["contacts"], sum(e["meme_lettre"]))) if ok else None
    bilan[pdb] = {"md_chaines": md_chaines, "candidats": cand, "choix": best, "essais": essais,
                  "interface_paires_5A_md_image0": man[pdb]["interface_paires_5A"]["structure_pdb_image0"],
                  "rcsb": man[pdb]["rcsb"]}
    b = best
    print(pdb, "choix", b["chaines"] if b else None, "contacts", b["contacts"] if b else None, "(MD image0:", man[pdb]["interface_paires_5A"]["structure_pdb_image0"], ")",
          "n_res", b.get("n_residus") if b else None, "refus", b["refus"] if b else None, flush=True)
(R / "strate_A").mkdir(exist_ok=True)
(R / "strate_A/choix_chaines.json").write_text(json.dumps(bilan, indent=1, ensure_ascii=False) + "\n")
c = collections.Counter(tuple(b["choix"]["refus"]) if b["choix"] else ("AUCUN",) for b in bilan.values())
print("BILAN refus par classe (choix retenu) :", dict(c), "n_systemes", len(bilan))
