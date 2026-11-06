#!/usr/bin/env python3
"""construire_taches.py <sortie.json> <pos|neg> <strate> [<strate> ...] [--max-tokens N] : liste des taches GPU (systemes PRETS des manifestes), triee par taille croissante.
Pour `neg`, la tache porte `site_chaine` et `site_residus` (test de SITE). Stdlib seulement. Imprime le compte et la taille cumulee."""
import json, sys
from pathlib import Path
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
args = sys.argv[1:]
mx = None
if "--max-tokens" in args:
    i = args.index("--max-tokens"); mx = int(args[i + 1]); del args[i:i + 2]
conds = None
if "--conds" in args:
    i = args.index("--conds"); conds = [[c.split(":")[0], float(c.split(":")[1])] for c in args[i + 1].split(",")]; del args[i:i + 2]
ids_only = None
if "--ids" in args:
    i = args.index("--ids"); ids_only = set(args[i + 1].split(",")); del args[i:i + 2]
niveaux = None
if "--niveaux" in args:
    i = args.index("--niveaux"); niveaux = {int(x) for x in args[i + 1].split(",")}; del args[i:i + 2]
hom = json.loads((R / "strate_B/homologie.json").read_text())["systemes"] if niveaux is not None else {}
sources_neg = None
if "--sources-neg" in args:       # les systemes dont la chaine A sert au temoin negatif NEG_A (comparaison appariee vrai partenaire / etranger)
    i = args.index("--sources-neg"); del args[i:i + 1]
    sources_neg = {e["source_A"] for e in json.loads((R / "strate_NEG_A/manifeste.json").read_text())["systemes"].values() if e["etat"] == "PRET"}
sortie, kind, strates = args[0], args[1], args[2:]
T = []
for nom in strates:
    m = json.loads((R / f"strate_{nom}/manifeste.json").read_text())
    for pid, e in sorted(m.get("systemes", m).items()):
        if e["etat"] != "PRET":
            continue
        tok = e.get("tokens") or (sum(len(s) for s in e["sequences"]) if e.get("sequences") else None)
        if mx and tok and tok > mx:
            continue
        if ids_only and pid not in ids_only:
            continue
        if sources_neg is not None and pid not in sources_neg:
            continue
        if niveaux is not None:
            h = hom[pid]
            nv = 3 if h["id90"]["entrees_avec_les_deux"] > 0 else 2 if h["id50"]["entrees_avec_les_deux"] > 0 else 1 if h["id30"]["entrees_avec_les_deux"] > 0 else 0
            if nv not in niveaux:
                continue
        t = {"id": pid, "strate": nom, "kind": kind, "x0": str(R / f"strate_{nom}/{pid}/x0.pdb"), "tokens": tok}
        if conds:
            t["conds"] = conds
        if kind == "neg":
            t.update(site_chaine=e["site_chaine"], site_residus=e["site_residus"])
        T.append(t)
T.sort(key=lambda t: (t["tokens"] or 0, t["id"]))
Path(sortie).write_text(json.dumps(T, indent=1) + "\n")
print(f"{len(T)} taches {kind} dans {sortie} · tokens min {min(t['tokens'] for t in T) if T else None} max {max(t['tokens'] for t in T) if T else None} somme {sum(t['tokens'] or 0 for t in T)}")
