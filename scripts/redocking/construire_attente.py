#!/usr/bin/env python3
"""construire_attente.py <sortie.json> <strate> [<strate> ...] : file d attente MSA = systemes PRETS des manifestes donnes, dans l ordre des strates donnees ; deux strates consecutives
marquees `+` sont ENTRELACEES (taille croissante dans chaque) pour qu une file interrompue laisse les deux strates equilibrees. Ex. : NEG_A A_corpus A_apparie+B. Stdlib seulement."""
import json, sys
from pathlib import Path
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
sortie, groupes = sys.argv[1], sys.argv[2:]


def prets(nom):
    m = json.loads((R / f"strate_{nom}/manifeste.json").read_text())
    L = [(e.get("tokens") or sum(len(s) for s in e["sequences"]), pid, nom, e["sequences"]) for pid, e in sorted(m.get("systemes", m).items()) if e["etat"] == "PRET"]
    return sorted(L)


file_ = []
for g in groupes:
    noms = g.split("+")
    listes = [prets(n) for n in noms]
    if len(listes) == 1:
        file_ += listes[0]
    else:
        for i in range(max(len(l) for l in listes)):
            for l in listes:
                if i < len(l):
                    file_.append(l[i])
vus, out = set(), []
for tok, pid, nom, seqs in file_:
    if pid in vus:
        continue
    vus.add(pid)
    out.append({"id": pid, "strate": nom, "tokens": tok, "sequences": seqs})
Path(sortie).write_text(json.dumps(out, indent=1) + "\n")
import collections
par = collections.Counter(o["strate"] for o in out)
print(f"{len(out)} systemes dans la file MSA {sortie} : {dict(par)}")
