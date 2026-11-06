#!/usr/bin/env python3
"""construire_attente2.py <sortie.json> — file d attente MSA ordonnee par CE QUE CHAQUE SYSTEME DECIDE (consigne du coordinateur, 2026-11-06) :
  1. B niveau 0 : AUCUNE entree anterieure a la coupure ne porte un homologue (>= 30 %) des deux chaines = les seuls inedits (toutes tailles : ce sont eux qui repondent a memoire / recherche)
  2. NEG_A (le temoin negatif)
  2b. les systemes (<= 500 residus) dont la chaine A sert au temoin negatif : mesure POSITIVE apparie au negatif\n  3. A_corpus restants (MSA en cache pour 8 sur 10)
  4. B niveau 1, 2, 3 (<= 500 residus) par ressemblance CROISSANTE a une paire deja vue : niveau 1 = paire homologue a 30-50 %, niveau 2 = 50-90 %, niveau 3 = >= 90 % -> le GRADIENT
  5. NEG_B
  6. A_apparie (<= 500 residus d abord, puis le reste)
Niveau d un systeme B = le plus haut seuil d identite (30 / 50 / 90 %) pour lequel une entree anterieure porte un homologue des deux chaines (strate_B/homologie.json). Stdlib seulement."""
import collections, json, sys
from pathlib import Path
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
sortie = sys.argv[1]
mB = json.loads((R / "strate_B/manifeste.json").read_text())["systemes"]
hom = json.loads((R / "strate_B/homologie.json").read_text())["systemes"]


def niveau(pid):
    h = hom[pid]
    if h["id90"]["entrees_avec_les_deux"] > 0:
        return 3
    if h["id50"]["entrees_avec_les_deux"] > 0:
        return 2
    if h["id30"]["entrees_avec_les_deux"] > 0:
        return 1
    return 0


def prets(nom):
    m = json.loads((R / f"strate_{nom}/manifeste.json").read_text())
    S = m.get("systemes", m)
    return [(e.get("tokens") or sum(len(s) for s in e["sequences"]), pid, e["sequences"]) for pid, e in sorted(S.items()) if e["etat"] == "PRET"]


B = [(t, pid, seqs, niveau(pid)) for t, pid, seqs in prets("B")]
file_ = []
for t, pid, seqs, nv in sorted(B):
    if nv == 0:
        file_.append((f"B niveau 0", pid, t, seqs))
for t, pid, seqs in sorted(prets("NEG_A")):
    file_.append(("NEG_A", pid, t, seqs))
# les systemes dont la chaine A sert au temoin negatif (NEG_A) : leur mesure POSITIVE (vrai partenaire) permet la comparaison APPARIEE vrai partenaire / etranger, chaine A par chaine A
mneg = json.loads((R / "strate_NEG_A/manifeste.json").read_text())["systemes"]
src = {e["source_A"] for e in mneg.values() if e["etat"] == "PRET"}
for nom in ("A_apparie", "A_corpus"):
    for t, pid, seqs in sorted(prets(nom)):
        if pid in src and t <= 500:
            file_.append((f"{nom} (source A de NEG_A)", pid, t, seqs))
for t, pid, seqs in sorted(prets("A_corpus")):
    file_.append(("A_corpus", pid, t, seqs))
for nv in (1, 2, 3):
    for t, pid, seqs, n_ in sorted(B):
        if n_ == nv and t <= 500:
            file_.append((f"B niveau {nv}", pid, t, seqs))
for t, pid, seqs in sorted(prets("NEG_B")):
    file_.append(("NEG_B", pid, t, seqs))
A = sorted(prets("A_apparie"))
for t, pid, seqs in [a for a in A if a[0] <= 500] + [a for a in A if a[0] > 500]:
    file_.append(("A_apparie", pid, t, seqs))
vus, out = set(), []
for etiq, pid, t, seqs in file_:
    if pid in vus:
        continue
    vus.add(pid); out.append({"id": pid, "strate": etiq, "tokens": t, "sequences": seqs})
Path(sortie).write_text(json.dumps(out, indent=1) + "\n")
niv = collections.Counter(niveau(p) for _, p, _, _ in [(0, b[1], 0, 0) for b in B])
print(f"{len(out)} systemes dans la file MSA {sortie} : {dict(collections.Counter(o['strate'] for o in out))}")
print("niveaux d homologie des 30 B :", dict(sorted(niv.items())), "· niveau 0 :", [(b[1], b[0]) for b in sorted(B) if b[3] == 0])
print("ordre des 14 premiers :", [(o["id"], o["strate"]) for o in out[:14]])
