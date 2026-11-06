#!/usr/bin/env python3
"""tokens_vram.py — les 28 classes par nombre de TOKENS, contre le plafond de VRAM MESURE.

Tokens = nombre de residus de la sequence donnee a BioEmu2 (un token par residu standard), lu dans les x0 cures
de l'audit (`audit_entrees.json`, `longueurs_chaines`), pas recompte. Deux mesures du RESERVE (grandeur qui
declenche l'OOM) donnees par le coordinateur : 390 tokens -> 11,36 Gio (6CWG, smoke) ; 777 tokens -> 22,08 Gio
(2AJF, bf16-mixed, un seul job sur la slice 1g.24gb de 23,62 Gio ; fp32 : OOM a 21,86 Gio alloues).
L'exposant et les reserves des AUTRES systemes sont des PROJECTIONS depuis ces deux points, pas des mesures.
"""
import json
import math
from pathlib import Path

SORTIE = Path("/gpfs/phherman/bioemu2_bench/prep_fix/curation")
SLICE_GIO = 23.62
P1, P2 = (390, 11.36), (777, 22.08)
EXPOSANT = math.log(P2[1] / P1[1]) / math.log(P2[0] / P1[0])
SEUIL_TOKENS = 600


def reserve(n):
    return P1[1] * (n / P1[0]) ** EXPOSANT


def main():
    R = json.load(open(SORTIE / "audit_entrees.json"))
    tok = {}
    for r in R:
        if r["lot"] == "perimetre" and r.get("kind") == "img0" and r.get("replique") == 1:
            if r["etat"] != "PRET":
                raise SystemExit(f"🔴 ARRET : {r['etiquette']} non PRET")
            tok[r["pdb"]] = (sum(r["longueurs_chaines"]), r["longueurs_chaines"])
    inp = {r["pdb"]: sum(r["longueurs_chaines"]) for r in R if r["lot"] == "perimetre" and r.get("kind") == "input" and r["etat"] == "PRET"}
    for pdb, n in inp.items():
        if tok[pdb][0] != n:
            raise SystemExit(f"🔴 ARRET : {pdb} : {tok[pdb][0]} tokens (image 0) contre {n} (input)")
    n90 = P1[0] * (0.9 * SLICE_GIO / P1[1]) ** (1 / EXPOSANT)
    n100 = P1[0] * (SLICE_GIO / P1[1]) ** (1 / EXPOSANT)
    print(f"exposant du RESERVE depuis les deux mesures : ln({P2[1]}/{P1[1]})/ln({P2[0]}/{P1[0]}) = {EXPOSANT:.3f}")
    print(f"tokens ou le reserve projete atteint 90 % de la slice ({0.9 * SLICE_GIO:.2f} Gio) : {n90:.0f} ; 100 % ({SLICE_GIO} Gio) : {n100:.0f} ; a 1030 tokens : {reserve(1030):.1f} Gio")
    print(f"{'systeme':<8}{'tokens':>7}  chaines            reserve Gio  % slice  statut")
    n_plus = []
    for pdb, (n, ch) in sorted(tok.items(), key=lambda kv: kv[1][0]):
        rv = reserve(n)
        if n in (P1[0], P2[0]):
            st = "MESURE"
        elif P1[0] < n < P2[0]:
            st = "interpole entre deux mesures"
        elif n < P1[0]:
            st = "projete (sous la plus petite mesure)"
        else:
            st = "EXTRAPOLE"
        marque = "  >~600 : UNE seule unite GPU a la fois sur la slice" if n > SEUIL_TOKENS else ""
        if n > SEUIL_TOKENS:
            n_plus.append(pdb)
        print(f"{pdb:<8}{n:>7}  {str(ch):<18} {rv:>8.2f}  {100 * rv / SLICE_GIO:>6.1f}%  {st}{marque}")
    print(f"\nsystemes au-dela de {SEUIL_TOKENS} tokens : {n_plus} ({len(n_plus)}/28) ; <= {SEUIL_TOKENS} : {28 - len(n_plus)}/28")
    print("tokens de l'image 0 == tokens de l'input pour tous les systemes qui ont un input :", True)
    (SORTIE / "tokens_vram.json").write_text(json.dumps({"exposant_reserve": EXPOSANT, "mesures": {"390": 11.36, "777": 22.08}, "slice_gio": SLICE_GIO,
                                                         "tokens": {k: v[0] for k, v in tok.items()}, "au_dela_de_600": n_plus}, indent=1) + "\n")


if __name__ == "__main__":
    main()
