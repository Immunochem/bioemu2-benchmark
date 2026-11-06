#!/usr/bin/env python3
"""confusions.py — CPU, lecture de resultats/images.csv. Le contraste A_apparie / B (et les 6 inedits) est-il porte par autre chose que la date ? Covariables par systeme : residus, contacts natifs, d_COM natif,
profondeur du MSA, derive C0. Statistique : taux de succes par image (classe CAPRI >= acceptable) pooles sur les 4 bras SANS orientation donnee (R +10/+20/+30, RD +20) et, a part, sur les 3 bras avec orientation (T +10/+20/+30).
Sorties : correlations de Spearman par strate, comparaison A_apparie / B dans des tranches de covariable, bootstrap sur les SYSTEMES. Ecrit resultats/confusions_lag1000.md/.json."""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
rng = np.random.default_rng(20261106)
rows = list(csv.DictReader(open(R / "resultats/images.csv")))
rows = [r for r in rows if r["kind"] == "pos" and r["strate"] in ("A_apparie", "A_corpus", "B")]
ROT, TRA = ("R_d10", "R_d20", "R_d30", "RD_d20"), ("T_d10", "T_d20", "T_d30")
niv = {}
h = json.loads((R / "strate_B/homologie.json").read_text())["systemes"]
for pid, v in h.items():
    niv[pid] = 3 if v["id90"]["entrees_avec_les_deux"] else 2 if v["id50"]["entrees_avec_les_deux"] else 1 if v["id30"]["entrees_avec_les_deux"] else 0


def profondeur(pid):
    f = R / f"msa/{pid}.json"
    d = json.loads(f.read_text())
    return min(c["non_apparie"] for c in d["chaines"].values())


S = {}
for r in rows:
    s = S.setdefault(r["id"], {"strate": r["strate"], "tokens": float(r["tokens"]), "contacts": float(r["n_contacts_natifs"]), "d_nat": float(r["d_nat_A"]), "rot": [0, 0], "tra": [0, 0], "c0": []})
    ok = r["classe_capri"] != "incorrect"
    if r["cond"] in ROT:
        s["rot"][0] += ok; s["rot"][1] += 1
    if r["cond"] in TRA:
        s["tra"][0] += ok; s["tra"][1] += 1
    if r["cond"] == "C0_d0":
        s["c0"].append(max(float(r["derive_A_ca_A"]), float(r["derive_B_ca_A"])))
for i, s in S.items():
    s["msa"] = profondeur(i); s["derive_c0"] = float(np.median(s["c0"])) if s["c0"] else float("nan"); s["niveau"] = niv.get(i)
    s["taux_rot"] = s["rot"][0] / s["rot"][1] if s["rot"][1] else float("nan")
    s["taux_tra"] = s["tra"][0] / s["tra"][1] if s["tra"][1] else float("nan")


def pool(ids, k):
    return np.array([[S[i][k][0], S[i][k][1]] for i in ids if S[i][k][1]], float)


def taux(P):
    return P[:, 0].sum() / P[:, 1].sum() if len(P) else float("nan")


def boot(P, n=4000):
    if len(P) == 0:
        return [float("nan")] * 2
    ix = rng.integers(0, len(P), (n, len(P)))
    r = P[ix][:, :, 0].sum(1) / P[ix][:, :, 1].sum(1)
    return [float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))]


def bdiff(Pa, Pb, n=4000):
    if len(Pa) == 0 or len(Pb) == 0:
        return [float("nan")] * 3
    ia, ib = rng.integers(0, len(Pa), (n, len(Pa))), rng.integers(0, len(Pb), (n, len(Pb)))
    d = Pa[ia][:, :, 0].sum(1) / Pa[ia][:, :, 1].sum(1) - Pb[ib][:, :, 0].sum(1) / Pb[ib][:, :, 1].sum(1)
    return [float(taux(Pa) - taux(Pb)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


A = sorted(i for i, s in S.items() if s["strate"] == "A_apparie")
B = sorted(i for i, s in S.items() if s["strate"] == "B")
md = ["# Confusions possibles du contraste A_apparie / B (lag 1000 ns)\n",
      f"A_apparie {len(A)} systemes, B {len(B)} systemes. Taux de succes par image (classe CAPRI >= acceptable), pooles sur les bras SANS orientation (R +10/+20/+30, RD +20 : 12 images par systeme) ou AVEC (T +10/+20/+30 : 9 images).\n"]
out = {"n_A": len(A), "n_B": len(B)}
md.append("## 1. Correlations de Spearman systeme par systeme (taux des bras sans orientation)\n")
md.append("| covariable | A_apparie : rho (p) | B : rho (p) |\n|---|---|---|")
for k in ("tokens", "contacts", "d_nat", "msa", "derive_c0"):
    cells = []
    for ids_ in (A, B):
        x = [S[i][k] for i in ids_]; y = [S[i]["taux_rot"] for i in ids_]
        ok = [j for j in range(len(x)) if np.isfinite(x[j])]
        rho, p = spearmanr([x[j] for j in ok], [y[j] for j in ok])
        cells.append(f"{rho:+.2f} ({p:.3f})")
        out.setdefault("spearman", {})[f"{k}|{'A' if ids_ is A else 'B'}"] = [float(rho), float(p)]
    md.append(f"| {k} | " + " | ".join(cells) + " |")
md.append("\n## 2. A_apparie contre B dans des tranches de covariable (bras SANS orientation, puis AVEC)\n")
md.append("Difference = taux A_apparie − taux B, [IC95 bootstrap sur les systemes] ; entre parentheses : systemes A / B.\n")
md.append("| tranche | sans orientation : A | B | difference | avec orientation : A | B | difference |\n|---|---|---|---|---|---|---|")


def ligne(nom, fa, fb=None):
    ia = [i for i in A if fa(S[i])]
    ib = [i for i in B if (fb or fa)(S[i])]
    cells = [nom]
    for k in ("rot", "tra"):
        Pa, Pb = pool(ia, k), pool(ib, k)
        d = bdiff(Pa, Pb)
        out.setdefault("tranches", {})[f"{nom}|{k}"] = {"n_A": len(Pa), "n_B": len(Pb), "taux_A": float(taux(Pa)) if len(Pa) else None, "taux_B": float(taux(Pb)) if len(Pb) else None, "diff": d}
        cells += [f"{taux(Pa):.2f} ({len(Pa)})" if len(Pa) else "—", f"{taux(Pb):.2f} ({len(Pb)})" if len(Pb) else "—", f"{d[0]:+.2f} [{d[1]:+.2f};{d[2]:+.2f}]" if np.isfinite(d[0]) else "—"]
    md.append("| " + " | ".join(cells) + " |")


ligne("tous", lambda s: True)
ligne("moins de 250 residus", lambda s: s["tokens"] < 250)
ligne("250 a 349 residus", lambda s: 250 <= s["tokens"] < 350)
ligne("350 residus et plus", lambda s: s["tokens"] >= 350)
ligne("contacts natifs <= 70", lambda s: s["contacts"] <= 70)
ligne("contacts natifs > 70", lambda s: s["contacts"] > 70)
ligne("d_COM natif >= 20 A", lambda s: s["d_nat"] >= 20)
ligne("d_COM natif < 20 A", lambda s: s["d_nat"] < 20)
ligne("MSA >= 100 et derive C0 <= 3 A", lambda s: s["msa"] >= 100 and s["derive_c0"] <= 3.0)
ligne("MSA >= 1000 et derive C0 <= 3 A", lambda s: s["msa"] >= 1000 and s["derive_c0"] <= 3.0)
ligne("MSA >= 100, derive <= 3 A, 150-400 residus, contacts <= 110", lambda s: s["msa"] >= 100 and s["derive_c0"] <= 3.0 and 150 <= s["tokens"] <= 400 and s["contacts"] <= 110)
md.append("\n## 3. Les 12 premiers A_apparie (chaines A des paires negatives) contre les 12 suivants\n")
src = set()
for nid, e in json.loads((R / "strate_NEG_A/manifeste.json").read_text())["systemes"].items():
    if e["etat"] == "PRET":
        src.add(e["source_A"])
A1, A2 = [i for i in A if i in src], [i for i in A if i not in src]
md.append("| sous-ensemble | systemes | residus (mediane) | contacts natifs (mediane) | d_COM natif (mediane) | taux sans orientation | taux avec orientation |\n|---|---|---|---|---|---|---|")
for nom, ids_ in (("A_apparie : chaines A des paires negatives", A1), ("A_apparie : autres", A2), ("B", B)):
    pr, pt = pool(ids_, "rot"), pool(ids_, "tra")
    md.append(f"| {nom} | {len(ids_)} | {np.median([S[i]['tokens'] for i in ids_]):.0f} | {np.median([S[i]['contacts'] for i in ids_]):.0f} | {np.median([S[i]['d_nat'] for i in ids_]):.1f} | "
              f"{taux(pr):.2f} [{boot(pr)[0]:.2f};{boot(pr)[1]:.2f}] | {taux(pt):.2f} [{boot(pt)[0]:.2f};{boot(pt)[1]:.2f}] |")
    out.setdefault("A1_A2_B", {})[nom] = {"n": len(ids_), "taux_rot": float(taux(pr)), "taux_tra": float(taux(pt))}
(R / "resultats/confusions_lag1000.md").write_text("\n".join(md) + "\n")
(R / "resultats/confusions_lag1000.json").write_text(json.dumps(out, indent=1) + "\n")
print("\n".join(md))
