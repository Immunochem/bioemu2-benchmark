#!/usr/bin/env python3
"""partenaire_contre_etranger.py <lag> <etiquette_neg> [<etiquette_neg> ...] — CPU. COMPARAISON APPARIEE, chaine A par chaine A : la meme chaine A (celle du complexe X) avec
  (a) son VRAI partenaire (systeme positif X, memes conditions R / RD, meme lag)  et  (b) un ETRANGER (paire NEG_<X>A_<Y>B).
Mesure : contacts, rappel du site natif de A (part des residus de A qui touchent le vrai partenaire dans le cristal et que touche la chaine B de la sortie), precision du site.
Si le modele ne faisait que « coller », (a) et (b) auraient les memes contacts ET le meme rappel ; s'il sait ou va le vrai partenaire, (a) a un rappel nettement plus haut avec des contacts du meme ordre.
Lit resultats/site_positifs_lag<L>.json (site_positifs.py) et sorties/<etiquette_neg>/*/lag<L>/mesures.json. Ecrit resultats/partenaire_contre_etranger_lag<L>.md et .json."""
import json, sys
from pathlib import Path
import numpy as np
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
lag = int(sys.argv[1]); etiqs = sys.argv[2:]
rng = np.random.default_rng(20261106)
pos = json.loads((R / f"resultats/site_positifs_lag{lag:04d}.json").read_text())["systemes"]
neg, grp = {}, {}
for nom in ("NEG_A", "NEG_B"):
    f = R / f"strate_{nom}/manifeste.json"
    if f.is_file():
        for nid, e in json.loads(f.read_text())["systemes"].items():
            neg[nid] = e; grp[nid] = nom
lignes = []     # (X, nid, cond, pos_rappel_img[], neg_rappel_img[], pos_contacts[], neg_contacts[], base)
for e in etiqs:
    for d in sorted((R / "sorties" / e).glob("*/")):
        mf = d / "meta.json"
        if not mf.is_file():
            continue
        meta = json.loads(mf.read_text())
        mp = d / f"lag{lag:04d}/mesures.json"
        if meta.get("statut") != "ok" or meta["kind"] != "neg" or not mp.is_file() or meta["id"] not in neg:
            continue
        X = neg[meta["id"]]["source_A"]
        if X not in pos:
            continue
        cond = json.loads((d / "conditions.json").read_text())["audit"]
        mes = json.loads(mp.read_text())["mesures"]
        for c in ("R_d10", "R_d20", "R_d30", "RD_d20"):
            ii = [a["indice_image"] for a in cond if a["accepte"] and a["cond"] == c]
            pp = [r for r in pos[X]["images"] if r["cond"] == c]
            if not ii or not pp:
                continue
            lignes.append({"groupe": grp[meta["id"]], "X": X, "neg": meta["id"], "cond": c, "neg_rappel": [mes["rappel_site"][i] for i in ii], "neg_contacts": [mes["n_contacts"][i] for i in ii],
                           "pos_rappel": [r["rappel_site"] for r in pp], "pos_contacts": [r["n_contacts"] for r in pp], "base": pos[X]["base_hasard_precision"]})


def boot_med_diff(v, n=4000):
    v = np.asarray(v, float)
    if len(v) == 0:
        return [None] * 3
    m = np.median(rng.choice(v, (n, len(v))), axis=1)
    return [float(np.median(v)), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


md = [f"# Vrai partenaire contre etranger, chaine A par chaine A — lag {lag} ns\n",
      "Chaque ligne = une paire (X, etranger) dont X a ete mesure avec son vrai partenaire. Valeurs = mediane sur les 3 images de la condition. Difference = vrai partenaire − etranger, IC95 = bootstrap sur les paires.\n",
      "| condition | n paires | rappel du site : vrai partenaire | rappel : etranger | difference rappel [IC95] | contacts : vrai | contacts : etranger | difference contacts [IC95] | paires ou le vrai a le plus haut rappel |\n|---|---|---|---|---|---|---|---|---|"]
out = {"lag": lag, "lignes": lignes, "resume": {}}
for c in ("R_d10", "R_d20", "R_d30", "RD_d20"):
    L = [x for x in lignes if x["cond"] == c]
    if not L:
        continue
    pr, nr = [np.median(x["pos_rappel"]) for x in L], [np.median(x["neg_rappel"]) for x in L]
    pc, nc = [np.median(x["pos_contacts"]) for x in L], [np.median(x["neg_contacts"]) for x in L]
    dr, dc = boot_med_diff(np.array(pr) - np.array(nr)), boot_med_diff(np.array(pc) - np.array(nc))
    wins = int(sum(a > b for a, b in zip(pr, nr)))
    out["resume"][c] = {"n": len(L), "rappel_vrai_med": float(np.median(pr)), "rappel_etranger_med": float(np.median(nr)), "diff_rappel": dr, "contacts_vrai_med": float(np.median(pc)),
                        "contacts_etranger_med": float(np.median(nc)), "diff_contacts": dc, "victoires_vrai": wins, "base_med": float(np.median([x["base"] for x in L]))}
    md.append(f"| {c} | {len(L)} | {np.median(pr):.2f} | {np.median(nr):.2f} | {dr[0]:+.2f} [{dr[1]:+.2f};{dr[2]:+.2f}] | {np.median(pc):.0f} | {np.median(nc):.0f} | {dc[0]:+.0f} [{dc[1]:+.0f};{dc[2]:+.0f}] | {wins}/{len(L)} |")
md.append("\n## Par groupe : chaines A d AVANT la coupure (NEG_A) contre chaines A POSTERIEURES a la coupure (NEG_B)\n")
md.append("Difference = rappel du site avec le vrai partenaire MOINS rappel avec l etranger, mediane des paires, IC95 = bootstrap sur les paires. « difference des differences » = NEG_A − NEG_B (tirages independants).\n")
md.append("| condition | groupe | paires | rappel vrai | rappel etranger | difference [IC95] | paires ou le vrai gagne | contacts vrai · etranger |\n|---|---|---|---|---|---|---|---|")
out["par_groupe"] = {}
for c in ("R_d10", "R_d20", "R_d30", "RD_d20"):
    dd = {}
    for g in ("NEG_A", "NEG_B"):
        L = [x for x in lignes if x["cond"] == c and x["groupe"] == g]
        if not L:
            continue
        pr, nr = np.array([np.median(x["pos_rappel"]) for x in L]), np.array([np.median(x["neg_rappel"]) for x in L])
        pc, nc = [np.median(x["pos_contacts"]) for x in L], [np.median(x["neg_contacts"]) for x in L]
        d = boot_med_diff(pr - nr); dd[g] = pr - nr
        out["par_groupe"][f"{c}|{g}"] = {"n": len(L), "rappel_vrai": float(np.median(pr)), "rappel_etranger": float(np.median(nr)), "diff": d, "victoires": int((pr > nr).sum()),
                                         "contacts_vrai": float(np.median(pc)), "contacts_etranger": float(np.median(nc))}
        md.append(f"| {c} | {g} | {len(L)} | {np.median(pr):.2f} | {np.median(nr):.2f} | {d[0]:+.2f} [{d[1]:+.2f};{d[2]:+.2f}] | {int((pr > nr).sum())}/{len(L)} | {np.median(pc):.0f} · {np.median(nc):.0f} |")
    if len(dd) == 2:
        a_, b_ = dd["NEG_A"], dd["NEG_B"]
        bs = np.median(rng.choice(a_, (4000, len(a_))), axis=1) - np.median(rng.choice(b_, (4000, len(b_))), axis=1)
        out["par_groupe"][f"{c}|diff_des_diff"] = [float(np.median(a_) - np.median(b_)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
        md.append(f"| {c} | **NEG_A − NEG_B** | | | | **{np.median(a_) - np.median(b_):+.2f} [{np.percentile(bs, 2.5):+.2f};{np.percentile(bs, 97.5):+.2f}]** | | |")
md.append("\nBase au hasard du rappel ~ part du site dans la chaine A (mediane %.2f). Une paire sans contraste (vrai ~ etranger) est une preuve que le modele colle sans savoir ou." % (np.median([x["base"] for x in lignes]) if lignes else float("nan")))
md.append("\nDetail par paire (R_d20) :\n\n| X | etranger | rappel vrai | rappel etranger | contacts vrai | contacts etranger |\n|---|---|---|---|---|---|")
for x in lignes:
    if x["cond"] == "R_d20":
        md.append(f"| {x['X']} | {x['neg']} | {np.median(x['pos_rappel']):.2f} | {np.median(x['neg_rappel']):.2f} | {np.median(x['pos_contacts']):.0f} | {np.median(x['neg_contacts']):.0f} |")
(R / f"resultats/partenaire_contre_etranger_lag{lag:04d}.md").write_text("\n".join(md) + "\n")
(R / f"resultats/partenaire_contre_etranger_lag{lag:04d}.json").write_text(json.dumps(out) + "\n")
print("\n".join(md))
