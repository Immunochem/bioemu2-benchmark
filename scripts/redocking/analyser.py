#!/usr/bin/env python3
"""analyser.py <etiquette> [<etiquette> ...] [--lag 50] — CPU. Relit sorties/<etiquette>/*/ (meta.json, conditions.json, lag<L>/mesures.json), le plancher (resultats/plancher/*.json) et les manifestes,
et ecrit resultats/images.csv (une ligne par image) + resultats/tableaux_lag<L>.md + resultats/tableaux_lag<L>.json.
UNITE D ANALYSE = LE SYSTEME : les K repetitions d une condition sont resumees par leur mediane AVANT toute statistique entre systemes ; IC95 par bootstrap sur les SYSTEMES (4 000 tirages).
Les systemes rejetes (MSA absent, OOM, mesure impossible) sont COMPTES dans chaque table. fnat n est jamais rapporte sans l iRMSD."""
import collections, csv, json, os, sys
from pathlib import Path
sys.dont_write_bytecode = True
import numpy as np

R = Path("/gpfs/phherman/bioemu2_bench/redocking")
args = sys.argv[1:]
LAG = 50
if "--lag" in args:
    i = args.index("--lag"); LAG = int(args[i + 1]); del args[i:i + 2]
ETIQ = args
rng = np.random.default_rng(20261106)
CL = ["incorrect", "acceptable", "medium", "high"]


def boot_med(v, n=4000):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) == 0:
        return [None, None]
    m = np.median(rng.choice(v, (n, len(v))), axis=1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def boot_frac(b, n=4000):
    b = np.asarray(b, float)
    if len(b) == 0:
        return [None, None]
    m = rng.choice(b, (n, len(b))).mean(axis=1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def boot_ratio(pairs, n=4000):
    """taux de succes POOLE par image : somme des succes / somme des images, IC95 par bootstrap sur les SYSTEMES (grappes). pairs = [(succes, images), ...] par systeme."""
    pairs = np.asarray(pairs, float)
    if len(pairs) == 0 or pairs[:, 1].sum() == 0:
        return [None, None, None]
    idx = rng.integers(0, len(pairs), (n, len(pairs)))
    r = pairs[idx][:, :, 0].sum(axis=1) / pairs[idx][:, :, 1].sum(axis=1)
    return [float(pairs[:, 0].sum() / pairs[:, 1].sum()), float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))]


def boot_diff_ratio(pa, pb, n=4000):
    pa, pb = np.asarray(pa, float), np.asarray(pb, float)
    if len(pa) == 0 or len(pb) == 0:
        return [None, None, None]
    ia, ib = rng.integers(0, len(pa), (n, len(pa))), rng.integers(0, len(pb), (n, len(pb)))
    d = pa[ia][:, :, 0].sum(axis=1) / pa[ia][:, :, 1].sum(axis=1) - pb[ib][:, :, 0].sum(axis=1) / pb[ib][:, :, 1].sum(axis=1)
    return [float(pa[:, 0].sum() / pa[:, 1].sum() - pb[:, 0].sum() / pb[:, 1].sum()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


def boot_diff_med(a, b, n=4000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return [None, None, None]
    d = np.median(rng.choice(a, (n, len(a))), axis=1) - np.median(rng.choice(b, (n, len(b))), axis=1)
    return [float(np.median(a) - np.median(b)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


# ─────────────── lecture ───────────────
systemes, rejets_lecture = {}, []
for e in ETIQ:
    for d in sorted((R / "sorties" / e).glob("*/")):
        if not (d / "meta.json").is_file():
            continue
        meta = json.loads((d / "meta.json").read_text())
        mf = d / f"lag{LAG:04d}" / "mesures.json"
        if meta.get("statut") != "ok" or not mf.is_file():
            rejets_lecture.append({"id": meta["id"], "motif": f"statut {meta.get('statut')} ou mesures absentes", "etiquette": e}); continue
        cond = json.loads((d / "conditions.json").read_text())
        mes = json.loads(mf.read_text())["mesures"]
        if "ERREUR_MESURE" in mes:
            rejets_lecture.append({"id": meta["id"], "motif": mes["ERREUR_MESURE"], "etiquette": e}); continue
        systemes[meta["id"]] = {"meta": meta, "cond": cond, "mes": mes, "etiquette": e}
man_B = json.loads((R / "strate_B/manifeste.json").read_text())["systemes"]
hom = {}
if (R / "strate_B/homologie.json").is_file():
    hom = json.loads((R / "strate_B/homologie.json").read_text())["systemes"]
plan = {}
for f in (R / "resultats/plancher").glob("*.json"):
    p = json.loads(f.read_text()); plan[p["id"]] = p

# ─────────────── images ───────────────
lignes = []
for sid, S in systemes.items():
    m, kind, strate = S["mes"], S["meta"]["kind"], S["meta"]["strate"]
    for a in S["cond"]["audit"]:
        if not a["accepte"]:
            continue
        i = a["indice_image"]
        L = {"id": sid, "strate": strate, "kind": kind, "cond": a["cond"], "arm": a["arm"], "delta_A": a["delta_A"], "rep": a["rep"], "tokens": S["meta"]["tokens"],
             "d_nat_A": S["meta"]["d_nat_A"], "n_contacts_natifs": S["meta"]["n_contacts_natifs"], "x0_d_com_A": a["d_com_A"], "x0_dmin_A": a["dmin_inter_A"]}
        for k in ("fnat", "irmsd_A", "lrmsd_A", "dockq", "classe_capri", "classe_dockq", "d_com_A", "n_contacts", "dmin_inter_A", "derive_A_ca_A", "derive_B_ca_A",
                  "rappel_site", "precision_site", "n_residus_touches_chaine_site"):
            if k in m:
                L[k if k != "dmin_inter_A" else "dmin_sortie_A"] = m[k][i]
        lignes.append(L)
cols = sorted({k for L in lignes for k in L})
with open(R / "resultats/images.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(lignes)

# ─────────────── tables par strate x condition (positifs) ───────────────
pos = [L for L in lignes if L["kind"] == "pos"]
strates = [s for s in ("A_apparie", "A_corpus", "B", "special") if any(L["strate"] == s for L in pos)]
conds = []
for L in pos:
    if L["cond"] not in conds:
        conds.append(L["cond"])
ordre = ["C0_d0", "T_d10", "T_d20", "T_d30", "R_d10", "R_d20", "R_d30", "RD_d20"]
conds = [c for c in ordre if c in conds] + sorted(c for c in conds if c not in ordre)
out = {"lag_ns": LAG, "etiquettes": ETIQ, "rejets_lecture": rejets_lecture, "cellules": {}, "differences": {}, "plancher": {}}


def par_systeme(cond, strate):
    """valeurs par systeme = mediane des repetitions."""
    rows = collections.defaultdict(list)
    for L in pos:
        if L["strate"] == strate and L["cond"] == cond:
            rows[L["id"]].append(L)
    res = {}
    for sid, rs in rows.items():
        med = lambda k: float(np.median([r[k] for r in rs]))
        f, ir, lr = med("fnat"), med("irmsd_A"), med("lrmsd_A")
        from_cl = "high" if f >= 0.5 and (lr <= 1 or ir <= 1) else "medium" if f >= 0.3 and (lr <= 5 or ir <= 2) else "acceptable" if f >= 0.1 and (lr <= 10 or ir <= 4) else "incorrect"
        res[sid] = {"fnat": f, "irmsd": ir, "lrmsd": lr, "dockq": med("dockq"), "d_com": med("d_com_A"), "n_contacts": med("n_contacts"), "classe": from_cl, "n_rep": len(rs),
                    "clash": float(np.mean([r["dmin_sortie_A"] < 2.0 for r in rs])), "derive_A": med("derive_A_ca_A"), "derive_B": med("derive_B_ca_A"),
                    "classe_images": [r["classe_capri"] for r in rs], "x0_fnat_contacts": float(np.median([r["n_contacts"] for r in rs]))}
    return res


md = [f"# Tableaux du redocking BioEmu2 — lag {LAG} ns\n", f"Etiquettes lues : {ETIQ}. Systemes lus : {len(systemes)} ; rejetes a la lecture : {len(rejets_lecture)} {[r['id'] + ':' + r['motif'][:40] for r in rejets_lecture]}.\n"]
md.append("Unite = le systeme (mediane des K repetitions), IC95 = bootstrap sur les systemes. fnat toujours avec iRMSD. Classe = criteres CAPRI d origine sur les medianes par systeme.\n")
md.append("## 1. Strate contaminee contre strate tenue a l ecart, par distance et par bras\n")
md.append("| condition | " + " | ".join(strates) + " |\n|---|" + "---|" * len(strates))
for c in conds:
    cells = []
    for s in strates:
        v = par_systeme(c, s)
        if not v:
            cells.append("—"); continue
        f = [x["fnat"] for x in v.values()]; ir = [x["irmsd"] for x in v.values()]; dq = [x["dockq"] for x in v.values()]
        cl = collections.Counter(x["classe"] for x in v.values())
        ok = np.mean([x["classe"] != "incorrect" for x in v.values()])
        cf, ci_ = boot_med(f), boot_med(ir)
        pa = [(sum(k != "incorrect" for k in x["classe_images"]), len(x["classe_images"])) for x in v.values()]
        pm = [(sum(k in ("medium", "high") for k in x["classe_images"]), len(x["classe_images"])) for x in v.values()]
        ra, rm = boot_ratio(pa), boot_ratio(pm)
        out["cellules"][f"{s}|{c}"] = {"images_acceptable_ou_mieux": {"k": int(sum(p[0] for p in pa)), "n": int(sum(p[1] for p in pa)), "taux": ra[0], "ic95": ra[1:]},
                                         "images_medium_ou_mieux": {"k": int(sum(p[0] for p in pm)), "n": int(sum(p[1] for p in pm)), "taux": rm[0], "ic95": rm[1:]}, "n": len(v), "fnat_med": float(np.median(f)), "fnat_ic95": cf, "irmsd_med": float(np.median(ir)), "irmsd_ic95": ci_, "dockq_med": float(np.median(dq)),
                                         "lrmsd_med": float(np.median([x["lrmsd"] for x in v.values()])), "d_com_med": float(np.median([x["d_com"] for x in v.values()])),
                                         "frac_acceptable_ou_mieux": float(ok), "frac_ic95": boot_frac([x["classe"] != "incorrect" for x in v.values()]), "classes": {k: cl.get(k, 0) for k in CL},
                                         "frac_images_clash_lt_2A": float(np.mean([x["clash"] for x in v.values()])), "par_systeme": {k: {kk: vv for kk, vv in x.items() if kk != "classe_images"} for k, x in v.items()}}
        cells.append(f"n={len(v)} · **images ≥acceptable {sum(p[0] for p in pa)}/{sum(p[1] for p in pa)} = {ra[0]:.2f} [{ra[1]:.2f};{ra[2]:.2f}] · ≥medium {sum(p[0] for p in pm)}/{sum(p[1] for p in pm)} = {rm[0]:.2f} [{rm[1]:.2f};{rm[2]:.2f}]** · "
                     f"fnat med {np.median(f):.3f} [{cf[0]:.3f};{cf[1]:.3f}] · iRMSD med {np.median(ir):.2f} [{ci_[0]:.2f};{ci_[1]:.2f}] A · systemes ≥acceptable {int(round(ok * len(v)))}/{len(v)} (" + "/".join(f"{k[:3]}{cl.get(k, 0)}" for k in CL) + ")")
    md.append(f"| {c} | " + " | ".join(cells) + " |")
md.append("\n(inc/acc/med/hig = nombre de systemes par classe.)\n")

md.append("## 2. Differences entre strates (taux de succes par image, IC95 bootstrap sur les systemes ; mediane des fnat et iRMSD par systeme en complement)\n")
paires = [("A_apparie", "B"), ("A_corpus", "B")]
md.append("| condition | " + " | ".join(f"{a} − {b}" for a, b in paires) + " |\n|---|" + "---|" * len(paires))
for c in conds:
    cells = []
    for a_, b_ in paires:
        va, vb = par_systeme(c, a_), par_systeme(c, b_)
        if not va or not vb:
            cells.append("—"); continue
        d = boot_diff_med([x["fnat"] for x in va.values()], [x["fnat"] for x in vb.values()])
        di = boot_diff_med([x["irmsd"] for x in va.values()], [x["irmsd"] for x in vb.values()])
        pa_ = [(sum(k != "incorrect" for k in x["classe_images"]), len(x["classe_images"])) for x in va.values()]
        pb_ = [(sum(k != "incorrect" for k in x["classe_images"]), len(x["classe_images"])) for x in vb.values()]
        dr = boot_diff_ratio(pa_, pb_)
        out["differences"][f"{a_}-{b_}|{c}"] = {"fnat": d, "irmsd": di, "taux_images_acceptable": dr, "n": [len(va), len(vb)]}
        cells.append(f"**taux d images ≥acceptable {dr[0]:+.3f} [{dr[1]:+.3f};{dr[2]:+.3f}]** · fnat med {d[0]:+.3f} [{d[1]:+.3f};{d[2]:+.3f}] · iRMSD med {di[0]:+.2f} [{di[1]:+.2f};{di[2]:+.2f}] A (n {len(va)}/{len(vb)})")
    md.append(f"| {c} | " + " | ".join(cells) + " |")

md.append("\n## 3. Plancher du hasard, plafond et excedent\n")
md.append("Plancher = F2 (placement aleatoire amene au contact) du MEME systeme : p99 et maximum de son fnat. Plafond = condition C0 (le cristal propage).\n")
md.append("| strate | n | fnat F2 p99 (mediane des systemes) | fnat F2 max (mediane) | fnat C0 plafond (mediane) | R_d20 : systemes dont le fnat depasse le MAX de leur propre plancher |\n|---|---|---|---|---|---|")
for s in strates:
    ids = sorted({L["id"] for L in pos if L["strate"] == s})
    p99 = [np.percentile(plan[i]["plancher"]["F2"]["fnat"], 99) for i in ids if i in plan and plan[i]["plancher"]["F2"]]
    mx = [max(plan[i]["plancher"]["F2"]["fnat"]) for i in ids if i in plan and plan[i]["plancher"]["F2"]]
    c0 = par_systeme("C0_d0", s)
    r20 = par_systeme("R_d20", s)
    dep = sum(1 for i, x in r20.items() if i in plan and plan[i]["plancher"]["F2"] and x["fnat"] > max(plan[i]["plancher"]["F2"]["fnat"]))
    out["plancher"][s] = {"n_avec_plancher": len(p99), "fnat_f2_p99_mediane": float(np.median(p99)) if p99 else None, "fnat_f2_max_mediane": float(np.median(mx)) if mx else None,
                          "fnat_c0_mediane": float(np.median([x["fnat"] for x in c0.values()])) if c0 else None, "r20_depasse_son_max": dep, "r20_n": len(r20)}
    md.append(f"| {s} | {len(ids)} | {np.median(p99):.3f} | {np.median(mx):.3f} | {np.median([x['fnat'] for x in c0.values()]) if c0 else float('nan'):.3f} | {dep}/{len(r20)} |")

# gradient par periode de depot (B)
md.append("\n## 4. Strate B par periode de depot (R_d20 et T_d20)\n")
md.append("Taux d images >= acceptable (images / images, systemes) et mediane des fnat par systeme. Si la vraie coupure d entrainement est PLUS TARDIVE que 2021-09-30, la premiere periode serait contaminee.\n")
md.append("| periode de depot | systemes | R_d20 : images >= acceptable | R_d20 fnat med | T_d20 : images >= acceptable | T_d20 fnat med |\n|---|---|---|---|---|---|")
per = collections.defaultdict(list)
for sid in {L["id"] for L in pos if L["strate"] == "B"}:
    per[man_B[sid]["periode"]].append(sid)
vR, vT = par_systeme("R_d20", "B"), par_systeme("T_d20", "B")
for p in sorted(per):
    cells = []
    for v, c in ((vR, "R_d20"), (vT, "T_d20")):
        ids_ = [i for i in per[p] if i in v]
        k_ = sum(sum(x != "incorrect" for x in v[i]["classe_images"]) for i in ids_); n_ = sum(len(v[i]["classe_images"]) for i in ids_)
        out.setdefault("periodes_B", {})[f"{p}|{c}"] = {"systemes": len(ids_), "k": int(k_), "images": int(n_)}
        cells.append(f"{k_}/{n_} ({len(ids_)} syst.) | {np.median([v[i]['fnat'] for i in ids_]):.3f}" if ids_ else "— | —")
    md.append(f"| {p} | {len(per[p])} | " + " | ".join(cells) + " |")
if hom:
    md.append("\n## 5. Strate B selon qu une paire homologue a ete deposee avant la coupure\n")
    md.append("| seuil d identite | avec paire homologue anterieure : n, fnat med (R_d20) | sans : n, fnat med (R_d20) |\n|---|---|---|")
    for k in ("id90", "id50", "id30"):
        oui = [i for i in vR if i in hom and hom[i][k]["entrees_avec_les_deux"] > 0]; non = [i for i in vR if i in hom and hom[i][k]["entrees_avec_les_deux"] == 0]
        f = lambda L_: f"n={len(L_)}, {np.median([vR[i]['fnat'] for i in L_]):.3f}" if L_ else "n=0"
        md.append(f"| {k[2:]} % | {f(oui)} | {f(non)} |")

# gradient d homologie (niveau 0 = aucune paire homologue anterieure >= 30 % ; 1 = 30-50 % ; 2 = 50-90 % ; 3 = >= 90 % ; 4 = deposes AVANT la coupure : A_apparie, A_corpus)
def niveau_hom(pid):
    h = hom[pid]
    return 3 if h["id90"]["entrees_avec_les_deux"] > 0 else 2 if h["id50"]["entrees_avec_les_deux"] > 0 else 1 if h["id30"]["entrees_avec_les_deux"] > 0 else 0


def rho_perm(x, y, n=10000):
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 4 or np.std(x) == 0 or np.std(y) == 0:
        return [None, None]
    r = spearmanr(x, y)[0]
    perm = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n // 10)])
    return [float(r), float((np.abs(perm) >= abs(r)).mean())]


if hom:
    md.append("\\n## 5c. GRADIENT : taux de succes par image selon la ressemblance a une paire deja vue (niveau d homologie)\\n")
    md.append("Niveau 0 = aucune entree anterieure a la coupure ne porte un homologue (>= 30 %) des DEUX chaines ; 1 = paire homologue a 30-50 % ; 2 = 50-90 % ; 3 = >= 90 % ; 4 = systemes deposes AVANT la coupure (A_apparie, A_corpus). "
              "Cellule = images >= acceptable / images, [IC95 bootstrap sur les systemes]. rho = Spearman entre niveau (0-3, B seul) et taux par systeme, p par permutation.\\n")
    md.append("| condition | niveau 0 | niveau 1 | niveau 2 | niveau 3 | niveau 4 (avant coupure) | rho (B, niveaux 0-3) |\\n|---|---|---|---|---|---|---|")
    for c in [x for x in conds if x in ("T_d20", "R_d20", "R_d30", "RD_d20")]:
        cells, xs, ys = [], [], []
        for nv in (0, 1, 2, 3, 4):
            if nv < 4:
                ids = [i for i in {L["id"] for L in pos if L["strate"] == "B" and L["cond"] == c} if niveau_hom(i) == nv]
            else:
                ids = sorted({L["id"] for L in pos if L["strate"] in ("A_apparie", "A_corpus") and L["cond"] == c})
            pr = []
            for i in ids:
                imgs = [L["classe_capri"] for L in pos if L["id"] == i and L["cond"] == c]
                pr.append((sum(k != "incorrect" for k in imgs), len(imgs)))
                if nv < 4:
                    xs.append(nv); ys.append(pr[-1][0] / pr[-1][1])
            if pr:
                rr_ = boot_ratio(pr)
                out.setdefault("gradient", {})[f"{c}|niveau{nv}"] = {"n_systemes": len(pr), "k": int(sum(p[0] for p in pr)), "images": int(sum(p[1] for p in pr)), "taux": rr_[0], "ic95": rr_[1:]}
                cells.append(f"n={len(pr)} · {int(sum(p[0] for p in pr))}/{int(sum(p[1] for p in pr))} = {rr_[0]:.2f} [{rr_[1]:.2f};{rr_[2]:.2f}]")
            else:
                cells.append("—")
        rp = rho_perm(xs, ys)
        out.setdefault("gradient", {})[f"{c}|rho"] = rp
        md.append(f"| {c} | " + " | ".join(cells) + f" | {('rho=%+.2f, p=%.3f (n=%d)' % (rp[0], rp[1], len(xs))) if rp[0] is not None else '—'} |")

# profondeur du MSA : confusion possible (un systeme sans homologue est AUSSI un systeme a MSA pauvre)
def profondeur(pid):
    f = R / f"msa/{pid}.json"
    if not f.is_file():
        return None
    d = json.loads(f.read_text())
    return min(c["non_apparie"] for c in d["chaines"].values()) if d.get("chaines") else None


md.append("\\n## 5d. Taux de succes par image selon la profondeur du MSA (le plus pauvre des deux alignements non apparies, requete comprise)\\n")
md.append("Un systeme sans homologue connu a aussi un MSA pauvre : sans cette table, « inedit » et « MSA pauvre » se confondent. Cellule = images >= acceptable / images (systemes).\\n")
md.append("| condition | groupe | MSA < 100 | 100 a 999 | >= 1000 |\\n|---|---|---|---|---|")
for c in [x for x in conds if x in ("T_d20", "R_d20", "RD_d20")]:
    for grp, noms in (("avant la coupure (A)", ("A_apparie", "A_corpus")), ("apres la coupure (B)", ("B",))):
        cells = []
        for lo, hi in ((0, 100), (100, 1000), (1000, 10 ** 9)):
            ids = sorted({L["id"] for L in pos if L["strate"] in noms and L["cond"] == c and profondeur(L["id"]) is not None and lo <= profondeur(L["id"]) < hi})
            k = n_ = 0
            for i in ids:
                imgs = [L["classe_capri"] for L in pos if L["id"] == i and L["cond"] == c]
                k += sum(x != "incorrect" for x in imgs); n_ += len(imgs)
            out.setdefault("profondeur_msa", {})[f"{c}|{grp}|{lo}-{hi}"] = {"systemes": len(ids), "k": int(k), "images": int(n_)}
            cells.append(f"{k}/{n_} ({len(ids)} syst.)" if n_ else "—")
        md.append(f"| {c} | {grp} | " + " | ".join(cells) + " |")

# sensibilite : B sans les systemes dont le CRISTAL propage (C0) derive deja de plus de 3 A, B a MSA >= 100, B sans les deux defauts ; comparaison A_apparie
md.append("\\n## 5e. Sensibilite : le contraste A_apparie / B tient-il quand on retire les systemes les plus durs de B ?\\n")
md.append("Retire de B : (i) les systemes dont une chaine derive de plus de 3 A MEME en partant du cristal (C0, mediane des images) ; (ii) ceux dont le MSA le plus pauvre a moins de 100 sequences ; (iii) les deux. Cellule = images >= acceptable / images (systemes) [IC95 bootstrap sur les systemes].\\n")
c0max = {}
for L in pos:
    if L["cond"] == "C0_d0":
        c0max.setdefault(L["id"], []).append(max(L["derive_A_ca_A"], L["derive_B_ca_A"]))
c0med = {i: float(np.median(v)) for i, v in c0max.items()}
ids_B = sorted({L["id"] for L in pos if L["strate"] == "B"})
ens = {"B (tous)": ids_B,
       "B sans derive C0 > 3 A": [i for i in ids_B if c0med.get(i, 0) <= 3.0],
       "B avec MSA >= 100": [i for i in ids_B if (profondeur(i) or 0) >= 100],
       "B sans derive ET MSA >= 100": [i for i in ids_B if c0med.get(i, 0) <= 3.0 and (profondeur(i) or 0) >= 100],
       "A_apparie sans derive ET MSA >= 100 (meme retrait que B)": [i for i in sorted({L["id"] for L in pos if L["strate"] == "A_apparie"}) if c0med.get(i, 0) <= 3.0 and (profondeur(i) or 0) >= 100],
       "A_apparie (reference)": sorted({L["id"] for L in pos if L["strate"] == "A_apparie"}),
       "A_corpus (reference)": sorted({L["id"] for L in pos if L["strate"] == "A_corpus"})}
md.append("| ensemble | systemes | T_d20 | R_d20 | R_d30 | RD_d20 |\\n|---|---|---|---|---|---|")
for nom_, ids_ in ens.items():
    cells = []
    for c in ("T_d20", "R_d20", "R_d30", "RD_d20"):
        pr = []
        for i in ids_:
            imgs = [L["classe_capri"] for L in pos if L["id"] == i and L["cond"] == c]
            if imgs:
                pr.append((sum(k != "incorrect" for k in imgs), len(imgs)))
        if pr:
            r_ = boot_ratio(pr)
            out.setdefault("sensibilite", {})[f"{nom_}|{c}"] = {"systemes": len(pr), "k": int(sum(p[0] for p in pr)), "images": int(sum(p[1] for p in pr)), "taux": r_[0], "ic95": r_[1:]}
            cells.append(f"{int(sum(p[0] for p in pr))}/{int(sum(p[1] for p in pr))} = {r_[0]:.2f} [{r_[1]:.2f};{r_[2]:.2f}] ({len(pr)} syst.)")
        else:
            cells.append("—")
    md.append(f"| {nom_} | {len(ids_)} | " + " | ".join(cells) + " |")

# differences A_apparie − B avec IC95, sur ensembles SYMETRIQUES (meme retrait des deux cotes)
md.append("\\n### 5e-bis. Difference A_apparie − B (taux d images >= acceptable) [IC95 bootstrap sur les systemes], ensembles symetriques\\n")
md.append("| comparaison | T_d10 | T_d20 | T_d30 | R_d10 | R_d20 | R_d30 | RD_d20 |\\n|---|---|---|---|---|---|---|---|")


def paires_img(ids_, c):
    pr = []
    for i in ids_:
        imgs = [L["classe_capri"] for L in pos if L["id"] == i and L["cond"] == c]
        if imgs:
            pr.append((sum(k != "incorrect" for k in imgs), len(imgs)))
    return pr


for etiq_, (ka, kb) in (("tous contre tous", ("A_apparie (reference)", "B (tous)")),
                        ("sans derive ET MSA >= 100, des deux cotes", ("A_apparie sans derive ET MSA >= 100 (meme retrait que B)", "B sans derive ET MSA >= 100"))):
    cells = []
    for c in ("T_d10", "T_d20", "T_d30", "R_d10", "R_d20", "R_d30", "RD_d20"):
        pa_, pb_ = paires_img(ens[ka], c), paires_img(ens[kb], c)
        dr_ = boot_diff_ratio(pa_, pb_)
        out.setdefault("diff_symetrique", {})[f"{etiq_}|{c}"] = {"n_A": len(pa_), "n_B": len(pb_), "diff": dr_}
        cells.append(f"{dr_[0]:+.2f} [{dr_[1]:+.2f};{dr_[2]:+.2f}] ({len(pa_)}/{len(pb_)})" if dr_[0] is not None else "—")
    md.append(f"| {etiq_} | " + " | ".join(cells) + " |")

# deformation des chaines (repere : 3 A, celui du garde recepteur G8 du pilote 6CWG ; a 1000 ns le pilote mesurait 5,30 A de derive du recepteur)
md.append("\n## 5b. Deformation des chaines : RMSD Calpha de chaque chaine a sa forme cristalline (superposition sur elle-meme)\n")
md.append("Un succes CAPRI obtenu en DEFORMANT les chaines serait trompeur. Repere 3 A (garde recepteur G8 du pilote 6CWG, qui mesurait 5,30 A a 1000 ns). Mediane sur les images, puis sur les systemes.\n")
md.append("| strate | condition | derive chaine A, mediane (A) | derive chaine B, mediane (A) | images >= acceptable avec une chaine > 3 A | toutes images avec une chaine > 3 A |\n|---|---|---|---|---|---|")
for s_ in strates:
    for c in conds:
        rs = [L for L in pos if L["strate"] == s_ and L["cond"] == c]
        if not rs:
            continue
        dA = float(np.median([L["derive_A_ca_A"] for L in rs])); dB = float(np.median([L["derive_B_ca_A"] for L in rs]))
        ok_ = [L for L in rs if L["classe_capri"] != "incorrect"]
        gr = lambda LL: sum(1 for L in LL if max(L["derive_A_ca_A"], L["derive_B_ca_A"]) > 3.0)
        out.setdefault("deformation", {})[f"{s_}|{c}"] = {"derive_A_med": dA, "derive_B_med": dB, "succes_avec_deformation_gt3A": gr(ok_), "succes": len(ok_), "images_avec_deformation_gt3A": gr(rs), "images": len(rs)}
        md.append(f"| {s_} | {c} | {dA:.2f} | {dB:.2f} | {gr(ok_)}/{len(ok_)} | {gr(rs)}/{len(rs)} |")

# temoin negatif
neg = [L for L in lignes if L["kind"] == "neg"]
if neg:
    md.append("\n## 6. Temoin negatif (paires sans bonne reponse)\n")
    md.append("| condition | n paires | d_COM entree (A) | d_COM sortie (A) | contacts residu-residu sortie (mediane) | distance atomique min sortie (A) | images avec >= 25 contacts | rappel du site de la chaine A | precision du site | base au hasard |\n|---|---|---|---|---|---|---|---|---|---|")
    for c in ["R_d10", "R_d20", "R_d30", "RD_d20"]:
        rows = collections.defaultdict(list)
        for L in neg:
            if L["cond"] == c:
                rows[L["id"]].append(L)
        if not rows:
            continue
        P = {sid: {k: float(np.median([r[k] for r in rs if r.get(k) is not None])) for k in ("x0_d_com_A", "d_com_A", "n_contacts", "dmin_sortie_A", "rappel_site", "precision_site") if any(r.get(k) is not None for r in rs)} for sid, rs in rows.items()}
        frac25 = float(np.mean([r["n_contacts"] >= 25 for rs in rows.values() for r in rs]))
        base = float(np.median([systemes[sid]["mes"]["base_hasard_precision"] for sid in rows]))
        g = lambda k: [p[k] for p in P.values() if k in p]
        out.setdefault("negatif", {})[c] = {"n": len(P), "d_com_entree": float(np.median(g("x0_d_com_A"))), "d_com_sortie": float(np.median(g("d_com_A"))), "contacts_med": float(np.median(g("n_contacts"))),
                                            "contacts_ic95": boot_med(g("n_contacts")), "dmin_med": float(np.median(g("dmin_sortie_A"))), "frac_images_ge25_contacts": frac25,
                                            "rappel_site_med": float(np.median(g("rappel_site"))) if g("rappel_site") else None, "precision_site_med": float(np.median(g("precision_site"))) if g("precision_site") else None, "base": base,
                                            "par_paire": P}
        md.append(f"| {c} | {len(P)} | {np.median(g('x0_d_com_A')):.1f} | {np.median(g('d_com_A')):.1f} | {np.median(g('n_contacts')):.0f} {boot_med(g('n_contacts'))} | {np.median(g('dmin_sortie_A')):.2f} | {frac25:.2f} | "
                  f"{np.median(g('rappel_site')):.2f} | {np.median(g('precision_site')):.2f} | {base:.2f} |")
    md.append("\n### 6b. Temoin negatif par groupe : chaines A d AVANT la coupure (NEG_A) contre POSTERIEURES a la coupure (NEG_B)\n")
    md.append("| condition | groupe | paires | images >= 25 contacts | contacts (mediane) | d_COM de sortie (A) | rappel du site de A | precision du site | base au hasard |\n|---|---|---|---|---|---|---|---|---|")
    for c in ["R_d10", "R_d20", "R_d30", "RD_d20"]:
        for g in ("NEG_A", "NEG_B"):
            rows = collections.defaultdict(list)
            for L in neg:
                if L["cond"] == c and L["strate"] == g:
                    rows[L["id"]].append(L)
            if not rows:
                continue
            P = {sid: {k: float(np.median([r[k] for r in rs if r.get(k) is not None])) for k in ("d_com_A", "n_contacts", "rappel_site", "precision_site") if any(r.get(k) is not None for r in rs)} for sid, rs in rows.items()}
            fr = float(np.mean([r["n_contacts"] >= 25 for rs in rows.values() for r in rs]))
            gg = lambda k: [p[k] for p in P.values() if k in p]
            base_ = float(np.median([systemes[sid]["mes"]["base_hasard_precision"] for sid in rows]))
            out.setdefault("negatif_groupes", {})[f"{c}|{g}"] = {"paires": len(P), "frac_images_ge25": fr, "contacts_med": float(np.median(gg("n_contacts"))), "d_com_sortie": float(np.median(gg("d_com_A"))),
                                                                  "rappel_site": float(np.median(gg("rappel_site"))), "precision_site": float(np.median(gg("precision_site"))), "base": base_}
            md.append(f"| {c} | {g} | {len(P)} | {fr:.2f} | {np.median(gg('n_contacts')):.0f} | {np.median(gg('d_com_A')):.1f} | {np.median(gg('rappel_site')):.2f} | {np.median(gg('precision_site')):.2f} | {base_:.2f} |")
    # comparaison avec les contacts NATIFS des positifs
    nat = [S["meta"]["n_contacts_natifs"] for S in systemes.values() if S["meta"]["kind"] == "pos"]
    if nat:
        md.append(f"\nRepere : contacts natifs des complexes positifs, mediane {np.median(nat):.0f} couples de residus (min {min(nat)}, max {max(nat)}).")
# annexe : un systeme par ligne
man_Aa = json.loads((R / "strate_A_apparie/manifeste.json").read_text())["systemes"] if (R / "strate_A_apparie/manifeste.json").is_file() else {}
man_Ac = json.loads((R / "strate_A_corpus/manifeste.json").read_text()) if (R / "strate_A_corpus/manifeste.json").is_file() else {}


def depot_de(sid, strate):
    if strate == "B":
        return man_B[sid]["depot"]
    if strate == "A_apparie":
        return man_Aa[sid]["depot"]
    if strate == "A_corpus":
        return man_Ac[sid]["date"]["deposit_date"][:10]
    return ""


md.append("\n## 7. Annexe : un systeme par ligne (positifs)\n")
md.append("Succes = images >= acceptable / images. Derive C0 = RMSD Calpha maximal des deux chaines (mediane des images C0). Niveau = niveau d homologie (B seul). MSA = profondeur du plus pauvre des deux alignements.\n")
md.append("| strate | id | depot | niveau | residus | MSA | d_COM natif (A) | contacts natifs | derive C0 (A) | T_d20 | R_d20 | RD_d20 |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
for s_ in strates:
    for sid in sorted({L["id"] for L in pos if L["strate"] == s_}):
        rs_ = [L for L in pos if L["id"] == sid]
        cs = lambda c: (lambda im: f"{sum(x != 'incorrect' for x in im)}/{len(im)}" if im else "—")([L["classe_capri"] for L in rs_ if L["cond"] == c])
        c0 = [max(L["derive_A_ca_A"], L["derive_B_ca_A"]) for L in rs_ if L["cond"] == "C0_d0"]
        niv_ = niveau_hom(sid) if (hom and s_ == "B" and sid in hom) else "—"
        md.append(f"| {s_} | {sid} | {depot_de(sid, s_) if s_ != 'special' else ''} | {niv_} | {rs_[0]['tokens']} | {profondeur(sid)} | {rs_[0]['d_nat_A']:.1f} | {rs_[0]['n_contacts_natifs']} | {np.median(c0):.2f} | {cs('T_d20')} | {cs('R_d20')} | {cs('RD_d20')} |" if c0 else f"| {s_} | {sid} | | | {rs_[0]['tokens']} | {profondeur(sid)} | | | — | {cs('T_d20')} | {cs('R_d20')} | {cs('RD_d20')} |")
(R / f"resultats/tableaux_lag{LAG:04d}.md").write_text("\n".join(md) + "\n")
(R / f"resultats/tableaux_lag{LAG:04d}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str) + "\n")
print("\n".join(md))
print(f"BILAN_ANALYSE : {len(systemes)} systemes lus, {len(rejets_lecture)} rejetes a la lecture (REJETS = {len(rejets_lecture)}), {len(lignes)} images")
