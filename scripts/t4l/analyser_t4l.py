#!/usr/bin/env python3
"""analyser_t4l.py — populations des quatre etats (SES CV, SES seuils), distribution de d1, plancher, couverture, comparaison.

Lit : <racine>/<unite>/cv.npz (d1, d2 en nm ; p ; physique), les references de l'auteur (Zenodo), 5LZM.
Ecrit : resultats/analyse_<racine>.json + figs/fig_comparaison_<racine>.png
Les populations sont donnees sur (a) TOUTES les structures, (b) les structures qui passent la physicalite ; les rejetees sont comptees, jamais retirees en silence.
"""
import json, os, sys
import numpy as np, mdtraj as md
from scipy import stats
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = "/gpfs/phherman/bioemu2_bench/t4l"
Z = f"{T}/donnees/zenodo/extrait/T4L-bechmarking"
RACINE = sys.argv[1] if len(sys.argv) > 1 else "unites"
ROOT = f"{T}/{RACINE}"
SES_CHIFFRES = {  # lus dans le papier (Fig. 2C, 2D), pas recalcules : l'EMD (27,6 us) n'est pas publie
    "EMD": {"exposed/open": 39.1, "exposed/closed": 3.9, "buried/open": 54.0, "buried/closed": 2.9},
    "AMS": {"exposed/open": 38.3, "exposed/closed": 0.7, "buried/open": 58.7, "buried/closed": 2.3}}


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    ph = k / n; den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den; h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return [round(100 * max(0, c - h), 2), round(100 * min(1, c + h), 2)]


def pops(d1, p):
    e = C.etat(d1, p); n = len(e)
    out = {"n": int(n), "n_indefini": int((e == -1).sum()), "etats": {}}
    for k, nom in enumerate(C.NOMS):
        c = int((e == k).sum()); out["etats"][nom] = {"n": c, "pct": round(100 * c / max(1, n - out["n_indefini"]), 2), "IC95_wilson_pct": wilson(c, n - out["n_indefini"])}
    return out


def d1_distribution(d1):
    q = np.percentile(d1, [0, 5, 25, 50, 75, 95, 100])
    return {"quantiles_0_5_25_50_75_95_100_nm": [round(float(x), 3) for x in q],
            "ferme_lt_2_5": int((d1 < 2.5).sum()), "ferme_robuste_lt_2_4": int((d1 < 2.4).sum()),
            "ambigu_2_4_a_2_6": int(((d1 >= 2.4) & (d1 < 2.6)).sum()), "distance_mediane_au_seuil_nm": round(float(np.median(d1) - 2.5), 3),
            "plus_proche_du_seuil_nm": round(float(np.min(np.abs(d1 - 2.5))), 3),
            "histogramme_0_25_nm": {f"{a:.2f}-{a+0.25:.2f}": int(((d1 >= a) & (d1 < a + 0.25)).sum()) for a in np.arange(1.5, 5.5, 0.25)}}


def couverture(d1, d2, p):
    """cases occupees sur des grilles de 0,25 nm (d1,d2) et de 0,25 (p) ; entropie normalisee de l'histogramme 1D de d1 et de p."""
    g = lambda x, lo, hi, w: np.floor((np.clip(x, lo, hi - 1e-9) - lo) / w).astype(int)
    cases = set(zip(g(d1, 1.0, 6.0, 0.25).tolist(), g(d2, 0.0, 6.0, 0.25).tolist()))
    def ent(x, lo, hi, w):
        h = np.bincount(g(x, lo, hi, w), minlength=int(round((hi - lo) / w))).astype(float); h = h[h > 0] / h.sum()
        return round(float(-(h * np.log(h)).sum() / np.log(int(round((hi - lo) / w)))), 3)
    pf = p[np.isfinite(p)]
    return {"cases_d1_d2_occupees_(0,25nm)": len(cases), "entropie_norm_d1": ent(d1, 1.0, 6.0, 0.25), "entropie_norm_p": ent(pf, -3.0, 3.0, 0.25) if pf.size else None,
            "d1_min_max": [round(float(d1.min()), 3), round(float(d1.max()), 3)], "d2_min_max": [round(float(d2.min()), 3), round(float(d2.max()), 3)],
            "p_min_med_max": [round(float(pf.min()), 3), round(float(np.median(pf)), 3), round(float(pf.max()), 3)] if pf.size else None,
            "p_proche_de_0_abs_lt_0_05": int((np.abs(pf) < 0.05).sum()) if pf.size else None}


def bruit(d1, sd=0.1, B=4000, seed=0):
    """Plancher : la population « fermee » que produirait un bruit gaussien de +/- sd nm sur d1, applique a CE MEME ensemble (compare a l'observe).
    Si le ferme observe n'est pas au-dessus de [2,5 ; 97,5] % de cette distribution, il ne depasse pas le bruit de mesure."""
    rng = np.random.default_rng(seed)
    f = np.array([(d1 + rng.normal(0, sd, d1.size) < 2.5).mean() for _ in range(B)]) * 100
    return {"sd_nm": sd, "ferme_observe_pct": round(100 * float((d1 < 2.5).mean()), 3), "ferme_sous_bruit_pct_2_5_50_97_5": [round(float(x), 3) for x in np.percentile(f, [2.5, 50, 97.5])]}


# ── references de l'auteur, recalculees avec MON code ───────────────────────────
def charger_auteur():
    out = {}
    for nom, (d, x, t) in {"AlphaFlow": ("AlphaFlow", "combined.xtc", "prot-masses.pdb"), "ESMFlow": ("ESMFlow", "combined.xtc", "prot-masses.pdb"),
                           "AF-cluster": ("AF-cluster", "combined.xtc", "prot-masses.pdb"), "ConforFold": ("ConforFold", "combined.xtc", "prot-masses.pdb"),
                           "ConfRover": ("ConfRover", "T4L.xtc", "T4L.pdb"), "rMSA-AF2": ("rMSA-AF2", "combined.xtc", "prot-masses.pdb"),
                           "rMSA-AF2-recyc": ("rMSA-AF2-recyc", "combined.xtc", "prot-masses.pdb"), "BioEmu_v1": ("BioEmu", "samples.xtc", "topology.pdb")}.items():
        r = C.cv(md.load(f"{Z}/Genai-models/{d}/{x}", top=f"{Z}/Genai-models/{d}/{t}")); out[nom] = r
    return out


def main():
    unites = sorted(d for d in os.listdir(ROOT) if os.path.isfile(f"{ROOT}/{d}/cv.npz"))
    res, rejets = {"racine": RACINE, "unites": {}}, []
    U = {}
    for u in unites:
        z = np.load(f"{ROOT}/{u}/cv.npz"); m = json.load(open(f"{ROOT}/{u}/meta.json"))
        d1, d2, p, ph = z["d1"], z["d2"], z["p"], z["physique"].astype(bool)
        U[u] = (d1, d2, p, ph)
        e = {"type": m["type"], "systeme": m["systeme"], "lag_ns": m.get("lag_ns"), "n": int(len(d1)), "physiques": int(ph.sum()), "rejetees_physicalite": int((~ph).sum()),
             "s_par_structure": m["s_par_structure"], "slurm_job": m.get("slurm_job"), "sampler": m["sampler"],
             "toutes": pops(d1, p), "physiques_seules": pops(d1[ph], p[ph]) if ph.any() else None,
             "d1": d1_distribution(d1), "couverture": couverture(d1, d2, p), "bruit_d1": bruit(d1)}
        res["unites"][u] = e
    A = charger_auteur()
    pool = ["AlphaFlow", "ESMFlow", "AF-cluster", "ConforFold", "ConfRover", "rMSA-AF2", "rMSA-AF2-recyc"]
    d1g = np.concatenate([A[k]["d1"] for k in pool]); pg = np.concatenate([A[k]["p"] for k in pool]); d2g = np.concatenate([A[k]["d2"] for k in pool])
    res["auteur"] = {"GenAI_sept_methodes_vrai_p": {"pops": pops(d1g, pg), "d1": d1_distribution(d1g), "couverture": couverture(d1g, d2g, pg)},
                     "BioEmu_v1": {"d1": d1_distribution(A["BioEmu_v1"]["d1"]), "couverture_d1_d2": couverture(A["BioEmu_v1"]["d1"], A["BioEmu_v1"]["d2"], np.array([np.nan])),
                                   "p_ring": "NON DEFINI (cycle de Phe4 absent)", "pops_proxy_CB": pops(A["BioEmu_v1"]["d1"], A["BioEmu_v1"]["p_cb_proxy"])},
                     "EMD_AMS_chiffres_du_papier": SES_CHIFFRES}
    # comparaison de d1 : KS et Wasserstein BioEmu2 (equilibrium s164) contre BioEmu v1 et contre le pool GenAI
    cmp = {}
    if "eq164_n500" in U:
        d = U["eq164_n500"][0]
        for nom, ref in (("BioEmu_v1", A["BioEmu_v1"]["d1"]), ("GenAI_pool", d1g)):
            ks = stats.ks_2samp(d, ref)
            cmp[f"d1_eq164_contre_{nom}"] = {"KS": round(float(ks.statistic), 3), "p": float(ks.pvalue), "wasserstein_nm": round(float(stats.wasserstein_distance(d, ref)), 3)}
    res["comparaison_d1"] = cmp
    # vrais positifs de l'auteur (221), points de reference de la figure
    ref_closed = {}
    import glob
    for lab in ("Exposed-closed", "Buried-closed"):
        a = [C.cv(md.load(f)) for f in sorted(glob.glob(f"{Z}/Closed-state/{lab}/*.pdb"))]
        ref_closed[lab] = (np.array([r["d1"][0] for r in a]), np.array([r["d2"][0] for r in a]), np.array([r["p"][0] for r in a]))
    cris = C.cv(md.load(f"{T}/x0/x0_5lzm.pdb"))
    res["cristal_5LZM"] = {"d1_nm": round(float(cris["d1"][0]), 4), "d2_nm": round(float(cris["d2"][0]), 4), "p": round(float(cris["p"][0]), 4)}
    json.dump(res, open(f"{T}/resultats/analyse_{RACINE}.json", "w"), indent=1, ensure_ascii=False)

    # ── figure ──
    os.makedirs(f"{T}/figs", exist_ok=True)
    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 5)
    couleurs = ["#DD8452", "#4C72B0", "#C44E52", "#55A868"]   # exp/open, exp/closed, bur/open, bur/closed (ses couleurs)
    def camembert(ax, pops_, titre):
        v = [pops_[n] for n in C.NOMS]
        w = [x for x in v if x > 0]; c = [couleurs[i] for i, x in enumerate(v) if x > 0]; lab = [C.NOMS[i] for i, x in enumerate(v) if x > 0]
        ax.pie(w, colors=c, labels=lab, autopct=lambda q: f"{q:.1f}%" if q > 0.05 else "", textprops={"fontsize": 8}); ax.set_title(titre, fontsize=10)
    eq = res["unites"].get("eq164_n500")
    if eq:
        camembert(fig.add_subplot(gs[0, 0]), {n: eq["toutes"]["etats"][n]["pct"] for n in C.NOMS}, f"BioEmu2 equilibrium (n={eq['n']})")
    ax = fig.add_subplot(gs[0, 1]); b = res["auteur"]["BioEmu_v1"]["pops_proxy_CB"]
    camembert(ax, {n: b["etats"][n]["pct"] for n in C.NOMS}, "BioEmu v1 — p = PROXY CB\n(p du papier non defini)")
    camembert(fig.add_subplot(gs[0, 2]), {n: res["auteur"]["GenAI_sept_methodes_vrai_p"]["pops"]["etats"][n]["pct"] for n in C.NOMS}, "GenAI (7 methodes), SON vrai p\n(son camembert = 100 % : lock = d2)")
    camembert(fig.add_subplot(gs[0, 3]), SES_CHIFFRES["AMS"], "AMS (chiffres du papier)")
    camembert(fig.add_subplot(gs[0, 4]), SES_CHIFFRES["EMD"], "EMD 27,6 us (chiffres du papier)")
    ax = fig.add_subplot(gs[1, 0:2])
    bins = np.arange(1.5, 5.6, 0.1)
    ax.hist(A["BioEmu_v1"]["d1"], bins=bins, histtype="step", lw=2, color="grey", label=f"BioEmu v1 (n={len(A['BioEmu_v1']['d1'])})", density=True)
    cyc = ["#C44E52", "#DD8452", "#55A868", "#4C72B0", "#8172B3", "#937860"]
    liste = (("eq164_n500", "#C44E52"), ("eq162_n250", "#DD8452"), ("prop162_lag0001", "#55A868"), ("prop162_lag0050", "#4C72B0"), ("prop162_lag1000", "#8172B3")) if RACINE == "unites" else tuple((u, cyc[i % 6]) for i, u in enumerate(U))
    for u, col in liste:
        if u in U:
            ax.hist(U[u][0], bins=bins, histtype="step", lw=1.6, color=col, label=f"{u} (n={len(U[u][0])})", density=True)
    ax.axvline(2.5, color="k", ls="--"); ax.axvline(cris["d1"][0], color="g", ls=":", label="5LZM cristal")
    ax.set_xlabel("d1 = Ca Ser44 - Ile150 (nm)"); ax.set_ylabel("densite"); ax.legend(fontsize=7); ax.set_title("distribution de d1 — le seuil 2,5 nm est la ligne tiretee")
    ax = fig.add_subplot(gs[1, 2:4])
    for k in ("Exposed-closed", "Buried-closed"):
        ax.scatter(ref_closed[k][0], ref_closed[k][2], s=8, marker="x", label=f"auteur : {k} (n={len(ref_closed[k][0])})")
    if "eq164_n500" in U:
        ax.scatter(U["eq164_n500"][0], U["eq164_n500"][2], s=6, alpha=.5, color="#C44E52", label="BioEmu2 eq164")
    ax.scatter(d1g, pg, s=3, alpha=.25, color="grey", label="GenAI (7), son vrai p")
    ax.scatter(cris["d1"], cris["p"], s=90, marker="*", color="g", label="5LZM")
    ax.axvline(2.5, color="k", ls="--"); ax.axhline(0, color="k", ls="--"); ax.set_xlabel("d1 (nm)"); ax.set_ylabel("p"); ax.set_ylim(-3, 3); ax.legend(fontsize=7); ax.set_title("plan (d1, p)")
    ax = fig.add_subplot(gs[1, 4])
    for k in ("Exposed-closed", "Buried-closed"):
        ax.scatter(ref_closed[k][0], ref_closed[k][1], s=8, marker="x")
    if "eq164_n500" in U:
        ax.scatter(U["eq164_n500"][0], U["eq164_n500"][1], s=6, alpha=.5, color="#C44E52")
    ax.scatter(A["BioEmu_v1"]["d1"], A["BioEmu_v1"]["d2"], s=4, alpha=.3, color="grey")
    ax.axvline(2.5, color="k", ls="--"); ax.set_xlabel("d1 (nm)"); ax.set_ylabel("d2 (nm)"); ax.set_title("plan (d1, d2)")
    fig.tight_layout(); fig.savefig(f"{T}/figs/fig_comparaison_{RACINE}.png", dpi=130)
    print(json.dumps({u: {"n": e["n"], "phys": e["physiques"], "pops_toutes_pct": {k: v["pct"] for k, v in e["toutes"]["etats"].items()}, "ferme": e["d1"]["ferme_lt_2_5"],
                         "ferme_robuste": e["d1"]["ferme_robuste_lt_2_4"], "d1_quantiles": e["d1"]["quantiles_0_5_25_50_75_95_100_nm"]} for u, e in res["unites"].items()}, indent=1))
    print("comparaison d1 :", json.dumps(cmp))
    print(f"REJETS : {len(rejets)} · unites lues : {len(unites)}")


if __name__ == "__main__":
    main()
