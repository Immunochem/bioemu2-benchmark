#!/usr/bin/env python3
"""h_agreger.py — agregation DESCRIPTIVE des tables DynaPIN des piles U3 (aucune interpretation) :
  dynapin/agregat_U3_par_image.tsv  : une ligne par (unite, image >= 1) : DockQ, iRMSD, lRMSD, fnat, fnonnat, clashes, RMSD complexe/chaines, et le NOMBRE de
                                      contacts par type d interaction (une colonne par itype rencontre)
  dynapin/agregat_U3_par_unite.tsv  : une ligne par unite : code DynaPIN, garde d empilement (image 0 : DockQ et fnat = 1), mediane/min/max des metriques, contacts moyens
Une pile sans dynapin.code == 0 ou sans table est une LIGNE `echec` du bilan, jamais un saut. Ecrit aussi le bilan au manifeste (type dynapin_pile)."""
import csv, json, os, statistics, sys
from collections import Counter, defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L

US = os.environ.get("H_USAGE", "U3")
HI = L.G / os.environ.get("H_RACINE", "hydrogenes")
PILES = HI / "piles"
OUT = HI / "dynapin"
OUT.mkdir(exist_ok=True)
JOB = os.environ.get("SLURM_JOB_ID")
METR = ["DockQ", "iRMSD", "lRMSD", "fnat", "fnonnat", "clashes", "Complex RMSD", "Chain A RMSD", "Chain B RMSD"]
unites = sorted(d.name for d in PILES.iterdir() if d.is_dir())
itypes = set()
images, par_unite, echecs = [], [], []
for u in unites:
    d = PILES / u
    pdb, usage, regime = u.split("__", 2)
    code = (d / "dynapin.code").read_text().strip() if (d / "dynapin.code").is_file() else None
    qc, ip = d / "dynapin/tables/qc_trajectory_metrics.csv", d / "dynapin/tables/int_pairwise_trajectory.csv"
    ligne = {"type": "dynapin_pile", "unite": u, "systeme": pdb, "usage": usage, "regime": "dynapin_v2b", "jobid": JOB, "pile": str(d)}
    if code != "0" or not qc.is_file() or not ip.is_file():
        raison = f"dynapin.code = {code!r}, qc present {qc.is_file()}, contacts present {ip.is_file()}" + (" ; pile refusee par les gardes" if (d / "pile_a_ne_pas_juger.txt").is_file() else "")
        ligne.update(statut="echec", erreur=raison)
        echecs.append({"unite": u, "raison": raison})
        L.ajouter_manifeste(ligne)
        continue
    rows = list(csv.DictReader(open(qc)))
    cnt = defaultdict(Counter)
    for r in csv.DictReader(open(ip)):
        cnt[int(r["Frame"])][r["itype"]] += 1
        itypes.add(r["itype"])
    z = rows[0]
    garde = abs(float(z["DockQ"]) - 1.0) < 1e-6 and abs(float(z["fnat"]) - 1.0) < 1e-6
    for r in rows[1:]:
        f = int(r["Frame"])
        images.append({"unite": u, "systeme": pdb, "regime": regime, "image": f, **{m: r.get(m) for m in METR}, **{f"n_{t}": cnt[f].get(t, 0) for t in sorted(itypes)}})
    sub = rows[1:]
    med = lambda m: round(statistics.median(float(r[m]) for r in sub), 4)
    resume = {m: {"mediane": med(m), "min": round(min(float(r[m]) for r in sub), 4), "max": round(max(float(r[m]) for r in sub), 4)} for m in ("DockQ", "iRMSD", "lRMSD", "fnat", "clashes")}
    cm = {t: round(sum(cnt[int(r["Frame"])].get(t, 0) for r in sub) / len(sub), 2) for t in sorted(itypes)}
    par_unite.append({"unite": u, "systeme": pdb, "regime": regime, "dynapin_code": code, "garde_empilement_image0_reference": garde, "n_images": len(sub),
                      "contacts_moyens_par_image": cm, **{f"{m}_{k}": v for m, d_ in resume.items() for k, v in d_.items()}})
    ligne.update(statut="ok" if garde else "echec", garde_empilement_image0_reference=garde, n_images=len(sub), resume_metriques=resume, contacts_moyens_par_image=cm,
                 tables=str(d / "dynapin/tables"))
    if not garde:
        echecs.append({"unite": u, "raison": f"garde d empilement : image 0 DockQ {z['DockQ']} fnat {z['fnat']}"})
        ligne["erreur"] = "image 0 n est pas la reference"
    L.ajouter_manifeste(ligne)
tl = sorted(itypes)
with open(OUT / f"agregat_{US}_par_image.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["unite", "systeme", "regime", "image"] + METR + [f"n_{t}" for t in tl])
    for r in images:
        w.writerow([r["unite"], r["systeme"], r["regime"], r["image"]] + [r[m] for m in METR] + [r.get(f"n_{t}", 0) for t in tl])
with open(OUT / f"agregat_{US}_par_unite.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    cols = ["unite", "systeme", "regime", "dynapin_code", "garde_empilement_image0_reference", "n_images"] + [f"{m}_{k}" for m in ("DockQ", "iRMSD", "lRMSD", "fnat", "clashes") for k in ("mediane", "min", "max")] + [f"contacts_moyens_{t}" for t in tl]
    w.writerow(cols)
    for r in par_unite:
        w.writerow([r[c] for c in cols[:-len(tl)]] + [r["contacts_moyens_par_image"].get(t, 0) for t in tl])
(OUT / "bilan.json").write_text(json.dumps({"n_piles": len(unites), "n_ok": len(par_unite) - sum(1 for e in echecs if e["unite"] in {p["unite"] for p in par_unite}),
                                           "n_echecs": len(echecs), "echecs": echecs, "itypes": tl, "n_images_agregees": len(images)}, indent=1, ensure_ascii=False) + "\n")
print(f"H_AGREGER {len(unites)} piles, {len(par_unite)} agregees, {len(images)} images, types {tl}, echecs = {len(echecs)} {echecs[:5]}")
sys.exit(0 if not echecs else 1)
