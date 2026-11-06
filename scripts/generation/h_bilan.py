#!/usr/bin/env python3
"""h_bilan.py <racine_out> <run> — UNE ligne `etage_H` par systeme au manifeste : statut, images, mur, debit, signature, ecart des lourds restitues, concurrence
pendant CE systeme (fenetre du sampler), empreintes de l outil. Un systeme sans resume.json est un ECHEC, avec sa premiere ligne ARRET/ECHEC."""
import json, os, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L

racine, run = Path(sys.argv[1]), Path(sys.argv[2])
Hd = Path("/gpfs/phherman/bioemu2_bench/hydrogenes")
HR, US = os.environ.get("H_RACINE", "hydrogenes"), os.environ.get("H_USAGE", "U3")
sha_outil = dict((l.split()[1], l.split()[0]) for l in (run / "sha_outil.txt").read_text().splitlines() if l.strip())
temps = {l.split("\t")[0]: l.split("\t")[1:] for l in (run / "systemes_temps.tsv").read_text().splitlines() if l.strip()}
n_ok = n_ko = 0
for pdb in sorted(temps):
    t0, t1, code = float(temps[pdb][0]), float(temps[pdb][1]), int(temps[pdb][2])
    d = racine / pdb
    liste = [l for l in (L.G / HR / "listes" / f"{pdb}.txt").read_text().splitlines() if l.strip()]
    log = (racine / f"{pdb}.log").read_text() if (racine / f"{pdb}.log").is_file() else ""
    m = re.search(r"(ARRET|ECHEC)[^\n]*", log)
    res = json.loads((d / "resume.json").read_text()) if (d / "resume.json").is_file() else None
    r3 = list(d.glob("H_t10_i2000/R3/*/img_*.pdb"))
    ok = code == 0 and res is not None and len(r3) == len(liste) and res.get("images") == len(liste)
    ligne = {"type": "etage_H", "unite": f"{pdb}__{US}__etage_H", "systeme": pdb, "usage": US, "regime": "reconstruction_H_gpu",
             "statut": "ok" if ok else "echec", "code_retour_outil": code, "premiere_erreur": m.group(0)[:200] if m else None,
             "images_attendues": len(liste), "images_R3_sur_disque": len(r3), "mur_job_s": round(t1 - t0, 1), "jobid": run.name.rsplit("_", 1)[-1],
             "sortie": str(d), "jobs_concurrents_gpu": L.concurrence_fenetre(t0, t1, run), "sha256_outil_16": sha_outil,
             "outil": "hydrogenes/scripts/reconstruire_H_gpu.py (restitution ORDINALE ; JAMAIS L0c/restituer_lourds.py)"}
    if res:
        ligne.update({k: res.get(k) for k in ("images", "workers", "mur_s", "s_par_image_effectif", "atomes_par_image", "H_par_image", "signature_unique",
                                               "variantes", "ecart_lourds_restitues_max_A", "precision", "tol", "maxit", "cuda_visible_devices", "temps_moyens_par_image_s")})
    L.ajouter_manifeste(ligne)
    n_ok += ok; n_ko += (not ok)
    print(f"H_BILAN {pdb} {ligne['statut']} images {len(r3)}/{len(liste)} mur {ligne['mur_job_s']} s H/image {ligne.get('H_par_image')} ecart {ligne.get('ecart_lourds_restitues_max_A')} {ligne['premiere_erreur'] or ''}")
print(f"H_BILAN_FIN ok={n_ok} echec={n_ko} (rejets = {n_ko})")
sys.exit(0 if n_ko == 0 else 1)
