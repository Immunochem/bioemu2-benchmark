#!/usr/bin/env python3
"""site_positifs.py <lag> <etiquette> [<etiquette> ...] — CPU. Pour chaque image d un systeme POSITIF : le « site » de la chaine 0 (A) = residus de la chaine A qui touchent la chaine B dans le CRISTAL (couples < 5 A) ;
rappel_site = part du site que la chaine B touche dans la sortie ; precision_site = part des residus de A touches qui sont dans le site ; n_contacts = couples de residus en contact.
Memes definitions que `redock_cond.mesurer_neg` (temoin negatif) -> comparaison APPARIEE, chaine A par chaine A : vrai partenaire contre etranger. Appariement ordinal (Reference). Ecrit resultats/site_positifs_lag<L>.json."""
import json, os, sys
from pathlib import Path
os.environ["CUDA_VISIBLE_DEVICES"] = ""
sys.dont_write_bytecode = True
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
sys.path.insert(0, str(R / "scripts"))
import numpy as np, mdtraj as md
import redock_lib as RL, redock_bioemu as RB
from redock_lib import Q

lag = int(sys.argv[1]); etiqs = sys.argv[2:]
out, rejets = {}, []
for e in etiqs:
    for d in sorted((R / "sorties" / e).glob("*/")):
        mf = d / "meta.json"
        if not mf.is_file():
            continue
        meta = json.loads(mf.read_text())
        if meta.get("statut") != "ok" or meta["kind"] != "pos" or not (d / f"lag{lag:04d}/samples.xtc").is_file():
            continue
        try:
            info = RB.systeme_depuis_x0(meta["x0"])
            Rf = RL.Reference(info["ref"])
            tr = md.load_xtc(str(d / f"lag{lag:04d}/samples.xtc"), top=str(d / f"lag{lag:04d}/topology.cif"))
            if [a.name for a in tr.topology.atoms] != [a.name for a in Rf.ref.topology.atoms]:
                raise ValueError("noms d atomes differents de la reference : appariement ordinal impossible")
            S = {p[0] for p in Rf.contacts}
            n_res = len(Rf.res[0])
            ct = Q.contacts_residus(tr, Q.residus_par_chaine(tr.topology))
            cond = json.loads((d / "conditions.json").read_text())["audit"]
            idx = {a["indice_image"]: a for a in cond if a["accepte"]}
            rows = []
            for i, c in enumerate(ct):
                touches = {p[0] for p in c}
                a = idx[i]
                rows.append({"cond": a["cond"], "rep": a["rep"], "n_contacts": len(c), "rappel_site": round(len(S & touches) / len(S), 4), "precision_site": round(len(S & touches) / len(touches), 4) if touches else None})
            out[meta["id"]] = {"strate": meta["strate"], "taille_site": len(S), "n_residus_chaine_A": n_res, "base_hasard_precision": len(S) / n_res, "images": rows, "etiquette": e}
        except Exception as ex:
            rejets.append({"id": meta["id"], "erreur": f"{type(ex).__name__}: {str(ex)[:200]}"})
(R / "resultats").mkdir(exist_ok=True)
(R / f"resultats/site_positifs_lag{lag:04d}.json").write_text(json.dumps({"lag": lag, "systemes": out, "rejets": rejets}) + "\n")
print(f"SITE_POSITIFS lag {lag} : {len(out)} systemes, {len(rejets)} rejets (REJETS = {len(rejets)}) {rejets[:3]}")
