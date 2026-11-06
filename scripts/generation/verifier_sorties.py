#!/usr/bin/env python3
"""verifier_sorties.py — relit les ARTEFACTS, pas les codes de retour : pour chaque dossier de sorties/, recalcule le sha256 de chaque fichier,
compte les images du xtc, relit physicality.json, et compare au manifeste (derniere ligne `unite` ok de cette unite). Ecrit verification_sorties.json.
Un dossier qui echoue est une ENTREE du bilan, jamais un saut."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L
import mdtraj as md

S = L.G / "sorties"
man = {}
for i, l in enumerate(L.MANIF.read_text().splitlines(), 1):
    if l.strip():
        d = json.loads(l)
        if d.get("type") == "unite" and d.get("statut") == "ok":
            man[d["unite"]] = (i, d)
bilan, anom = {"dossiers": 0, "ok": 0}, []
vus = set()
for d in sorted(S.iterdir()):
    if not d.is_dir():
        continue
    if d.name.endswith(".tmp") or ".perime_" in d.name:
        anom.append({"dossier": d.name, "probleme": "dossier .tmp ou perime present"}); continue
    bilan["dossiers"] += 1
    vus.add(d.name)
    probs = []
    try:
        meta = json.loads((d / "meta.json").read_text())
        for n, h in meta["sha256_sorties"].items():
            if L.sha256_fichier(d / n) != h:
                probs.append(f"sha256 de {n} different de meta.json")
        with md.formats.XTCTrajectoryFile(str(d / "samples.xtc")) as f:
            nf = len(f)
        ph = json.loads((d / "physicality.json").read_text())
        passed = sum(1 for x in ph if x["status"] == "passed")
        if nf != meta["n"]:
            probs.append(f"{nf} images dans le xtc, n = {meta['n']}")
        if len(ph) != nf or passed != meta["retenus"]:
            probs.append(f"physicality.json : {len(ph)} entrees / {passed} passed, meta dit {meta['retenus']}")
        m = man.get(d.name)
        if m is None:
            probs.append("aucune ligne `unite` ok dans le manifeste")
        elif m[1]["sha256_sorties"] != meta["sha256_sorties"]:
            probs.append(f"sha256 du manifeste (ligne {m[0]}) != meta.json")
    except Exception as e:
        probs.append(f"{type(e).__name__}: {str(e)[:200]}")
    if probs:
        anom.append({"dossier": d.name, "problemes": probs})
    else:
        bilan["ok"] += 1
manque = sorted(set(man) - vus)
bilan.update(anomalies=anom, n_anomalies=len(anom), unites_du_manifeste_sans_dossier=manque, n_unites_manifeste_ok=len(man))
(L.G / "verification_sorties.json").write_text(json.dumps(bilan, indent=1, ensure_ascii=False) + "\n")
print("VERIFICATION", json.dumps({k: v for k, v in bilan.items() if k != "anomalies"}, ensure_ascii=False), f"(anomalies = {len(anom)})")
for a in anom[:20]:
    print("  ", a)
sys.exit(0 if not anom and not manque else 1)
