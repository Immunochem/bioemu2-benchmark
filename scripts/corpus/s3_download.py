"""J3 — telechargement des 3 repliques (xtc + pdb + psf) des systemes retenus.

Protocole d'integrite (le serveur ne publie AUCUN hash) :
  passe 1  -> ecrit le fichier, octets compares a la taille DECLAREE par le serveur (filenotes) ;
  passe 2  -> second telechargement, haché en flux sans rien ecrire ;
  egal     = sha256(passe1) == sha256(passe2) ET octets == declare ;
  sinon    -> passe 3 ; s'il reste incoherent : verifie=False (donnee, jamais un echec silencieux).
Le numero `.n` de l'API est la POSITION dans `mds` (piege L0a : 4M5Z = repliques 1 et 3) ;
le dossier porte le numero de reference « replica k/3 », pas `.n`.
"""
import json
import os
import re
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import (API, CACHE, L0A, ROOT, charger_projets_dynabench, declared_length,  # noqa: E402
                        http_get, jload, now_iso, sha256_file)

FICHIERS = ["trajectory.xtc", "structure.pdb", "topology.psf"]
VERROU = threading.Lock()
LOG = CACHE + "/download_log.jsonl"


def deja_verifie():
    ok = {}
    if os.path.exists(LOG):
        for line in open(LOG, encoding="utf-8"):
            r = json.loads(line)
            if r.get("verifie"):
                ok[r["chemin"]] = r
    return ok


def tache(t, deja):
    acc, pdb, n, k, f = t["acc"], t["pdb"], t["n"], t["k"], t["f"]
    d = "%s/systemes/%s/rep%d" % (ROOT, pdb, k)
    os.makedirs(d, exist_ok=True)
    dest = "%s/%s" % (d, f)
    decl = declared_length(acc, n, f)
    url = "%s/projects/%s.%d/files/%s" % (API, acc, n, f)
    rec = {"pdb": pdb, "accession": acc, "api_index": n, "replique": k, "fichier": f, "chemin": dest,
           "declare_serveur": decl, "url": url}
    if dest in deja and os.path.exists(dest) and os.path.getsize(dest) == decl and sha256_file(dest) == deja[dest]["sha256"]:
        rec.update({k2: deja[dest][k2] for k2 in ("octets", "sha256", "sha256_2", "sha256_3", "egal", "verifie", "date")})
        rec["reprise"] = True
        return rec
    r1 = http_get(url, dest=dest)
    r2 = http_get(url, dest=None)
    rec.update({"octets": r1["octets"], "sha256": r1["sha256"], "sha256_2": r2["sha256"], "sha256_3": None,
                "statut_1": r1["status"], "statut_2": r2["status"], "secondes_1": r1["secondes"], "secondes_2": r2["secondes"],
                "debit_Mo_s_1": round(r1["octets"] / 1e6 / r1["secondes"], 2) if r1["secondes"] else None,
                "erreur_1": r1["erreur"], "erreur_2": r2["erreur"], "reprise": False})
    egal = (r1["status"] == 200 and r2["status"] == 200 and r1["sha256"] == r2["sha256"]
            and decl is not None and r1["octets"] == decl and r2["octets"] == decl)
    if not egal:
        r3 = http_get(url, dest=None)
        rec["sha256_3"] = r3["sha256"]
        rec["note"] = "passes 1/2 incoherentes (sha egal=%s, octets=%s/%s, declare=%s) ; passe 3 sha egal a la passe 1: %s, a la passe 2: %s" % (
            r1["sha256"] == r2["sha256"], r1["octets"], r2["octets"], decl, r3["sha256"] == r1["sha256"], r3["sha256"] == r2["sha256"])
        if r3["status"] == 200 and r3["octets"] == decl and r3["sha256"] == r1["sha256"] and r1["octets"] == decl:
            egal = True
        elif r3["status"] == 200 and r3["octets"] == decl and r3["sha256"] == r2["sha256"]:
            r4 = http_get(url, dest=dest)      # la copie ecrite etait la mauvaise : on la remplace
            egal = (r4["sha256"] == r3["sha256"] and r4["octets"] == decl)
            rec["sha256"] = r4["sha256"]
            rec["octets"] = r4["octets"]
    rec["egal"] = egal
    rec["verifie"] = bool(egal and os.path.exists(dest) and os.path.getsize(dest) == decl)
    rec["date"] = now_iso()
    return rec


def main():
    sel = jload(ROOT + os.environ.get("CORPUS_SELECTION", "/selection.json"))["selection"]
    projs = {q["accession"]: q for q in charger_projets_dynabench()["projets"]}
    taches = []
    for s in sel:
        q = projs[s["accession"]]
        for n, name in enumerate(q["mds"], start=1):
            m = re.match(r"replica (\d+)/(\d+)", name)
            k = int(m.group(1)) if m else n
            if not m:
                print("ALERTE mdName inattendu %s %s: %r (dossier rep%d)" % (s["pdb"], s["accession"], name, k), flush=True)
            for f in FICHIERS:
                taches.append({"acc": s["accession"], "pdb": s["pdb"], "n": n, "k": k, "f": f, "mdName": name})
    tot_decl = sum((declared_length(t["acc"], t["n"], t["f"]) or 0) for t in taches)
    print("TACHES %d fichiers pour %d systemes ; octets declares=%d (%.2f Go) ; x2 telechargements=%.2f Go" % (
        len(taches), len(sel), tot_decl, tot_decl / 1e9, 2 * tot_decl / 1e9), flush=True)
    deja = deja_verifie()
    t0 = time.time()
    ok = ko = 0
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs = [ex.submit(tache, t, deja) for t in taches]
        with open(LOG, "a", encoding="utf-8") as fo:
            for i, fu in enumerate(futs, start=1):
                rec = fu.result()
                with VERROU:
                    fo.write(json.dumps(rec, ensure_ascii=False) + "\n")
                    fo.flush()
                if rec["verifie"]:
                    ok += 1
                else:
                    ko += 1
                print("DL %3d/%d %s rep%d %-15s octets=%s declare=%s egal=%s verifie=%s reprise=%s debit=%s Mo/s %s" % (
                    i, len(taches), rec["pdb"], rec["replique"], rec["fichier"], rec.get("octets"), rec["declare_serveur"],
                    rec.get("egal"), rec["verifie"], rec.get("reprise"), rec.get("debit_Mo_s_1"), rec.get("note", "")), flush=True)
    # filenotes (provenance, aucune requete : copie du cache L0a)
    for s in sel:
        q = projs[s["accession"]]
        for n, name in enumerate(q["mds"], start=1):
            m = re.match(r"replica (\d+)/(\d+)", name)
            k = int(m.group(1)) if m else n
            src = "%s/meta_all/%s.%d.filenotes.json" % (L0A, s["accession"], n)
            dst = "%s/systemes/%s/rep%d/filenotes.json" % (ROOT, s["pdb"], k)
            shutil.copyfile(src, dst)
    dt = time.time() - t0
    print("BILAN fichiers_verifies=%d  non_verifies=%d  secondes=%.1f" % (ok, ko, dt), flush=True)
    print("FINI s3 %s" % now_iso(), flush=True)
    sys.exit(0 if ko == 0 else 4)


if __name__ == "__main__":
    main()
