"""Bibliotheque commune du corpus DynaBench / BioEmu2.

Machine : staging-server. Dossier d'ecriture unique : /gpfs/phherman/bioemu2_bench/corpus/.
Lecture seule sur /gpfs/phherman/bioemu2_bench/L0a_work/ (cache de l'API, deja interrogee).

Regles portees par ce fichier :
  - un rejet est une ENTREE, jamais un saut silencieux ;
  - tout rapport de chargement dit combien d'items sont rejetes ;
  - on ne recalcule rien qui existe deja : le cache L0a est lu tel quel.
"""
import hashlib
import json
import os
import time
import urllib.error
import urllib.request

ROOT = "/gpfs/phherman/bioemu2_bench/corpus"
L0A = "/gpfs/phherman/bioemu2_bench/L0a_work"
API = "https://dynarepo.inria.fr/api/rest/v1"
CACHE = ROOT + "/cache"
LOGS = ROOT + "/logs"
STAGING_INPUTS = "/gpfs/phherman/evomd-data/Trajectory_Inputs/Trajectory_Inputs"
STAGING12 = "1B6C 1BKD 1FQ1 1M10 1ZM4 2AJF 2C0L 2OT3 3K75 5JMO 5VNW 6CWG".split()

MAX_RES = 1030
MIN_PAIRS = 25
CONTACT_A = 5.0


def now_iso():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256_file(path, bloc=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bloc)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def jload(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _numpy_en_python(o):
    """Seuls les scalaires numpy sont convertis ; tout autre type inattendu leve (jamais str() en silence)."""
    if hasattr(o, "item") and getattr(o, "shape", None) == ():
        return o.item()
    raise TypeError("type non serialisable: %s" % type(o).__name__)


def jdump(obj, path):
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, default=_numpy_en_python)
    os.replace(tmp, path)


def http_get(url, dest=None, timeout=900, tries=3, pause=2.0):
    """Telecharge `url`. Si `dest` est None, hache en flux sans rien ecrire.

    Rend TOUJOURS un dict (jamais d'exception) : status, octets, sha256, secondes, erreur.
    Un echec est une donnee.
    """
    last = None
    for t in range(1, tries + 1):
        t0 = time.time()
        h = hashlib.sha256()
        n = 0
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tlr4-corpus/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                status = r.status
                out = open(dest + ".part", "wb") if dest else None
                try:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        h.update(b)
                        n += len(b)
                        if out:
                            out.write(b)
                finally:
                    if out:
                        out.close()
            if dest:
                os.replace(dest + ".part", dest)
            return {"url": url, "status": status, "octets": n, "sha256": h.hexdigest(),
                    "secondes": round(time.time() - t0, 3), "essai": t, "erreur": None}
        except urllib.error.HTTPError as e:
            last = {"url": url, "status": e.code, "octets": n, "sha256": None,
                    "secondes": round(time.time() - t0, 3), "essai": t, "erreur": "HTTPError %s" % e.code}
            if e.code in (400, 401, 403, 404):
                return last
        except Exception as e:  # noqa: BLE001 - l'erreur est consignee, pas avalee
            last = {"url": url, "status": None, "octets": n, "sha256": None,
                    "secondes": round(time.time() - t0, 3), "essai": t,
                    "erreur": "%s: %s" % (type(e).__name__, e)}
        time.sleep(pause * t)
    return last


def declared_length(acc, n, filename):
    """Taille declaree par le serveur (cache L0a, aucune requete)."""
    p = "%s/meta_all/%s.%d.filenotes.json" % (L0A, acc, n)
    if not os.path.exists(p):
        return None
    for x in jload(p):
        if x.get("filename") == filename:
            return x.get("length")
    return None


def charger_projets_dynabench():
    """Les 210 projets DynaBench lus dans le cache de l'API (aucune requete).

    Retourne un dict : projets (liste), n_total_cache, n_dynabench, sha256 du cache.
    """
    p = L0A + "/all_projects.json"
    allp = jload(p)
    dyn = [q for q in allp if q["metadata"].get("COLLECTIONS") == ["dynabench"]]
    return {"fichier": p, "sha256": sha256_file(p), "n_total_cache": len(allp),
            "n_dynabench": len(dyn), "projets": dyn}
