#!/usr/bin/env python3
"""redock_msa.py — MSA ColabFold du redocking, UNE REQUETE A LA FOIS (docs/msa.md : « query serially from a single IP »).

DERIVE de generation/scripts/msa_serie.py (meme instrumentation, memes defauts de get_msa) ; ECARTS DECLARES : (1) la liste des systemes vient d un fichier d attente
(msa_attente.json : id + sequences, dans l ordre de priorite) au lieu du corpus de 28 ; (2) cache = redocking/msa_cache (amorce par COPIE des caches existants, jamais ecrit dans
les originaux) ; (3) journal et profils sous redocking/msa/.
  - un flock exclusif sur redocking/msa/.verrou_serie avant toute requete : un second processus meurt (code 40) ;
  - une pause de --pause s avant chaque requete sauf la premiere du processus ; un chien de garde de --timeout s par requete ;
  - aucune thread, aucun pool ; un echec est une ENTREE (statut=ecarte), jamais un saut ; une passe de reprise unique ;
  - une passe de CLASSIFICATION d abord, sans reseau : « dans le cache » / « a chercher ».
Ce script N ECRIT PAS le MSA : get_msa(System) tel quel (memes requetes, memes cles de cache, memes defauts que BioEmu2)."""
import argparse, datetime, fcntl, hashlib, json, os, signal, socket, sys, time, traceback
from pathlib import Path

sys.dont_write_bytecode = True
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
ap = argparse.ArgumentParser()
ap.add_argument("--attente", default=str(R / "msa/msa_attente.json"))
ap.add_argument("--pause", type=float, default=20.0)
ap.add_argument("--timeout", type=int, default=900)
ap.add_argument("--reprise-pause", type=float, default=120.0)
ap.add_argument("--max", type=int, default=0, help="n' envoyer des requetes que pour les N premiers systemes a chercher (0 = tous)")
ap.add_argument("--classifier-seulement", action="store_true")
a = ap.parse_args()

M = R / "msa"
M.mkdir(exist_ok=True)
if os.environ.get("BIOEMU2_MSA_CACHE", "") != str(R / "msa_cache"):
    sys.exit(f"ARRET : BIOEMU2_MSA_CACHE doit valoir {R / 'msa_cache'}, valeur = {os.environ.get('BIOEMU2_MSA_CACHE')!r}")
LOCK = open(M / ".verrou_serie", "a+")
try:
    fcntl.flock(LOCK, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError:
    sys.exit("ARRET 40 : un autre processus MSA du redocking tient le verrou serie.")
LOCK.seek(0); LOCK.truncate()
LOCK.write(f"{os.getpid()} job={os.environ.get('SLURM_JOB_ID')} {datetime.datetime.now().astimezone().isoformat()}\n"); LOCK.flush()

import bioemu2.msa.colabfold as cf
from bioemu2 import System, get_msa

JOURNAL = M / "journal_requetes.jsonl"
ETAT = {"n": 0, "systeme": None}
ORIG = cf.run_colabfold_msa
JOURNAL_PROCESSUS = []


class ChienDeGarde(BaseException):
    pass


class ReseauInterdit(Exception):
    pass


def maintenant():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def _alarme(signum, frame):
    raise ChienDeGarde(f"chien de garde : requete > {a.timeout} s")


def _terme(signum, frame):
    with open(JOURNAL, "a") as f:
        f.write(json.dumps({"evenement": "signal", "signal": signum, "systeme": ETAT["systeme"], "t": maintenant()}) + "\n")
    sys.exit(143)


signal.signal(signal.SIGALRM, _alarme)
signal.signal(signal.SIGTERM, _terme)


def run_instrumente(sequences, **kw):
    ETAT["n"] += 1
    if ETAT["n"] > 1:
        time.sleep(a.pause)
    ent = {"n_dans_le_processus": ETAT["n"], "systeme": ETAT["systeme"], "job": os.environ.get("SLURM_JOB_ID"), "mode": "pair" if kw.get("use_pairing") else "msa",
           "use_env": kw.get("use_env"), "pairing_strategy": kw.get("pairing_strategy"), "longueurs": [len(s) for s in sequences], "t0": maintenant()}
    t0 = time.time()
    signal.alarm(int(a.timeout))
    try:
        res = ORIG(sequences, **kw)
        ent.update(statut="ok", n_a3m=len(res), octets=[len(r) for r in res])
        return res
    except BaseException as e:
        ent.update(statut="echec", erreur=f"{type(e).__name__}: {str(e)[:300]}")
        raise
    finally:
        signal.alarm(0)
        ent["duree_s"] = round(time.time() - t0, 2); ent["t1"] = maintenant()
        JOURNAL_PROCESSUS.append(ent)
        with open(JOURNAL, "a") as f:
            f.write(json.dumps(ent, ensure_ascii=False) + "\n"); f.flush(); os.fsync(f.fileno())


def entetes(a3m):
    return sum(1 for l in a3m.splitlines() if l.startswith(">"))


def profil(msa):
    prof = {}
    for c in msa.chain_ids:
        u = msa.a3m(c)
        d = {"non_apparie": entetes(u), "sha256_non_apparie": hashlib.sha256(u.encode()).hexdigest()[:16]}
        try:
            p = msa.a3m(c, paired=True)
            d["apparie"] = entetes(p); d["sha256_apparie"] = hashlib.sha256(p.encode()).hexdigest()[:16]
        except ValueError:
            d["apparie"] = None
        prof[c] = d
    return prof


def traiter(ent, tentative):
    ETAT["systeme"] = ent["id"]
    rec = {"id": ent["id"], "tentative": tentative, "job": os.environ.get("SLURM_JOB_ID"), "hote": socket.gethostname(), "date": maintenant(), "longueurs": [len(s) for s in ent["sequences"]]}
    n0 = len(JOURNAL_PROCESSUS)
    t0 = time.time()
    try:
        system = System.from_description(proteins=ent["sequences"])
        msa = get_msa(system)
        rec.update(statut="ok", raison=None, chaines=profil(msa), temps_get_msa_s=round(time.time() - t0, 1), requetes=JOURNAL_PROCESSUS[n0:],
                   depuis_cache_integralement=len(JOURNAL_PROCESSUS[n0:]) == 0)
    except (Exception, ChienDeGarde) as e:
        rec.update(statut="ecarte", raison=f"{type(e).__name__}: {str(e)[:400]}", trace=traceback.format_exc()[-1200:], requetes=JOURNAL_PROCESSUS[n0:], temps_get_msa_s=round(time.time() - t0, 1))
    (M / f"{ent['id']}.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n")
    return rec


attente = json.loads(Path(a.attente).read_text())
# ── 1. classification SANS reseau
def _interdit(*x, **k):
    raise ReseauInterdit("classification hors reseau")
cf.run_colabfold_msa = _interdit
cache, chercher, ecarte = [], [], {}
for ent in attente:
    try:
        get_msa(System.from_description(proteins=ent["sequences"]))
        cache.append(ent["id"])
    except ReseauInterdit:
        chercher.append(ent)
    except Exception as e:
        ecarte[ent["id"]] = f"{type(e).__name__}: {str(e)[:300]}"
print(f"CLASSIFICATION : {len(cache)} dans le cache, {len(chercher)} a chercher, {len(ecarte)} ecartes (sur {len(attente)})", flush=True)
(M / "classification.json").write_text(json.dumps({"cache": cache, "a_chercher": [e["id"] for e in chercher], "ecarte": ecarte}, indent=1) + "\n")
if a.classifier_seulement:
    sys.exit(0)
cf.run_colabfold_msa = run_instrumente
bilan = {"ok": [], "ecartes": dict(ecarte)}
# ── 2. les systemes deja en cache : profil, aucune requete
for ent in attente:
    if ent["id"] in cache:
        r = traiter(ent, 1)
        if r["statut"] == "ok" and r["depuis_cache_integralement"]:
            bilan["ok"].append(ent["id"])
        else:
            bilan["ecartes"][ent["id"]] = r.get("raison") or "classe « dans le cache » mais des requetes ont ete emises"
print(f"CACHE : {len(bilan['ok'])} profils ecrits sans requete", flush=True)
# ── 3. les autres, en SERIE, dans l ordre de la file d attente
if a.max:
    chercher = chercher[:a.max]
echecs = []
for ent in chercher:
    t0 = time.time()
    r = traiter(ent, 1)
    print(f"MSA {ent['id']} statut={r['statut']} t={time.time() - t0:.1f}s requetes={[(q['mode'], q['duree_s'], q['statut']) for q in r.get('requetes', [])]} raison={r.get('raison')}", flush=True)
    (bilan["ok"].append(ent["id"]) if r["statut"] == "ok" else echecs.append(ent))
if echecs:
    print(f"REPRISE (une seule passe) de {[e['id'] for e in echecs]} apres {a.reprise_pause} s", flush=True)
    time.sleep(a.reprise_pause)
    for ent in list(echecs):
        r = traiter(ent, 2)
        print(f"MSA(reprise) {ent['id']} statut={r['statut']} raison={r.get('raison')}", flush=True)
        if r["statut"] == "ok":
            bilan["ok"].append(ent["id"]); echecs.remove(ent)
        else:
            bilan["ecartes"][ent["id"]] = r.get("raison")
n_rej = len(bilan["ecartes"])
(M / "bilan.json").write_text(json.dumps({"ok": sorted(bilan["ok"]), "ecartes": bilan["ecartes"], "n_ok": len(bilan["ok"]), "n_ecartes": n_rej,
                                          "n_requetes_dans_ce_processus": ETAT["n"], "fin": maintenant()}, indent=2, ensure_ascii=False) + "\n")
print(f"MSA_BILAN ok={len(bilan['ok'])}/{len(attente)} ecartes={n_rej} (rejets = {n_rej}) ; requetes envoyees dans ce processus = {ETAT['n']} ; fin {maintenant()}", flush=True)
sys.exit(0 if n_rej == 0 else 1)
