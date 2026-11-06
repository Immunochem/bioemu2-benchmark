#!/usr/bin/env python3
"""msa_serie.py — les MSA ColabFold des 28 systemes, UNE REQUETE A LA FOIS (docs/msa.md : « query serially from a single IP »).

DERIVE de `msa/scripts/chercher_msa.py` (qui ecrit sous msa/ et exige BIOEMU2_MSA_CACHE == msa/cache : on ne peut pas
l utiliser sans ecrire dans un dossier interdit). On reprend SON instrumentation, ligne pour ligne dans l esprit :
  - un flock exclusif sur generation/msa/.verrou_serie AVANT toute requete : un second processus meurt (code 40) ;
  - une pause de --pause s avant chaque requete sauf la premiere du processus ;
  - un chien de garde de --timeout s par requete (BaseException : tenacity ne la rattrape pas) ;
  - aucune thread, aucun pool ; un echec est une ENTREE (statut=ecarte), jamais un saut ; une passe de reprise unique.
ECARTS DECLARES : (1) cache = generation/msa_cache, AMORCE par copie des caches existants (msa/cache, 6CWG/msa_cache,
vram_1030/msa_cache) — jamais ecrit dans l original ; (2) systemes lus sur les x0 CURES (gen_lib.construire_systeme), pas sur
msa/<PDB>/sequences.json ; (3) une passe de CLASSIFICATION d abord, sans reseau (run_colabfold_msa remplace par une
fonction qui leve) : elle dit, systeme par systeme, « dans le cache » ou « a chercher » ; (4) ordre des requetes :
les systemes de la courbe de capacite d abord (--priorite), puis tokens croissants.
Ce script N ECRIT PAS le MSA : get_msa(System) tel quel (memes requetes, memes cles de cache, memes defauts que BioEmu2).
"""
import argparse, datetime, fcntl, json, os, signal, socket, sys, time, traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L

ap = argparse.ArgumentParser()
ap.add_argument("--pause", type=float, default=20.0)
ap.add_argument("--timeout", type=int, default=900)
ap.add_argument("--reprise-pause", type=float, default=120.0)
ap.add_argument("--priorite", nargs="*", default=["1GCQ", "5E5M", "3RJQ"])
ap.add_argument("--classifier-seulement", action="store_true", help="classification sans reseau, puis sortie : aucune requete")
a = ap.parse_args()

M = L.G / "msa"
M.mkdir(exist_ok=True)
if os.environ.get("BIOEMU2_MSA_CACHE", "") != str(L.G / "msa_cache"):
    sys.exit(f"ARRET : BIOEMU2_MSA_CACHE doit valoir {L.G / 'msa_cache'}, valeur = {os.environ.get('BIOEMU2_MSA_CACHE')!r}")
LOCK = open(M / ".verrou_serie", "a+")
try:
    fcntl.flock(LOCK, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError:
    sys.exit("ARRET 40 : un autre processus MSA tient le verrou serie (une seule requete a la fois, jamais en parallele).")
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
    ent = {"n_dans_le_processus": ETAT["n"], "systeme": ETAT["systeme"], "job": os.environ.get("SLURM_JOB_ID"),
           "mode": "pair" if kw.get("use_pairing") else "msa", "use_env": kw.get("use_env"),
           "pairing_strategy": kw.get("pairing_strategy"), "longueurs": [len(s) for s in sequences], "t0": maintenant()}
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


def traiter(pdb, tentative):
    d = M / pdb
    d.mkdir(exist_ok=True)
    rec = {"pdb": pdb, "tentative": tentative, "job": os.environ.get("SLURM_JOB_ID"), "hote": socket.gethostname(), "date": maintenant(),
           "cache": os.environ["BIOEMU2_MSA_CACHE"]}
    try:
        system, seqs, tokens, x0info = L.construire_systeme(pdb)
    except L.Rejet as e:
        rec.update(statut="ecarte", raison=f"x0/systeme : {e}")
        (d / "profondeur.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n"); return rec
    rec.update(tokens=tokens, longueurs=[len(s) for s in seqs], x0=x0info["x0"], sha256_x0=x0info["sha256_x0"])
    ETAT["systeme"] = pdb
    n0 = len(JOURNAL_PROCESSUS)
    avant = set(os.listdir(os.environ["BIOEMU2_MSA_CACHE"]))
    t0 = time.time()
    try:
        msa, prof = L.msa_et_profil(system)
    except (Exception, ChienDeGarde) as e:
        rec.update(statut="ecarte", raison=f"{type(e).__name__}: {str(e)[:400]}", trace=traceback.format_exc()[-1200:],
                   requetes=JOURNAL_PROCESSUS[n0:], temps_get_msa_s=round(time.time() - t0, 1))
        (d / "profondeur.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n"); return rec
    rec.update(statut="ok", raison=None, chaines=prof, temps_get_msa_s=round(time.time() - t0, 1),
               requetes=JOURNAL_PROCESSUS[n0:], depuis_cache_integralement=len(JOURNAL_PROCESSUS[n0:]) == 0,
               fichiers_cache_nouveaux=sorted(set(os.listdir(os.environ["BIOEMU2_MSA_CACHE"])) - avant),
               note="profondeur = lignes d en-tete '>' du A3M, requete comprise ; apparie=None : une seule sequence unique")
    (d / "profondeur.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n")
    return rec


# ── 1. classification SANS reseau ────────────────────────────────────────────
cf.run_colabfold_msa = lambda *x, **k: (_ for _ in ()).throw(L.ReseauInterdit("classification hors reseau"))
tous = L.corpus_28()
tok = L.tokens_attendus()
classes = {"cache": [], "a_chercher": [], "ecarte": {}}
for pdb in sorted(tous):
    try:
        system, seqs, tokens, _ = L.construire_systeme(pdb)
        get_msa(system)
        classes["cache"].append(pdb)
    except L.ReseauInterdit:
        classes["a_chercher"].append(pdb)
    except L.Rejet as e:
        classes["ecarte"][pdb] = str(e)
    except Exception as e:
        classes["ecarte"][pdb] = f"{type(e).__name__}: {str(e)[:300]}"
print(f"CLASSIFICATION : {len(classes['cache'])} dans le cache, {len(classes['a_chercher'])} a chercher, {len(classes['ecarte'])} ecartes", flush=True)
print("  dans le cache :", classes["cache"], flush=True)
print("  a chercher    :", classes["a_chercher"], flush=True)
print("  ecartes       :", classes["ecarte"], flush=True)
(M / "classification.json").write_text(json.dumps(classes, indent=2, ensure_ascii=False) + "\n")
if a.classifier_seulement:
    print("CLASSIFICATION SEULE : aucune requete, aucun profil ecrit")
    sys.exit(0 if not classes["ecarte"] else 1)

# ── 2. les systemes deja en cache : profil, aucune requete ───────────────────
cf.run_colabfold_msa = run_instrumente
bilan = {"ok": [], "ecartes": dict(classes["ecarte"])}
for pdb in classes["cache"]:
    r = traiter(pdb, 1)
    if r["statut"] == "ok" and not r["depuis_cache_integralement"]:
        bilan["ecartes"][pdb] = "classe « dans le cache » mais des requetes ont ete emises : incoherence"
    (bilan["ok"].append(pdb) if r["statut"] == "ok" and r["depuis_cache_integralement"] else bilan["ecartes"].setdefault(pdb, r.get("raison")))
print(f"CACHE : {len(bilan['ok'])} profils ecrits sans requete", flush=True)

# ── 3. les autres, en SERIE ──────────────────────────────────────────────────
ordre = [p for p in a.priorite if p in classes["a_chercher"]] + sorted((p for p in classes["a_chercher"] if p not in a.priorite), key=lambda p: tok[p])
print("ORDRE DES REQUETES :", ordre, flush=True)
echecs = []
for pdb in ordre:
    t0 = time.time()
    r = traiter(pdb, 1)
    print(f"MSA {pdb} statut={r['statut']} t={time.time() - t0:.1f}s requetes={[(q['mode'], q['duree_s'], q['statut']) for q in r.get('requetes', [])]} raison={r.get('raison')}", flush=True)
    (bilan["ok"].append(pdb) if r["statut"] == "ok" else echecs.append(pdb))
if echecs:
    print(f"REPRISE (une seule passe) de {echecs} apres {a.reprise_pause} s", flush=True)
    time.sleep(a.reprise_pause)
    for pdb in list(echecs):
        r = traiter(pdb, 2)
        print(f"MSA(reprise) {pdb} statut={r['statut']} raison={r.get('raison')}", flush=True)
        if r["statut"] == "ok":
            bilan["ok"].append(pdb); echecs.remove(pdb)
        else:
            bilan["ecartes"][pdb] = r.get("raison")
n_rej = len(bilan["ecartes"])
(M / "bilan.json").write_text(json.dumps({"ok": sorted(bilan["ok"]), "ecartes": bilan["ecartes"], "n_ok": len(bilan["ok"]), "n_ecartes": n_rej,
                                          "n_requetes_dans_ce_processus": ETAT["n"], "fin": maintenant()}, indent=2, ensure_ascii=False) + "\n")
print(f"MSA_BILAN ok={len(bilan['ok'])}/28 ecartes={n_rej} (rejets = {n_rej}) ; requetes envoyees dans ce processus = {ETAT['n']} ; fin {maintenant()}", flush=True)
sys.exit(0 if n_rej == 0 and len(bilan["ok"]) == 28 else 1)
