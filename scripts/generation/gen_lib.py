#!/usr/bin/env python3
"""gen_lib.py — bibliotheque commune de la generation BioEmu2 du benchmark (capacite, U1, U2, U3).

Ce qu elle NE reecrit PAS (importe tel quel) : `concurrence` et `arbre_sha` (L0b/outils), get_msa / EmuSession /
check_physicality / save_frames (bioemu2). Elle ne porte que ce qui manquait : le corpus, la lecture des x0 curés
(avec leur sha256), le manifeste append-only, la fenetre de concurrence par unite, le releve des pics memoire.

Sous-commande :  gen_lib.py fin_job --run DOSSIER --code N [--unite U]   (ajoute la ligne `fin_job` au manifeste)

🔴 Aucun `2>/dev/null`, aucune erreur convertie en valeur rassurante : un rejet est une ENTREE du manifeste.
"""
import fcntl, hashlib, json, os, resource, socket, sys, time, traceback, zlib
from pathlib import Path

BENCH = Path("/gpfs/phherman/bioemu2_bench")
G = BENCH / "generation"
CUR = BENCH / "prep_fix/curation/x0"
L0B = BENCH / "L0b/outils"
MANIF = Path(os.environ.get("GEN_MANIFESTE", str(G / "manifeste_generation.json")))   # JSON-lines append-only ; GEN_MANIFESTE : fumee seulement
SHA_UPSTREAM_ATTENDU = "06be369dd555849165a737472d42d936959e46fc"
STD = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
HORS_PERIMETRE = {"1ZM4": "DDE (histidine + adduit covalent de 10 atomes)", "1FQ1": "TPO (phosphothreonine)"}
GIO = 2 ** 30
sys.path.insert(0, str(L0B))


class Rejet(Exception):
    """Un rejet est une donnee : l appelant l ecrit au manifeste, jamais de saut silencieux."""


# ───────────────────────────── corpus et x0 ─────────────────────────────
def corpus_28():
    """28 systemes = statut actif_v2 du manifeste, moins 1ZM4 et 1FQ1 (residus non standard). Leve si le compte differe."""
    L = [json.loads(l) for l in open(BENCH / "corpus/manifeste_corpus.json") if l.strip()]
    A = [x for x in L if x.get("statut") == "actif_v2" and x["pdb"] not in HORS_PERIMETRE]
    pdbs = sorted(x["pdb"] for x in A)
    if len(pdbs) != 28 or len(set(pdbs)) != 28:
        raise SystemExit(f"ARRET : corpus = {len(pdbs)} systemes distincts {len(set(pdbs))}, 28 attendus")
    return {x["pdb"]: x for x in A}


def tokens_attendus():
    return json.loads((BENCH / "prep_fix/curation/tokens_vram.json").read_text())["tokens"]


def sha256_fichier(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def chemin_x0(pdb, rep=1):
    return CUR / f"{pdb}__img0_rep{rep}.pdb"


def verifier_x0(pdb, rep=1):
    """sha256 du x0 relu contre `sha256_x0` du JSON compagnon. Leve Rejet si absent ou different."""
    p, j = chemin_x0(pdb, rep), CUR / f"{pdb}__img0_rep{rep}.json"
    if not p.is_file() or not j.is_file():
        raise Rejet(f"x0 ou JSON compagnon absent : {p}")
    attendu = json.loads(j.read_text()).get("sha256_x0")
    reel = sha256_fichier(p)
    if attendu != reel:
        raise Rejet(f"sha256 du x0 {p.name} = {reel[:16]}, JSON compagnon dit {str(attendu)[:16]}")
    return {"x0": str(p), "sha256_x0": reel}


def construire_systeme(pdb, rep=1):
    """(system, seqs, tokens, x0_info). Meme chemin que vram_1zm4.py:46-53 (to_fasta -> System.from_description) et
    qu il refuse les memes residus non standard ; la lecture se fait sur le x0 CURE (contrat d entree §10)."""
    import mdtraj as md
    from bioemu2 import System
    if pdb in HORS_PERIMETRE:
        raise Rejet(f"{pdb} hors perimetre : {HORS_PERIMETRE[pdb]}")
    info = verifier_x0(pdb, rep)
    ref = md.load(info["x0"])
    nonstd = sorted({r.name for ch in ref.topology.chains for r in ch.residues if r.name not in STD})
    if nonstd:
        raise Rejet(f"residus non standard {nonstd} — hors perimetre BioEmu2")
    seqs = list(ref.topology.to_fasta())
    tokens = sum(len(s) for s in seqs)
    att = tokens_attendus().get(pdb)
    if att is not None and att != tokens:
        raise Rejet(f"tokens = {tokens}, tokens_vram.json dit {att}")
    return System.from_description(proteins=seqs), seqs, tokens, info


# ───────────────────────────── MSA ─────────────────────────────
class ReseauInterdit(RuntimeError):
    pass


def interdire_reseau_msa():
    """Hors du job MSA serie, AUCUNE requete ColabFold : un cache manquant leve au lieu d interroger le serveur public."""
    import bioemu2.msa.colabfold as cf

    def _interdit(*a, **k):
        raise ReseauInterdit("requete ColabFold interdite hors du job MSA serie (docs/msa.md : query serially)")
    cf.run_colabfold_msa = _interdit


def _entetes(a3m):
    return sum(1 for l in a3m.splitlines() if l.startswith(">"))


def msa_et_profil(system):
    """get_msa(system) tel quel (defauts use_env=True, greedy) + profondeur par chaine, non appariee ET appariee (G11)."""
    from bioemu2 import get_msa
    msa = get_msa(system)
    prof = {}
    for c in msa.chain_ids:
        u = msa.a3m(c)
        d = {"non_apparie": _entetes(u), "sha256_non_apparie": hashlib.sha256(u.encode()).hexdigest()[:16]}
        try:
            p = msa.a3m(c, paired=True)
            d["apparie"] = _entetes(p)
            d["sha256_apparie"] = hashlib.sha256(p.encode()).hexdigest()[:16]
        except ValueError:
            d["apparie"] = None
            d["apparie_note"] = "aucun alignement apparie (une seule sequence unique)"
        prof[c] = d
    return msa, prof


# ───────────────────────────── ponts disulfure ─────────────────────────────
REPERE_PONT_ROMPU_A = 3.0      # un REPERE (celui du pilote 6CWG), pas une verite : la reference est la plage mesuree dans la MD


def charger_ponts(pdb):
    """Ponts du systeme (ponts_reference.py). Leve Rejet si le releve n existe pas ou a echoue pour ce systeme."""
    f = G / "ponts/reference.json"
    if not f.is_file():
        raise Rejet("ponts/reference.json absent : lancer ponts_reference.py avant la generation")
    s = json.loads(f.read_text())["systemes"].get(pdb)
    if s is None or s.get("statut") != "ok":
        raise Rejet(f"releve des ponts de {pdb} absent ou ecarte : {None if s is None else s.get('raison')}")
    return s


def indices_sg(topology, ponts):
    """[(id, i, j)] : indices d atomes SG des deux cysteines de chaque pont, residus pris par ORDINAL dans la chaine (A, B…)."""
    chaines = {chr(ord("A") + k): ch for k, ch in enumerate(topology.chains)}
    out = []
    for p in ponts:
        ids = []
        for k in ("i", "j"):
            res = list(chaines[p[f"chaine_{k}"]].residues)[p[f"ordinal_{k}"] - 1]
            sg = [a.index for a in res.atoms if a.name == "SG"]
            if res.name != "CYS" or len(sg) != 1:
                raise Rejet(f"pont {p['id']} : le residu ordinal {p[f'ordinal_{k}']} de la chaine {p[f'chaine_{k}']} est {res.name} avec {len(sg)} SG")
            ids.append(sg[0])
        out.append((p["id"], ids[0], ids[1]))
    return out


def resume_ponts(traj, sg, ponts_ref, passed):
    """Distances SG–SG (A) pour chaque pont et chaque structure, et ce que check_physicality en dit sur les MEMES structures."""
    import numpy as np
    n = int(traj.n_frames)
    d, par, rompu_any = {}, {}, np.zeros(n, dtype=bool)
    pas = np.array([bool(x) for x in passed])
    for (pid, i, j), p in zip(sg, ponts_ref):
        dist = np.linalg.norm(traj.xyz[:, i] - traj.xyz[:, j], axis=1) * 10.0
        d[pid] = [round(float(x), 3) for x in dist]
        rompu = dist > REPERE_PONT_ROMPU_A
        rompu_any |= rompu
        hors = dist > p["plage_toutes_images_A"][1]
        par[pid] = {"mediane_A": round(float(np.median(dist)), 3), "min_A": round(float(dist.min()), 3), "max_A": round(float(dist.max()), 3),
                    "n_rompus_repere_3A": int(rompu.sum()), "n_au_dela_du_max_des_images_MD": int(hors.sum()),
                    "n_rompus_repere_3A_et_retenues_par_check_physicality": int((rompu & pas).sum()),
                    "ref_md_mediane_des_medianes_A": p["mediane_des_medianes_par_replique_A"],
                    "ref_md_plage_des_medianes_par_replique_A": p["plage_des_medianes_par_replique_A"],
                    "ref_md_plage_toutes_images_A": p["plage_toutes_images_A"], "ref_md_image0_rep1_A": p["md_img0_rep1_A"]}
    return d, {"n_ponts": len(sg), "repere_pont_rompu_A": REPERE_PONT_ROMPU_A, "par_pont": par,
               "n_structures_avec_un_pont_rompu": int(rompu_any.sum()),
               "n_structures_rompues_et_retenues_par_check_physicality": int((rompu_any & pas).sum()),
               "n_structures_retenues": int(pas.sum())}


# ───────────────────────────── memoire, temps, concurrence ─────────────────────────────
def rss_max_gio():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 20), 2)   # Linux : ko


def seed_unite(unite):
    return zlib.crc32(unite.encode()) % (2 ** 31)


def concurrence_instant():
    import concurrence
    return concurrence.un_echantillon("staging")


def concurrence_fenetre(t0, t1, run, intervalle=15.0):
    """min/max du nombre de calculs GPU echantillonnes (sampler de garde_gpu) sur [t0-intervalle, t1+intervalle]."""
    f = Path(run) / "concurrence.jsonl"
    vals, err = [], 0
    if f.is_file():
        for l in f.read_text().splitlines():
            if not l.strip():
                continue
            d = json.loads(l)
            if "erreur" in d:
                err += 1
            elif t0 - intervalle <= d["t"] <= t1 + intervalle:
                vals.append(d["total"])
    return {"min": min(vals) if vals else None, "max": max(vals) if vals else None,
            "n_echantillons": len(vals), "n_erreurs_echantillon": err,
            "note": "total = jobs SLURM RUNNING avec gres gpu (moi compris) + calculs hors slurm ; BORNE BASSE (concurrence.py)"}


def lire_fiche(run):
    p = Path(run) / "fiche.json"
    return json.loads(p.read_text()) if p.is_file() else {}


_ARBRE = {}


def empreinte_moteur():
    """sha256 de l arbre bioemu2 REELLEMENT installe (site-packages), independant de git (arbre_sha.py)."""
    if "v" not in _ARBRE:
        import importlib.util, arbre_sha
        site = os.path.dirname(importlib.util.find_spec("bioemu2").origin)
        a = arbre_sha.arbre(site, ".py")
        _ARBRE["v"] = {"sha256_arbre_installe": arbre_sha.empreinte(a), "n_fichiers_py": len(a), "site": site}
    return _ARBRE["v"]


def contexte_job(run):
    f = lire_fiche(run)
    gen = G / "scripts"
    return {
        "machine": "staging", "partition": os.environ.get("SLURM_JOB_PARTITION"), "jobid": os.environ.get("SLURM_JOB_ID"),
        "noeud": os.environ.get("SLURMD_NODENAME") or socket.gethostname(),
        "time_demande": f.get("time_demande"), "cpus_par_tache": os.environ.get("SLURM_CPUS_PER_TASK"),
        "mem_demandee_Mo_par_noeud": os.environ.get("SLURM_MEM_PER_NODE"),
        "sha_commit_moteur_execute": f.get("sha_commit_execute"), "sha_commit_moteur_attendu": f.get("sha_commit_attendu"),
        "sha256_arbre_moteur_installe": empreinte_moteur()["sha256_arbre_installe"],
        "sha256_contourner_nvml": (f.get("sha256_fichiers_moteur") or {}),
        "profil_gpu_vu": (f.get("gpu") or {}).get("profil_mig_vu"), "cuda_visible_devices": f.get("cuda_visible_devices"),
        "jobs_concurrents_prevus": f.get("jobs_concurrents_prevus"),
        "sha256_scripts": {n: sha256_fichier(gen / n)[:16] for n in ("gen_lib.py", "commun_gen.sh") if (gen / n).is_file()},
        "sha256_commun_gpu_6cwg": sha256_fichier(BENCH / "6CWG/scripts/commun_gpu.sh")[:16],
    }


# ───────────────────────────── manifeste append-only ─────────────────────────────
def ajouter_manifeste(e):
    """Une ligne JSON, flock exclusif, fsync. JAMAIS de reecriture d une ligne existante."""
    e = dict(e)
    e.setdefault("schema", "generation-bioemu2/1")
    e["ecrit_le_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    ligne = (json.dumps(e, ensure_ascii=False) + "\n").encode()
    fd = os.open(MANIF, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o664)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        os.write(fd, ligne)
        os.fsync(fd)
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


def trace_courte(n=900):
    return traceback.format_exc()[-n:]


# ───────────────────────────── fin de job ─────────────────────────────
def cmd_fin_job(run, code, unite):
    f = lire_fiche(run)
    lignes = [json.loads(l) for l in MANIF.read_text().splitlines() if l.strip()] if MANIF.is_file() else []
    jid = os.environ.get("SLURM_JOB_ID")
    miens = [l for l in lignes if l.get("jobid") == jid and l.get("type") in ("unite", "essai_capacite", "rejet")]
    comptes = {}
    for l in miens:
        comptes[l.get("statut")] = comptes.get(l.get("statut"), 0) + 1
    ajouter_manifeste({
        "type": "fin_job", "unite": unite, "jobid": jid, "code_retour_python": code,
        "time_demande": f.get("time_demande"), "elapsed_job_s": f.get("elapsed_job_s"),
        "jobs_concurrents_gpu": f.get("jobs_concurrents_gpu"),
        "depot_stable_pendant_le_job": f.get("depot_stable_pendant_le_job"),
        "n_lignes_du_job": len(miens), "statuts": comptes,
        "note": "un code de retour n est pas une mesure : verifier les artefacts (xtc, sha256, physicality.json)",
    })
    print(f"fin_job : {len(miens)} ligne(s) du job {jid} : {comptes}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    f = sp.add_parser("fin_job")
    f.add_argument("--run", required=True)
    f.add_argument("--code", type=int, required=True)
    f.add_argument("--unite", default="")
    a = ap.parse_args()
    cmd_fin_job(a.run, a.code, a.unite)
