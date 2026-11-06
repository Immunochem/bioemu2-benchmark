#!/usr/bin/env python3
"""generer.py — production BioEmu2 du benchmark : U3 (propagation), U2 (ensemble d equilibre), U1 (prediction de structure).

DERIVE de `6CWG/scripts/propager.py` (qu on ne modifie pas) : meme discipline — une EmuSession par mode (embedding calcule UNE
fois), x0 = image 0 de la replique 1, reset_peak_memory_stats avant chaque unite, check_physicality(high_quality) sur chaque lot,
reprise (unite deja produite = sautee, et DITE), une erreur n arrete pas les autres unites.
ECARTS DECLARES :
  1. x0 CURES `<PDB>__img0_rep1.pdb` (sha256 relu), pas `6CWG/x0/x0_md0.pdb` ; 28 systemes, pas un.
  2. mode ET step_scale ne sont JAMAIS ecrits ici : le step_scale est LU sur le sampler que BioEmu2 construit pour le mode
     (`sampler_from_preset("edm", mode=mode).step_scale`) et compare a la valeur codee en dur a inference.py:174 (1,5 en
     `structure`, 1,0 sinon) ; un ecart arrete l unite. La propagation utilise `mode="equilibrium"` (inference.py:222).
  3. l embedding est calcule HORS chronometre (`EmuSession._get_embedding`, la meme methode que sample/propagate appellent) :
     son temps et ses pics sont une ligne `session_embedding`, et le debit des unites est un debit de regime permanent.
  4. l OOM est une donnee : ligne `rejet` (batch, pics), puis nouvel essai au batch inferieur (moitie, puis 1).
  5. reservation atomique `mkdir` par systeme (claims/<vague>/<PDB>) ; un worker ne reclame un systeme que si
     temps ecoule + 1,35 x cout predit x max(1, echelle observee) tient dans son budget — sinon il LE DIT et passe au suivant.
  6. sortie ecrite dans `<unite>.tmp` puis renommee : un dossier sans .tmp est complet.
  Le manifeste (une ligne par unite, JSON-lines, append-only) porte tous les champs demandes ; les rejets en sont des entrees.
"""
import argparse, gc, json, os, subprocess, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L

ap = argparse.ArgumentParser()
ap.add_argument("--vague", required=True, help="nom de la vague (dossier de reservation claims/<vague>)")
ap.add_argument("--usages", nargs="+", required=True, choices=["U1", "U2", "U3", "U3D"])
ap.add_argument("--taches", required=True, help="fichier : un PDB par ligne, ordre de reclamation (plus gros d abord)")
ap.add_argument("--batch", type=int, required=True)
ap.add_argument("--batch-json", default="", help="{PDB: batch} qui prime sur --batch")
ap.add_argument("--n", type=int, default=25)
ap.add_argument("--lag-ss", type=float, default=50.0, help="lag UNIQUE du bras disulfure DECLARE (U3D)")
ap.add_argument("--lags", type=float, nargs="+", default=[1, 10, 50, 200, 1000])
ap.add_argument("--budget-s", type=float, default=0.0, help="duree du job - marge ; 0 = pas de garde de budget")
ap.add_argument("--cout-json", default="", help="{PDB: secondes predites pour TOUS les usages de la vague}")
ap.add_argument("--etiquette", required=True)
ap.add_argument("--sortie-racine", default=str(L.G / "sorties"))
ap.add_argument("--a-blanc", action="store_true")
a = ap.parse_args()

RUN = os.environ.get("GEN_RUN", "")
if not a.a_blanc and not RUN:
    sys.exit("ARRET : GEN_RUN non pose (le job doit passer par garde_gen)")
if not os.environ.get("BIOEMU2_MSA_CACHE", "").startswith(str(L.G)):
    sys.exit(f"ARRET : BIOEMU2_MSA_CACHE doit etre sous {L.G}, valeur = {os.environ.get('BIOEMU2_MSA_CACHE')!r}")
if os.environ.get("BIOEMU2_EMBEDS_CACHE"):
    sys.exit("ARRET : BIOEMU2_EMBEDS_CACHE est pose ; il rendrait faux le temps d embedding (garde_gen le defait)")

import mdtraj as md, numpy as np, torch
from bioemu2 import EmuSession, check_physicality, load_model, map_cif_atom_order, save_frames
from bioemu2.sampling import sampler_from_preset
L.interdire_reseau_msa()

JOB = os.environ.get("SLURM_JOB_ID", "hors_slurm")
PRECISION, MATMUL = "bf16-mixed", "high"
SORTIES = Path(a.sortie_racine)
CLAIMS = L.G / "claims" / a.vague
GIO_ = 1.0 / L.GIO
T_JOB0 = time.time()
TACHES = [l.strip() for l in Path(a.taches).read_text().splitlines() if l.strip()]
BATCH_PAR = json.loads(Path(a.batch_json).read_text()) if a.batch_json else {}
COUT = json.loads(Path(a.cout_json).read_text()) if a.cout_json else {}
STEP_ATTENDU = {"structure": 1.5, "equilibrium": 1.0}      # inference.py:174


def step_scale_du_mode(mode):
    s = sampler_from_preset("edm", mode=mode).step_scale
    if s != STEP_ATTENDU[mode]:
        raise L.Rejet(f"step_scale lu = {s} pour le mode {mode}, inference.py:174 en attend {STEP_ATTENDU[mode]}")
    return s


# ── réservation ──────────────────────────────────────────────────────────────
def vivants():
    r = subprocess.run(["squeue", "-h", "-t", "R,PD", "-o", "%i"], capture_output=True, text=True)
    if r.returncode != 0:
        return None                                     # une erreur n est jamais convertie en « aucun vivant »
    return {l.strip() for l in r.stdout.splitlines() if l.strip()}


def reclamer(pdb):
    CLAIMS.mkdir(parents=True, exist_ok=True)
    d = CLAIMS / pdb
    try:
        os.mkdir(d)
    except FileExistsError:
        if (d / "done.json").is_file():
            return False
        try:
            own = json.loads((d / "owner.json").read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            return False                                # reservation en cours d ecriture par un autre : on n y touche pas
        viv = vivants()
        if viv is None or str(own.get("jobid")) in viv or own.get("jobid") == JOB:
            return False
        os.rename(d, CLAIMS / f"{pdb}.perime_{own.get('jobid')}_{int(time.time())}")
        try:
            os.mkdir(d)
        except FileExistsError:
            return False
        print(f"RECLAME_PERIME {pdb} : le job {own.get('jobid')} n est plus vivant et n avait pas fini", flush=True)
    (d / "owner.json").write_text(json.dumps({"jobid": JOB, "t": time.time(), "vague": a.vague}) + "\n")
    return True


# ── préparation CPU d un système ─────────────────────────────────────────────
def preparer(pdb):
    system, seqs, tokens, x0info = L.construire_systeme(pdb)
    msa, prof = L.msa_et_profil(system)
    x0m = map_cif_atom_order(system, x0info["x0"])
    ref = L.charger_ponts(pdb)
    sg = L.indices_sg(x0m.topology, ref["ponts"])
    return dict(system=system, seqs=seqs, tokens=tokens, x0info=x0info, msa=msa, prof=prof, x0m=x0m, ponts=ref["ponts"], sg=sg, ponts_ref=ref)


def deja_produit(d):
    m = d / "meta.json"
    if not (d.is_dir() and m.is_file()):
        return False
    try:
        meta = json.loads(m.read_text())
        return meta.get("statut") == "ok" and L.sha256_fichier(d / "samples.xtc") == meta["sha256_sorties"]["samples.xtc"]
    except (OSError, KeyError, json.JSONDecodeError):
        return False


def n_frames_xtc(p):
    with md.formats.XTCTrajectoryFile(str(p)) as f:
        return len(f)


# ── une unité ────────────────────────────────────────────────────────────────
def produire(P, pdb, usage, regime, mode, step, lag, fn, emu_mode, declares=()):
    unite = f"{pdb}__{usage}__{regime}"
    d = SORTIES / unite
    if deja_produit(d):
        print(f"SAUTE {unite} : deja produit (meta.json + sha256 du xtc relus)", flush=True)
        return "saute"
    seed = L.seed_unite(unite)
    bs0 = int(BATCH_PAR.get(pdb, a.batch))
    echelle = [bs0] + [b for b in (bs0 // 2, 1) if 0 < b < bs0]
    echelle = list(dict.fromkeys(echelle))
    base = {"unite": unite, "systeme": pdb, "usage": usage, "regime": regime, "mode": mode, "step_scale": step,
            "step_scale_origine": "lu sur sampler_from_preset('edm', mode) ; inference.py:174 : 1,5 en structure, 1,0 sinon",
            "lag_ns": lag, "precision": PRECISION, "matmul_precision": MATMUL, "seed": seed, "n": a.n, "tokens": P["tokens"],
            "longueurs": [len(s) for s in P["seqs"]], "msa_profondeur": P["prof"], "x0": P["x0info"]["x0"],
            "sha256_x0": P["x0info"]["sha256_x0"], "x0_atomes": P["x0m"].n_atoms, "jobid": JOB, "etiquette_concurrence": a.etiquette,
            "batch_demande": bs0, **CTX}
    tentatives = []
    for bs in echelle:
        av_t, av_c = time.time(), L.concurrence_instant()["total"]
        try:
            torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
            t = time.time()
            r = fn(bs, seed)
            torch.cuda.synchronize(); dt = time.time() - t
            pk, rs = torch.cuda.max_memory_allocated() * GIO_, torch.cuda.max_memory_reserved() * GIO_
        except torch.cuda.OutOfMemoryError as e:
            ent = {"batch_size": bs, "statut": "oom", "erreur": str(e)[:300],
                   "pic_alloue_Gio": round(torch.cuda.max_memory_allocated() * GIO_, 2), "pic_reserve_Gio": round(torch.cuda.max_memory_reserved() * GIO_, 2)}
            tentatives.append(ent)
            L.ajouter_manifeste({"type": "rejet", **base, "batch_size": bs, "statut": "oom", "motif": "OOM, nouvel essai au batch inferieur",
                                 "erreur": ent["erreur"], "pic_alloue_Gio": ent["pic_alloue_Gio"], "pic_reserve_Gio": ent["pic_reserve_Gio"],
                                 "jobs_concurrents_gpu": {"instant_debut": av_c, **L.concurrence_fenetre(av_t, time.time(), RUN)}})
            print(f"OOM {unite} batch {bs} : alloue {ent['pic_alloue_Gio']} reserve {ent['pic_reserve_Gio']} Gio — batch inferieur", flush=True)
            gc.collect(); torch.cuda.empty_cache()
            continue
        except Exception as e:
            L.ajouter_manifeste({"type": "rejet", **base, "batch_size": bs, "statut": "echec", "erreur": f"{type(e).__name__}: {str(e)[:300]}",
                                 "trace": L.trace_courte()})
            print(f"ECHEC {unite} : {type(e).__name__}: {str(e)[:200]}", flush=True)
            gc.collect(); torch.cuda.empty_cache()
            return "echec"
        # ── l artefact, pas le code de retour ──
        tmp = SORTIES / f"{unite}.tmp"
        if tmp.exists():
            import shutil
            shutil.rmtree(tmp)
        rap = check_physicality(r, "high_quality", system=P["system"])
        ok = int(sum(bool(x) for x in rap.passed))
        save_frames(r, tmp, physicality=rap)
        noms = ["samples.xtc", "topology.cif", "physicality.json"]
        passed = [bool(x) for x in rap.passed]
        raisons = {}
        for x in rap.reasons:
            if x is not None:
                raisons[str(x)] = raisons.get(str(x), 0) + 1
        if P["sg"]:
            dist, resume_p = L.resume_ponts(r, P["sg"], P["ponts"], passed)
            (tmp / "disulfures.json").write_text(json.dumps({"unite": unite, "repere_pont_rompu_A": L.REPERE_PONT_ROMPU_A, "distances_SG_SG_A": dist,
                                                              "retenue_par_check_physicality": passed, "ponts_declares_au_modele": list(declares)}, ensure_ascii=False) + "\n")
            noms.append("disulfures.json")
        else:
            resume_p = {"n_ponts": 0}
        manque = [n for n in noms if not (tmp / n).is_file() or (tmp / n).stat().st_size == 0]
        nf = n_frames_xtc(tmp / "samples.xtc") if "samples.xtc" not in manque else -1
        sha = {n: L.sha256_fichier(tmp / n) for n in noms if n not in manque}
        verdict = "ok" if (not manque and nf == a.n and int(r.n_frames) == a.n and bool(rap.checked)) else "artefact_incomplet"
        ligne = {"type": "unite", **base, "statut": verdict, "ponts_declares_au_modele": list(declares), "disulfures": resume_p,
                 "raisons_check_physicality_des_rejetees": raisons, "batch_size": bs, "tentatives": tentatives, "temps_s": round(dt, 1),
                 "debit_s_par_structure": round(dt / a.n, 2), "embedding_hors_chrono": True, "pic_alloue_Gio": round(pk, 2),
                 "pic_reserve_Gio": round(rs, 2), "retenus": ok, "retenus_sur": a.n, "physicalite_verifiee": bool(rap.checked),
                 "check_physicality": "high_quality", "n_frames_xtc_relu": nf, "fichiers_manquants": manque,
                 "jobs_concurrents_gpu": {"instant_debut": av_c, "instant_fin": L.concurrence_instant()["total"], **L.concurrence_fenetre(av_t, time.time(), RUN)},
                 "rss_max_Gio": L.rss_max_gio(), "elapsed_job_s_a_la_fin_de_l_unite": round(time.time() - T_JOB0, 1),
                 "sortie": str(d), "sha256_sorties": sha}
        if verdict == "ok":
            if d.exists():                                  # reste d une unite incomplete : on la met de cote, on ne l ecrase pas
                os.rename(d, SORTIES / f"{unite}.perime_{int(time.time())}")
            (tmp / "meta.json").write_text(json.dumps({k: v for k, v in ligne.items() if k != "type"} | {"statut": "ok"}, indent=2, ensure_ascii=False) + "\n")
            os.rename(tmp, d)
        else:
            ligne["sortie"] = str(tmp)
        L.ajouter_manifeste(ligne)
        print(f"UNITE {unite} {verdict} : {dt/a.n:.2f} s/struct (batch {bs}, n={a.n}), alloue {pk:.2f}, reserve {rs:.2f} Gio, retenus {ok}/{a.n}, concurrence {ligne['jobs_concurrents_gpu']['min']}-{ligne['jobs_concurrents_gpu']['max']}", flush=True)
        return verdict
    L.ajouter_manifeste({"type": "rejet", **base, "statut": "oom_tous_batches", "motif": f"OOM a tous les batchs {echelle}", "tentatives": tentatives})
    print(f"REJET {unite} : OOM a tous les batchs {echelle}", flush=True)
    return "rejet"


def session_embedding(P, pdb, usage_label, mode):
    emu = EmuSession(P["system"], model=MODEL, mode=mode, msa=P["msa"])
    gc.collect(); torch.cuda.empty_cache(); torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
    av_t, av_c = time.time(), L.concurrence_instant()["total"]
    t = time.time()
    emu._get_embedding()          # methode privee, MAIS c est exactement ce que sample()/propagate() appellent en premier
    torch.cuda.synchronize(); dt = time.time() - t
    L.ajouter_manifeste({"type": "session_embedding", "unite": f"{pdb}__{usage_label}__embedding_{mode}", "systeme": pdb, "mode": mode,
                         "tokens": P["tokens"], "temps_s": round(dt, 1), "pic_alloue_Gio": round(torch.cuda.max_memory_allocated() * GIO_, 2),
                         "pic_reserve_Gio": round(torch.cuda.max_memory_reserved() * GIO_, 2), "precision": PRECISION, "jobid": JOB,
                         "jobs_concurrents_gpu": {"instant_debut": av_c, **L.concurrence_fenetre(av_t, time.time(), RUN)}, "statut": "ok",
                         "etiquette_concurrence": a.etiquette, "msa_profondeur": P["prof"]})
    print(f"EMBEDDING {pdb} {mode} : {dt:.1f} s, alloue {torch.cuda.max_memory_allocated()*GIO_:.2f}, reserve {torch.cuda.max_memory_reserved()*GIO_:.2f} Gio", flush=True)
    return emu


def declarer_ponts(P, pdb):
    """Meme systeme + ponts DECLARES (docs/system.md « Declare every disulfide »), ceux de la MD (PSF), cysteines prises par ORDINAL.
    Derive de 6CWG/scripts/propager.py:systeme_ss, qui les cherchait par geometrie (SG–SG < 2,5 A) ; ici c est la MD qui les declare."""
    from bioemu2 import Atom, Bond, Protein, System
    if not P["ponts"]:
        L.ajouter_manifeste({"type": "rejet", "unite": f"{pdb}__U3D__aucun_pont", "systeme": pdb, "statut": "sans_pont",
                             "motif": "aucun pont dans la MD : rien a declarer, le bras n a pas de sens", "jobid": JOB, "etiquette_concurrence": a.etiquette})
        print(f"SANS_PONT {pdb} : bras U3D sans objet", flush=True)
        return None
    asym = {a_ for e in P["system"].entities for a_ in e.asym_ids}
    lettres = {chr(ord("A") + k) for k in range(len(list(P["x0m"].topology.chains)))}
    if asym != lettres:
        raise L.Rejet(f"asym_ids du systeme {sorted(asym)} != lettres de chaine {sorted(lettres)}")
    bonds = [Bond(atom1=Atom(name="SG", res_name="CYS", res_id=p["ordinal_i"], asym_id=p["chaine_i"]),
                  atom2=Atom(name="SG", res_name="CYS", res_id=p["ordinal_j"], asym_id=p["chaine_j"])) for p in P["ponts"]]
    sys_ss = System([Protein(sequence=e.sequence, asym_ids=list(e.asym_ids)) for e in P["system"].entities], bonds=bonds)
    x0_ss = map_cif_atom_order(sys_ss, P["x0info"]["x0"])
    P2 = dict(P); P2["system"] = sys_ss; P2["x0m"] = x0_ss
    P2["sg"] = L.indices_sg(x0_ss.topology, P["ponts"])
    return P2


def traiter_systeme(pdb):
    bilan = {"ok": 0, "saute": 0, "echec": 0, "rejet": 0, "artefact_incomplet": 0}
    try:
        P = preparer(pdb)
    except (L.Rejet, L.ReseauInterdit) as e:
        L.ajouter_manifeste({"type": "rejet", "unite": f"{pdb}__{'+'.join(a.usages)}__preparation", "systeme": pdb, "statut": "rejet_preparation",
                             "motif": str(e), "jobid": JOB, "etiquette_concurrence": a.etiquette})
        print(f"REJET_PREPARATION {pdb} : {e}", flush=True)
        bilan["rejet"] += 1
        return bilan
    x0m = P["x0m"]
    rep = lambda n: md.Trajectory(xyz=np.repeat(x0m.xyz[:1], n, axis=0), topology=x0m.topology)
    if "U3" in a.usages or "U2" in a.usages:
        emu = session_embedding(P, pdb, "U3" if "U3" in a.usages else "U2", "equilibrium")
        step_eq = step_scale_du_mode("equilibrium")
        if "U3" in a.usages:
            for lag in a.lags:
                r = produire(P, pdb, "U3", f"eq_lag{int(lag):04d}ns", "equilibrium", step_eq, lag,
                             lambda bs, seed, lag=lag: emu.propagate(rep(a.n), lag_time_ns=lag, batch_size=bs, seed=seed, progress=False), "equilibrium")
                bilan[r] = bilan.get(r, 0) + 1
        if "U2" in a.usages:
            r = produire(P, pdb, "U2", "eq", "equilibrium", step_eq, None,
                         lambda bs, seed: emu.sample(num_samples=a.n, batch_size=bs, seed=seed, progress=False), "equilibrium")
            bilan[r] = bilan.get(r, 0) + 1
        del emu; gc.collect(); torch.cuda.empty_cache()
    if "U3D" in a.usages:
        P2 = declarer_ponts(P, pdb)
        if P2 is None:
            return bilan
        emu = session_embedding(P2, pdb, "U3D", "equilibrium")
        step_eq = step_scale_du_mode("equilibrium")
        rep2 = lambda n: md.Trajectory(xyz=np.repeat(P2["x0m"].xyz[:1], n, axis=0), topology=P2["x0m"].topology)
        r = produire(P2, pdb, "U3D", f"eq_lag{int(a.lag_ss):04d}ns_ssdecl", "equilibrium", step_eq, a.lag_ss,
                     lambda bs, seed: emu.propagate(rep2(a.n), lag_time_ns=a.lag_ss, batch_size=bs, seed=seed, progress=False), "equilibrium",
                     declares=[p["id"] for p in P["ponts"]])
        bilan[r] = bilan.get(r, 0) + 1
        del emu; gc.collect(); torch.cuda.empty_cache()
    if "U1" in a.usages:
        emu = session_embedding(P, pdb, "U1", "structure")
        step_st = step_scale_du_mode("structure")
        r = produire(P, pdb, "U1", "struct", "structure", step_st, None,
                     lambda bs, seed: emu.sample(num_samples=a.n, batch_size=bs, seed=seed, progress=False), "structure")
        bilan[r] = bilan.get(r, 0) + 1
        del emu; gc.collect(); torch.cuda.empty_cache()
    return bilan


# ── programme ────────────────────────────────────────────────────────────────
if a.a_blanc:
    for pdb in TACHES:
        try:
            P = preparer(pdb)
            print(f"A BLANC {pdb} : {P['tokens']} tokens, x0 {P['x0m'].n_atoms} atomes, MSA {P['prof']}", flush=True)
        except (L.Rejet, L.ReseauInterdit) as e:
            print(f"A BLANC REJET {pdb} : {e}", flush=True)
    print("step_scale lu :", {m: step_scale_du_mode(m) for m in ("structure", "equilibrium")})
    sys.exit(0)

CTX = L.contexte_job(RUN)
MODEL = load_model(device="cuda", precision=PRECISION, matmul_precision=MATMUL)
print(f"MODELE charge : {torch.cuda.memory_allocated()*GIO_:.2f} Gio, slice {torch.cuda.get_device_properties(0).total_memory*GIO_:.2f} Gio, "
      f"GPU {torch.cuda.get_device_name(0)}, batch {a.batch}, usages {a.usages}, vague {a.vague}", flush=True)
obs, pred, resume = 0.0, 0.0, {}


def msa_prete(pdb):
    """Le job MSA serie a ecrit msa/<PDB>/profondeur.json avec statut ok (l entree, pas le code de retour du job)."""
    p = L.G / "msa" / pdb / "profondeur.json"
    try:
        return json.loads(p.read_text()).get("statut") == "ok"
    except (OSError, json.JSONDecodeError):
        return False


def msa_serie_vivante():
    """Vrai si un processus tient le verrou du job MSA serie (flock non bloquant ; on le relache aussitot)."""
    import fcntl
    p = L.G / "msa" / ".verrou_serie"
    if not p.is_file():
        return False
    with open(p, "a+") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(f, fcntl.LOCK_UN)
        return False


passes = 0
while True:
    passes += 1
    progres, attente = 0, []
    for pdb in TACHES:
        if (CLAIMS / pdb / "done.json").is_file():
            continue
        if not msa_prete(pdb):
            attente.append(pdb)
            continue
        p_s = float(COUT.get(pdb, 0.0))
        ech = max(1.0, obs / pred) if pred > 0 else 1.0
        ecoule = time.time() - T_JOB0
        if a.budget_s > 0 and p_s > 0 and ecoule + 1.35 * ech * p_s > a.budget_s:
            print(f"NON_RECLAME {pdb} : {ecoule:.0f} s ecoules + 1,35 x {ech:.2f} x {p_s:.0f} s predits > budget {a.budget_s:.0f} s — laisse a un autre worker / a la reprise", flush=True)
            continue
        if not reclamer(pdb):
            continue
        t0 = time.time()
        print(f"RECLAME {pdb} (job {JOB}) · {time.time() - T_JOB0:.0f} s ecoules", flush=True)
        b = traiter_systeme(pdb)
        dur = time.time() - t0
        resume[pdb] = {"duree_s": round(dur, 1), "predit_s": p_s, **b}
        if p_s > 0:
            obs += dur; pred += p_s
        (CLAIMS / pdb / "done.json").write_text(json.dumps({"jobid": JOB, "duree_s": round(dur, 1), "bilan": b, "t": time.time()}) + "\n")
        print(f"TERMINE {pdb} en {dur:.0f} s (predit {p_s:.0f} s) : {b}", flush=True)
        progres += 1
    if attente and progres == 0 and msa_serie_vivante() and (a.budget_s <= 0 or time.time() - T_JOB0 < a.budget_s - 600):
        print(f"EN_ATTENTE_MSA passe {passes} : {attente} — le job MSA serie tourne encore, nouvelle passe dans 90 s", flush=True)
        time.sleep(90)
        continue
    if progres == 0:
        if attente:
            print(f"SANS_MSA {attente} : pas de profondeur.json ok et plus de job MSA serie vivant — non traites", flush=True)
        break
(L.G / "runs" / f"resume_{a.vague}_job{JOB}.json").write_text(json.dumps(resume, indent=2, ensure_ascii=False) + "\n")
tot = {k: sum(v.get(k, 0) for v in resume.values()) for k in ("ok", "saute", "echec", "rejet", "artefact_incomplet")}
print(f"FIN generer : {len(resume)} systeme(s) traites par ce worker : {tot} (rejets = {tot['echec'] + tot['rejet'] + tot['artefact_incomplet']})", flush=True)
sys.exit(0)
