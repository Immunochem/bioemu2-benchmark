#!/usr/bin/env python3
"""redock_gpu.py — REDOCKING BioEmu2 : pour chaque tache (systeme cristallin, ou paire NEGATIVE), construit les leurres (redock_cond.construire_leurres : memes graines que la verification CPU),
les propage (EmuSession.propagate, mode equilibrium, sampler par defaut, step_scale 1,0 laisse au harnais) et mesure contre le cristal. UNE EmuSession par systeme (embeddings calcules une fois).
  - une propagation par lag, toutes les images du systeme dans UN appel (x0 differents d une image a l autre : sampler edm par defaut, images independantes) ;
  - echelle d OOM batch 4 -> 2 -> 1 ; un OOM a tous les batchs est un REJET ecrit, pas un trou ;
  - MSA : CACHE SEUL (requete ColabFold interdite ici) ; cache absent = rejet `msa_absent` ecrit ;
  - sorties par tache/lag : samples.xtc + topology.cif (save_frames), x0_leurres.xtc, conditions.json (audit des leurres), mesures.json ; meta.json (statut, debit, concurrence).
APPARIEMENT ORDINAL : la reference est mappee par map_cif_atom_order dans l ordre du System, les echantillons sortent dans le meme ordre ; les noms d atomes sont compares avant toute mesure."""
import argparse, gc, json, os, sys, time, traceback
from pathlib import Path

sys.dont_write_bytecode = True
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
if os.environ.get("BIOEMU2_MSA_CACHE", "") != str(R / "msa_cache"):
    sys.exit(f"ARRET : BIOEMU2_MSA_CACHE doit valoir {R / 'msa_cache'}")
if not os.environ.get("REDOCK_RUN"):
    sys.exit("ARRET : REDOCK_RUN non pose (le job doit passer par garde_redock_gpu)")
sys.path.insert(0, str(R / "scripts"))

ap = argparse.ArgumentParser()
ap.add_argument("--taches", required=True)
ap.add_argument("--ids", nargs="*", default=[])
ap.add_argument("--lags", type=float, nargs="+", default=[50.0])
ap.add_argument("--batch", type=int, default=4)
ap.add_argument("--K", type=int, default=3)
ap.add_argument("--sortie", required=True)
ap.add_argument("--budget-s", type=float, default=0.0, help="ne lance une tache que si 1,35 x l estimation mesuree tient dans ce budget (0 = sans limite)")
ap.add_argument("--s-par-struct-prior", type=float, default=14.0, help="prior (pilote 6CWG : 12,4-13 s/structure a 390 residus, batch 4) avant toute mesure du job")
a = ap.parse_args()

import mdtraj as md, numpy as np, torch
import bioemu2.msa.colabfold as cf


def _interdit(*x, **k):
    raise RuntimeError("requete ColabFold interdite hors du job MSA serie (docs/msa.md : query serially)")
cf.run_colabfold_msa = _interdit
from bioemu2 import EmuSession, get_msa, load_model, save_frames
import redock_lib as RL, redock_cond as RC, redock_bioemu as RB

GIO = 2 ** 30
JOB = os.environ.get("SLURM_JOB_ID", "hors_slurm")
RUN = Path(os.environ["REDOCK_RUN"])
SORTIE = Path(a.sortie); SORTIE.mkdir(parents=True, exist_ok=True)
T_JOB0 = time.time()


def concurrence_fenetre(t0, t1, intervalle=15.0):
    f = RUN / "concurrence.jsonl"
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
    return {"min": min(vals) if vals else None, "max": max(vals) if vals else None, "n_echantillons": len(vals), "n_erreurs_echantillon": err}


taches = json.loads(Path(a.taches).read_text())
if a.ids:
    taches = [t for t in taches if t["id"] in set(a.ids)]
print(f"{len(taches)} tache(s) · lags {a.lags} · K={a.K} · batch {a.batch}", flush=True)
model = load_model(device="cuda", precision="bf16-mixed", matmul_precision="high")
GPU = torch.cuda.get_device_name(0)
mesures_debit = []          # (tokens, n, s_par_struct)
bilan, rejets = [], []

for t in taches:
    tid, kind = t["id"], t["kind"]
    d = SORTIE / tid
    meta_f = d / "meta.json"
    if meta_f.is_file():
        try:
            if json.loads(meta_f.read_text()).get("statut") == "ok":
                print(f"SAUTE {tid} : deja produit", flush=True); continue
        except (OSError, json.JSONDecodeError):
            pass
    try:
        info = RB.systeme_depuis_x0(t["x0"])
        Rf = RL.Reference(info["ref"])
        conds = [tuple(c) for c in t["conds"]] if t.get("conds") else (RC.COND_POS if kind == "pos" else RC.COND_NEG)   # `conds` : surcharge explicite (tache speciale)
        frames, audit = RC.construire_leurres(Rf, tid, conds, a.K)
        n = len(frames)
        refus_leurre = [x for x in audit if not x["accepte"]]
        if n == 0:
            raise RuntimeError("aucun leurre accepte")
        # budget : estimation mesuree (la plus lente rencontree, mise a l echelle ~ tokens^2) sinon prior
        if a.budget_s:
            if mesures_debit:
                est = max(s * (info["tokens"] / tk) ** 2 for tk, nn, s in mesures_debit) * n * len(a.lags)
            else:
                est = a.s_par_struct_prior * (info["tokens"] / 390.0) ** 2 * n * len(a.lags)
            if time.time() - T_JOB0 + 1.35 * est > a.budget_s:
                rejets.append({"id": tid, "motif": "non_lance_budget", "estimation_s": round(est)}); print(f"NON LANCE {tid} : 1,35 x {est:.0f} s ne tient pas dans le budget", flush=True); continue
        msa = get_msa(info["nu"])
    except Exception as e:
        motif = "msa_absent" if "interdite" in str(e) else "preparation"
        rejets.append({"id": tid, "motif": motif, "erreur": f"{type(e).__name__}: {str(e)[:300]}"})
        print(f"REJET {tid} ({motif}) : {type(e).__name__}: {str(e)[:200]}", flush=True)
        continue
    x0traj = RC.traj_de(Rf, frames)
    emu = EmuSession(info["systeme"], model=model, mode="equilibrium", msa=msa)
    d.mkdir(parents=True, exist_ok=True)
    x0traj.save_xtc(str(d / "x0_leurres.xtc"))
    (d / "conditions.json").write_text(json.dumps({"id": tid, "kind": kind, "strate": t["strate"], "K": a.K, "audit": audit, "refus_leurres": len(refus_leurre)}) + "\n")
    meta = {"id": tid, "strate": t["strate"], "kind": kind, "x0": t["x0"], "tokens": info["tokens"], "ponts_declares": info["ponts"], "n_ponts_inter_chaine": info["n_ponts_inter_chaine"],
            "d_nat_A": round(Rf.d_nat, 3), "n_contacts_natifs": Rf.n_contacts_nat, "n_images": n, "n_leurres_refuses": len(refus_leurre), "mode": "equilibrium", "step_scale": "defaut du harnais (1,0)",
            "precision": "bf16-mixed", "matmul_precision": "high", "gpu": GPU, "slurm_job": JOB, "lags": {}, "statut": None}
    ok_tous = True
    for lag in a.lags:
        seed = RL.seed_de(tid, "propagate", int(lag))
        echelle = list(dict.fromkeys([a.batch] + [b for b in (a.batch // 2, 1) if 0 < b < a.batch]))
        r, dt, bs_ok, tent = None, None, None, []
        for bs in echelle:
            t_av, c_av = time.time(), concurrence_fenetre(time.time(), time.time())
            try:
                torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats(); t0 = time.time()
                r = emu.propagate(x0traj, lag_time_ns=lag, batch_size=bs, seed=seed)
                torch.cuda.synchronize(); dt = time.time() - t0; bs_ok = bs
                break
            except torch.cuda.OutOfMemoryError as e:
                tent.append({"batch": bs, "statut": "oom", "pic_alloue_Gio": round(torch.cuda.max_memory_allocated() / GIO, 2), "pic_reserve_Gio": round(torch.cuda.max_memory_reserved() / GIO, 2)})
                print(f"OOM {tid} lag {lag} batch {bs}", flush=True); gc.collect(); torch.cuda.empty_cache()
            except Exception as e:
                tent.append({"batch": bs, "statut": "echec", "erreur": f"{type(e).__name__}: {str(e)[:300]}"}); print(f"ECHEC {tid} lag {lag}: {type(e).__name__}: {str(e)[:200]}", flush=True)
                gc.collect(); torch.cuda.empty_cache(); break
        if r is None:
            rejets.append({"id": tid, "lag": lag, "motif": "propagation", "tentatives": tent}); ok_tous = False; continue
        ld = d / f"lag{int(lag):04d}"
        save_frames(r, ld)
        tr = md.Trajectory(xyz=r.xyz, topology=Rf.ref.topology)
        try:
            mes = Rf.mesurer(tr) if kind == "pos" else RC.mesurer_neg(Rf, tr, t["site_chaine"], t["site_residus"])
        except Exception as e:
            mes = {"ERREUR_MESURE": f"{type(e).__name__}: {str(e)[:300]}"}; ok_tous = False
        with md.formats.XTCTrajectoryFile(str(ld / "samples.xtc")) as f_:
            nf = len(f_)
        manque = [x for x in ("samples.xtc", "topology.cif") if not (ld / x).is_file() or (ld / x).stat().st_size == 0]
        (ld / "mesures.json").write_text(json.dumps({"id": tid, "lag_ns": lag, "indices_image": {x["cond"] + f"_r{x['rep']}": x.get("indice_image") for x in audit if x["accepte"]}, "mesures": mes}) + "\n")
        fen = concurrence_fenetre(t_av, time.time())
        meta["lags"][str(int(lag))] = {"seed": seed, "batch_size": bs_ok, "tentatives": tent, "temps_s": round(dt, 1), "s_par_structure": round(dt / n, 2), "embedding_dans_le_chrono": True,
                                       "pic_alloue_Gio": round(torch.cuda.max_memory_allocated() / GIO, 2), "pic_reserve_Gio": round(torch.cuda.max_memory_reserved() / GIO, 2),
                                       "n_images_xtc": nf, "fichiers_manquants": manque, "jobs_concurrents_gpu": {"debut": c_av, **fen}}
        if nf != n or manque:
            ok_tous = False
        mesures_debit.append((info["tokens"], n, dt / n))
        print(f"UNITE {tid} ({t['strate']}, {kind}) lag {int(lag)} : {dt / n:.2f} s/struct (batch {bs_ok}, n={n}, tokens {info['tokens']}) · pic alloue {meta['lags'][str(int(lag))]['pic_alloue_Gio']} Gio · "
              f"concurrence {fen['min']}-{fen['max']} · refus de leurres {len(refus_leurre)}", flush=True)
        del r, tr
        gc.collect(); torch.cuda.empty_cache()
    meta["statut"] = "ok" if ok_tous else "artefact_incomplet"
    meta_f.write_text(json.dumps(meta, indent=1, ensure_ascii=False) + "\n")
    bilan.append({"id": tid, "statut": meta["statut"]})
    del emu
    gc.collect(); torch.cuda.empty_cache()

(SORTIE / f"bilan_{JOB}.json").write_text(json.dumps({"bilan": bilan, "rejets": rejets, "duree_s": round(time.time() - T_JOB0, 1)}, indent=1) + "\n")
print(f"FIN : {sum(b['statut'] == 'ok' for b in bilan)} tache(s) ok, {sum(b['statut'] != 'ok' for b in bilan)} incomplete(s), {len(rejets)} rejet(s) (REJETS = {len(rejets)}) · {time.time() - T_JOB0:.0f} s", flush=True)
for rj in rejets:
    print("REJET", json.dumps(rj, ensure_ascii=False)[:300], flush=True)
sys.exit(0 if not rejets and all(b["statut"] == "ok" for b in bilan) else 1)
