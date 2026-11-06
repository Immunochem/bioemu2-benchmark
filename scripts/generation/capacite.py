#!/usr/bin/env python3
"""capacite.py — courbe de capacite de BioEmu2 : VRAM (allouee ET reservee), debit, echelle de batch, par systeme.

DERIVE de `vram_1030/vram_1zm4.py` (que l on ne modifie pas ; il ecrit sous vram_1030/). On reprend SA mesure :
chrono() = synchronize + reset_peak_memory_stats + propagate(lag 50 ns, equilibrium, bf16-mixed, matmul high) +
max_memory_allocated ET max_memory_reserved ; premier appel (embedding) a part ; refus des residus non standard.

ECARTS DECLARES (chacun adosse a une mesure ou a une consigne) :
  1. x0 = x0 CURE `<PDB>__img0_rep1.pdb` (sha256 relu contre son JSON compagnon), pas `_input.pdb` brut : le contrat d entree
     (RAPPORT_CURATION §10) est deja applique par curer_entree.py (noms d atomes, chaines A,B…). Meme objet que l ancien chemin :
     2AJF -> 6319 atomes aux deux endroits (vram_2AJF.json:`x0_atomes` = 6319 = `n_atomes_x0` du JSON compagnon).
  2. `empty_cache()` AVANT reset_peak_memory_stats() a l echauffement de chaque batch. Raison, mesuree : dans vram_2AJF.json le
     reserve du « regime permanent » (22,08) egale celui du premier appel (22,07) alors que l alloue passe de 14,99 a 8,09 :
     le reserve du regime etait le cache HERITE de l embedding, pas le besoin du debruitage. On releve donc les DEUX :
     `pic_reserve_apres_vidage` (besoin du batch) et `pic_reserve_continu` (cache herite, regime de production).
  3. echelle de batch 1, 2, 4, 5, 8 (consigne du coordinateur ; 5 parce que 25 = 5 x 5), un echauffement JETE par batch,
     puis un chronometrage n = max(4, batch). Un OOM est une entree, pas un echec : on s arrete au premier OOM de l echelle.
  4. check_physicality sur chaque lot chronometre (high_quality), compte de retenus.
  5. chaque essai est une ligne du manifeste, avec jobs_concurrents_gpu releve (instant debut/fin + fenetre du sampler).
"""
import argparse, gc, json, os, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L

ap = argparse.ArgumentParser()
ap.add_argument("pdbs", nargs="+")
ap.add_argument("--batches", type=int, nargs="+", default=[1, 2, 4, 5, 8])
ap.add_argument("--n-chrono", type=int, default=4)
ap.add_argument("--lag-ns", type=float, default=50.0)
ap.add_argument("--etiquette", required=True, help="solo | c4 | … : le regime de concurrence PREVU, ecrit dans l unite")
ap.add_argument("--barriere-n", type=int, default=0)
ap.add_argument("--barriere-dir", default="")
ap.add_argument("--a-blanc", action="store_true", help="CPU seulement : systeme, MSA depuis le cache, x0 ; aucun modele")
a = ap.parse_args()
RUN = os.environ.get("GEN_RUN", "")
if not a.a_blanc and not RUN:
    sys.exit("ARRET : GEN_RUN non pose (le job doit passer par garde_gen)")
if not os.environ.get("BIOEMU2_MSA_CACHE", "").startswith(str(L.G)):
    sys.exit(f"ARRET : BIOEMU2_MSA_CACHE doit etre sous {L.G} (jamais msa/ ni 6CWG/), valeur = {os.environ.get('BIOEMU2_MSA_CACHE')!r}")

import mdtraj as md, numpy as np, torch
from bioemu2 import EmuSession, check_physicality, load_model, map_cif_atom_order
L.interdire_reseau_msa()

CIBLE = L.G / "capacite"
CIBLE.mkdir(exist_ok=True)
PRECISION = "bf16-mixed"
JOB = os.environ.get("SLURM_JOB_ID", "hors_slurm")
resultats = {"job": JOB, "etiquette": a.etiquette, "lag_ns": a.lag_ns, "precision": PRECISION, "systemes": {}}


def entree_manifeste(pdb, regime, base, **kw):
    e = {"type": "essai_capacite", "unite": f"{pdb}__CAP__{regime}", "systeme": pdb, "usage": "CAP", "regime": regime,
         "mode": "equilibrium", "step_scale": 1.0, "step_scale_origine": "defaut de sampler_from_preset(mode) — inference.py:174",
         "lag_ns": a.lag_ns, "precision": PRECISION, "matmul_precision": "high", "etiquette_concurrence": a.etiquette}
    e.update(base); e.update(kw)
    L.ajouter_manifeste(e)


def barriere():
    if a.barriere_n <= 0:
        return None
    d = Path(a.barriere_dir)
    d.mkdir(parents=True, exist_ok=True)
    (d / JOB).write_text(str(time.time()))
    t0 = time.time()
    while len(list(d.iterdir())) < a.barriere_n and time.time() - t0 < 900:
        time.sleep(2)
    n = len(list(d.iterdir()))
    r = {"atteinte": n >= a.barriere_n, "presents": n, "attendus": a.barriere_n, "attente_s": round(time.time() - t0, 1)}
    print("BARRIERE", json.dumps(r), flush=True)
    return r


# ── préparation de chaque système (CPU) ──────────────────────────────────────
prep = {}
for pdb in a.pdbs:
    try:
        system, seqs, tokens, x0info = L.construire_systeme(pdb)
        t0 = time.time()
        msa, prof = L.msa_et_profil(system)
        x0m = map_cif_atom_order(system, x0info["x0"])
        prep[pdb] = dict(system=system, seqs=seqs, tokens=tokens, x0info=x0info, msa=msa, prof=prof, x0m=x0m,
                         t_msa=round(time.time() - t0, 1))
        print(f"PREP {pdb} : {tokens} tokens, chaines {[len(s) for s in seqs]}, x0 {x0m.n_atoms} atomes, MSA {prof} ({prep[pdb]['t_msa']} s)", flush=True)
    except (L.Rejet, L.ReseauInterdit) as e:
        print(f"REJET {pdb} : {e}", flush=True)
        resultats["systemes"][pdb] = {"statut": "rejet_preparation", "raison": str(e)}
        if not a.a_blanc:
            L.ajouter_manifeste({"type": "rejet", "unite": f"{pdb}__CAP__preparation", "systeme": pdb, "usage": "CAP",
                                 "statut": "rejet", "raison": str(e), "jobid": JOB, "etiquette_concurrence": a.etiquette})
if a.a_blanc:
    (CIBLE / f"a_blanc_{JOB}.json").write_text(json.dumps(
        {p: {"tokens": v["tokens"], "x0_atomes": v["x0m"].n_atoms, "msa": v["prof"]} for p, v in prep.items()}
        | {"rejets": {p: r for p, r in resultats["systemes"].items()}}, indent=2, ensure_ascii=False) + "\n")
    print(f"A BLANC : {len(prep)} systeme(s) prets, {len(resultats['systemes'])} rejet(s)")
    sys.exit(0)

ctx = L.contexte_job(RUN)
G_ = 1.0 / L.GIO
model = load_model(device="cuda", precision=PRECISION, matmul_precision="high")
vram_slice = round(torch.cuda.get_device_properties(0).total_memory * G_, 2)
mem_modele = round(torch.cuda.memory_allocated() * G_, 2)
print(f"MODELE charge : {mem_modele} Gio alloues, slice {vram_slice} Gio, GPU {torch.cuda.get_device_name(0)}", flush=True)
resultats.update(vram_slice_Gio=vram_slice, modele_alloue_Gio=mem_modele, gpu=torch.cuda.get_device_name(0), contexte=ctx)
resultats["barriere"] = barriere()


def essai(emu, x0m, n, bs, seed, vider):
    """(trajectoire, dt_s, pic_alloue, pic_reserve, reserve_au_depart). Meme discipline que vram_1zm4.py:65-70."""
    torch.cuda.synchronize()
    if vider:
        gc.collect(); torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    res0 = torch.cuda.memory_reserved() * G_
    t = time.time()
    r = emu.propagate(md.Trajectory(xyz=np.repeat(x0m.xyz[:1], n, axis=0), topology=x0m.topology),
                      lag_time_ns=a.lag_ns, batch_size=bs, seed=seed, progress=False)
    torch.cuda.synchronize()
    return r, time.time() - t, torch.cuda.max_memory_allocated() * G_, torch.cuda.max_memory_reserved() * G_, res0


def avant():
    return time.time(), L.concurrence_instant()["total"]


def apres(av):
    """jobs_concurrents_gpu pour [debut de l essai, maintenant] : instant au debut, instant a la fin, fenetre du sampler."""
    t0, c0 = av
    return {"instant_debut": c0, "instant_fin": L.concurrence_instant()["total"], **L.concurrence_fenetre(t0, time.time(), RUN)}


for pdb, p in prep.items():
    res = {"tokens": p["tokens"], "longueurs": [len(s) for s in p["seqs"]], "x0_atomes": p["x0m"].n_atoms, "msa": p["prof"],
           "temps_get_msa_s": p["t_msa"], "essais": []}
    resultats["systemes"][pdb] = res
    base = {"systeme": pdb, "tokens": p["tokens"], "msa_profondeur": p["prof"], "x0": p["x0info"]["x0"],
            "sha256_x0": p["x0info"]["sha256_x0"], "jobid": JOB, "x0_atomes": p["x0m"].n_atoms, **ctx}
    emu = EmuSession(p["system"], model=model, mode="equilibrium", msa=p["msa"])
    try:
        # ── premier appel, JETE pour le debit : il contient l embedding (le pic qui a trompe la 1re courbe) ──
        av = avant()
        _, dt, pk, rs, rs0 = essai(emu, p["x0m"], 1, 1, 1, True)
        res["premier_appel"] = {"temps_s": round(dt, 1), "pic_alloue_Gio": round(pk, 2), "pic_reserve_Gio": round(rs, 2),
                                "note": "embedding inclus, JETE pour le debit"}
        print(f"{pdb} premier appel (embedding) : {dt:.1f} s, alloue {pk:.2f}, reserve {rs:.2f} Gio", flush=True)
        entree_manifeste(pdb, f"eq_lag{int(a.lag_ns):04d}ns_premier_appel_embedding", base, statut="ok_jete_pour_le_debit",
                         batch_size=1, n=1, seed=1, temps_s=round(dt, 1), pic_alloue_Gio=round(pk, 2), pic_reserve_Gio=round(rs, 2),
                         jobs_concurrents_gpu=apres(av), rss_max_Gio=L.rss_max_gio())
    except torch.cuda.OutOfMemoryError as e:
        res["premier_appel"] = {"statut": "oom", "erreur": str(e)[:300], "pic_alloue_avant_oom_Gio": round(torch.cuda.max_memory_allocated() * G_, 2),
                                "pic_reserve_avant_oom_Gio": round(torch.cuda.max_memory_reserved() * G_, 2)}
        print(f"{pdb} OOM AU PREMIER APPEL (embedding) : {res['premier_appel']}", flush=True)
        entree_manifeste(pdb, f"eq_lag{int(a.lag_ns):04d}ns_premier_appel_embedding", base, statut="oom", batch_size=1, n=1,
                         **{k: v for k, v in res["premier_appel"].items() if k != "statut"}, jobs_concurrents_gpu=apres(av))
        del emu; gc.collect(); torch.cuda.empty_cache()
        continue
    # ── échelle de batch ──
    for bs in a.batches:
        regime = f"eq_lag{int(a.lag_ns):04d}ns_bs{bs}"
        n = max(a.n_chrono, bs)
        ess = {"batch_size": bs, "n_chrono": n}
        try:
            av = avant()
            _, dtw, pkw, rsw, rsw0 = essai(emu, p["x0m"], bs, bs, 1000 + bs, True)      # echauffement JETE, apres vidage
            ess["echauffement_jete"] = {"temps_s": round(dtw, 1), "s_par_structure": round(dtw / bs, 2),
                                        "pic_alloue_Gio": round(pkw, 2), "pic_reserve_apres_vidage_Gio": round(rsw, 2)}
            r, dt, pk, rs, rs0 = essai(emu, p["x0m"], n, bs, 2000 + bs, False)           # chronometrage, cache herite
            rap = check_physicality(r, "high_quality", system=p["system"])
            ok = int(sum(bool(x) for x in rap.passed))
            ess.update(statut="ok", temps_s=round(dt, 1), s_par_structure=round(dt / n, 2), pic_alloue_Gio=round(pk, 2),
                       pic_reserve_continu_Gio=round(rs, 2), reserve_au_depart_Gio=round(rs0, 2),
                       pic_reserve_apres_vidage_Gio=round(rsw, 2), physiques=ok, physicalite_verifiee=bool(rap.checked), n_frames=int(r.n_frames))
            print(f"{pdb} bs{bs} : {dt/n:.2f} s/struct (n={n}), alloue {pk:.2f}, reserve continu {rs:.2f}, apres vidage {rsw:.2f}, "
                  f"retenus {ok}/{n}, echauffement {dtw/bs:.2f} s/struct", flush=True)
            entree_manifeste(pdb, regime, base, statut="ok", batch_size=bs, n=n, seed=2000 + bs, temps_s=ess["temps_s"],
                             debit_s_par_structure=ess["s_par_structure"], pic_alloue_Gio=ess["pic_alloue_Gio"],
                             pic_reserve_Gio=ess["pic_reserve_continu_Gio"], pic_reserve_apres_vidage_Gio=ess["pic_reserve_apres_vidage_Gio"],
                             echauffement_jete=ess["echauffement_jete"], retenus=ok, retenus_sur=n, physicalite_verifiee=ess["physicalite_verifiee"],
                             jobs_concurrents_gpu=apres(av), rss_max_Gio=L.rss_max_gio(), sorties=None,
                             note_sorties="trajectoires de capacite NON conservees : seule la mesure compte")
            res["essais"].append(ess)
        except torch.cuda.OutOfMemoryError as e:
            ess.update(statut="oom", erreur=str(e)[:300], pic_alloue_avant_oom_Gio=round(torch.cuda.max_memory_allocated() * G_, 2),
                       pic_reserve_avant_oom_Gio=round(torch.cuda.max_memory_reserved() * G_, 2))
            print(f"{pdb} bs{bs} OOM : alloue {ess['pic_alloue_avant_oom_Gio']} / reserve {ess['pic_reserve_avant_oom_Gio']} Gio — on s arrete a ce batch", flush=True)
            entree_manifeste(pdb, regime, base, statut="oom", batch_size=bs, n=n, erreur=ess["erreur"],
                             pic_alloue_Gio=ess["pic_alloue_avant_oom_Gio"], pic_reserve_Gio=ess["pic_reserve_avant_oom_Gio"],
                             jobs_concurrents_gpu=apres(av))
            res["essais"].append(ess)
            gc.collect(); torch.cuda.empty_cache()
            break
        except Exception as e:
            ess.update(statut="echec", erreur=f"{type(e).__name__}: {str(e)[:300]}", trace=L.trace_courte())
            print(f"{pdb} bs{bs} ECHEC {ess['erreur']}", flush=True)
            entree_manifeste(pdb, regime, base, statut="echec", batch_size=bs, n=n, erreur=ess["erreur"], trace=ess["trace"],
                             jobs_concurrents_gpu=apres(av))
            res["essais"].append(ess)
            gc.collect(); torch.cuda.empty_cache()
            break
    del emu
    gc.collect(); torch.cuda.empty_cache()
    (CIBLE / f"{pdb}__{a.etiquette}__job{JOB}.json").write_text(json.dumps({**{k: v for k, v in resultats.items() if k != "systemes"}, "systeme": pdb, **res},
                                                                           indent=2, ensure_ascii=False) + "\n")

(CIBLE / f"job{JOB}__{a.etiquette}.json").write_text(json.dumps(resultats, indent=2, ensure_ascii=False) + "\n")
n_rej = sum(1 for v in resultats["systemes"].values() if v.get("statut") == "rejet_preparation")
print(f"FIN capacite : {len(prep)} systeme(s) mesure(s), {n_rej} rejet(s) de preparation (rejets = {n_rej})", flush=True)
sys.exit(0)
