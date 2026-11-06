#!/usr/bin/env python3
"""bioemu2_t4l.py — BioEmu2 sur T4L : des UNITES (equilibrium / propagate / steer), les CV de Bhakat calculees EN LIGNE, reprise possible.

Une unite = un dossier `<racine>/<nom>/` : brut/ (samples.xtc + topology.cif, `save_frames`), meta.json, cv.npz (d1, d2 en nm ; p ; physique),
physicalite.json. Une unite deja produite est SAUTEE et DITE. Une erreur n'arrete pas les autres unites mais fait sortir 1.

Systemes :
  s164 : la sequence du PAPIER (164 residus, WT, C54 et C97 presentes)
  s162 : la sequence MODELISEE de 5LZM (ATOM 1-162 ; Asn163 et Leu164 absents du cristal, REMARK 465) = l'objet de la MD non biaisee de l'auteur
         (UnbiasedMD/prot-masses-noh.pdb : 162 residus). Le x0 des propagations est le cristal 5LZM, proteine seule, passe par curer_entree.py.
Mode : equilibrium partout, sampler "edm" -> step_scale 1,0 calcule PAR LE HARNAIS depuis le mode (sampling.py:54 ; inference.py:174). On le LIT et on le consigne.
"""
import argparse, json, os, sys, time, traceback
from collections import Counter
from pathlib import Path

sys.dont_write_bytecode = True
T = Path("/gpfs/phherman/bioemu2_bench/t4l")
sys.path.insert(0, str(T / "scripts"))
if not os.environ.get("BIOEMU2_MSA_CACHE", "").startswith("/gpfs/"):
    sys.exit("ARRET : BIOEMU2_MSA_CACHE doit pointer sous /gpfs (jamais ~).")

ap = argparse.ArgumentParser()
ap.add_argument("--plan", required=True, help="JSON : liste d'unites")
ap.add_argument("--racine", required=True, help="dossier des unites (sous t4l/)")
ap.add_argument("--msa-seulement", action="store_true", help="recuperer/relever les MSA puis sortir (aucun GPU)")
ap.add_argument("--budget-s", type=float, default=0.0, help="si > 0 : une unite dont 1,35 x n x s/struct depasse le budget restant est NON RECLAMEE, dite")
a = ap.parse_args()

import mdtraj as md, numpy as np
import cv_t4l as C
if not a.msa_seulement:
    import torch
from bioemu2.sampling import sampler_from_preset
from bioemu2 import (EDMSampler, EmuSession, GuidedFKSampler, Potential, System, check_physicality, get_msa, load_model, map_cif_atom_order,
                     save_frames)

RACINE = T / a.racine
PLAN = json.loads(Path(a.plan).read_text())
X0 = T / "x0/x0_5lzm.pdb"
T0 = time.time()

# ── systemes ───────────────────────────────────────────────────────────────────
x0_ref = md.load(str(X0))
seq162 = "".join(C.UN[r.name] for r in x0_ref.topology.residues)
assert seq162 == C.SEQ_PAPIER[:162], "la sequence modelisee de 5LZM n'est plus les 162 premiers residus du papier"
assert x0_ref.n_residues == 162
SYS = {"s164": System.from_description(proteins=[C.SEQ_PAPIER]), "s162": System.from_description(proteins=[seq162])}
MSA, PROF = {}, {}
for k, s in SYS.items():
    t = time.time()
    MSA[k] = get_msa(s)                      # UNE requete a la fois, serie ; cache sous BIOEMU2_MSA_CACHE
    a3m = MSA[k].a3m(MSA[k].chain_ids[0])
    PROF[k] = {"longueur": len(s.entities[0].sequence), "lignes_en_tete": sum(1 for l in a3m.splitlines() if l.startswith(">")), "temps_s": round(time.time() - t, 1)}
    print(f"MSA {k} : {PROF[k]}", flush=True)
(T / "msa").mkdir(exist_ok=True)
(T / "msa/profondeur.json").write_text(json.dumps({"cache": os.environ["BIOEMU2_MSA_CACHE"], "systemes": PROF}, indent=2) + "\n")
if a.msa_seulement:
    sys.exit(0)
X0M = map_cif_atom_order(SYS["s162"], str(X0))
print(f"x0 5LZM mappe : {X0M.n_atoms} atomes, {X0M.n_residues} residus", flush=True)

MODELES, SESSIONS = {}, {}


def emu_de(sk, prec):
    if prec not in MODELES:
        MODELES[prec] = load_model(device="cuda", precision=prec, matmul_precision="high" if prec == "bf16-mixed" else None)
    if (sk, prec) not in SESSIONS:
        SESSIONS[(sk, prec)] = EmuSession(SYS[sk], model=MODELES[prec], mode="equilibrium", msa=MSA[sk])
    return SESSIONS[(sk, prec)]


class FermerD1(Potential):
    """Mur UNILATERAL sur d1 : E = 0,5 * kappa * relu(d1 - T)^2, d1 = distance Ca44-Ca150 en A, atomes pris PAR ORDINAL (le 44e et le 150e CA du systeme).
    Patron : steering/scripts/steer_q2.py:89 (mode gd) et upstream/examples/tilt_a_sample_distribution_with_steering.py. kappa = 1*1/(1+1) = 0,5 : pour deux atomes
    seuls, la meme normalisation que la leur (na*nb/(na+nb)) fait combler une fraction `guidance_weight` de l'ecart par pas de guidage."""

    def __init__(self, target_angstrom, guidance_weight):
        self.T, self.gw, self.kappa, self.ia, self.ib = float(target_angstrom), float(guidance_weight), 0.5, None, None

    def prepare(self, system):
        ca = [x.index for x in system.topology.atoms if x.name == "CA"]
        names = [x.residue.name for x in system.topology.atoms if x.name == "CA"]
        assert names[43] == "SER" and names[149] == "ILE", (names[43], names[149])
        self.ia, self.ib = ca[43], ca[149]

    def compute_parameters(self, t):
        return {"guidance_weight": self.gw, "guidance_interval": 1, "resampling_weight": 1.0}

    def __call__(self, coords_angstrom, *, t):
        d = torch.linalg.vector_norm(coords_angstrom[:, self.ia] - coords_angstrom[:, self.ib], dim=-1)
        return 0.5 * self.kappa * torch.relu(d - self.T).pow(2)


def decrire_sampler(s):
    d = {"classe": type(s).__name__}
    for k in ("num_steps", "step_scale", "noise_scale", "resampling_interval", "num_gd_steps", "ess_threshold", "weight_update", "reheat"):
        if hasattr(s, k):
            d[k] = getattr(s, k)
    if hasattr(s, "potentials"):
        d["potentiels"] = [type(p).__name__ for p in s.potentials]
    return d


def faire(u, idx):
    nom = u["nom"]; d = RACINE / nom
    if (d / "cv.npz").is_file() and (d / "meta.json").is_file():
        print(f"SAUTE {nom} : deja produit", flush=True); return "saute"
    n, bs, prec, kind, sk = int(u["n"]), int(u["batch"]), u.get("precision", "bf16-mixed"), u["type"], u["systeme"]
    s_est = float(u.get("s_est", 0))
    if a.budget_s > 0 and s_est > 0 and (time.time() - T0) + 1.35 * n * s_est > a.budget_s:
        print(f"NON_RECLAME {nom} : {time.time()-T0:.0f} s ecoules + 1,35 x {n} x {s_est} s/struct > budget {a.budget_s:.0f} s", flush=True); return "non_reclame"
    emu = emu_de(sk, prec)
    seed = int(u["seed"])
    sampler_desc = None
    try:
        torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats(); t = time.time()
        if kind == "equilibrium":
            r = emu.sample(num_samples=n, batch_size=bs, seed=seed, progress=False)
            sampler_desc = decrire_sampler(sampler_from_preset("edm", mode=emu.mode))   # CE QUE LE HARNAIS CONSTRUIT (sampling.py:54), pas une copie ecrite a la main
        elif kind == "propagate":
            assert sk == "s162", "le x0 5LZM n'existe que pour s162"
            x0 = md.Trajectory(xyz=np.repeat(X0M.xyz[:1], n, axis=0), topology=X0M.topology)
            r = emu.propagate(x0, lag_time_ns=float(u["lag"]), batch_size=bs, seed=seed, progress=False)
            sampler_desc = decrire_sampler(sampler_from_preset("edm", mode="equilibrium"))   # inference.py : propagate() construit ce preset
        elif kind == "steer":
            from bioemu2.denoising.steering.potentials import get_default_potentials
            pot = FermerD1(float(u["cible_A"]), float(u["f"]))
            s = GuidedFKSampler(resampling_interval=None, potentials=[pot] + (list(get_default_potentials()) if u.get("phys", True) else []))
            r = emu.sample(num_samples=n, batch_size=bs, sampler=s, seed=seed, progress=False)
            sampler_desc = decrire_sampler(s)
        else:
            raise ValueError(kind)
        torch.cuda.synchronize(); dt = time.time() - t
        pic_a, pic_r = torch.cuda.max_memory_allocated() / 2**30, torch.cuda.max_memory_reserved() / 2**30
        save_frames(r, d / "brut")
        rap = check_physicality(r, "high_quality", system=SYS[sk])
        passed = np.array([bool(p) for p in rap.passed])
        cvr = C.cv(r)
        if not (1.0 < float(np.median(cvr["d1"])) < 7.0):
            raise ValueError(f"d1 median {np.median(cvr['d1']):.3f} hors de 1-7 nm : unite suspecte (Å pris pour nm ?)")
        np.savez(d / "cv.npz", d1=cvr["d1"], d2=cvr["d2"], p=cvr["p"], physique=passed)
        (d / "physicalite.json").write_text(json.dumps({"n": n, "passed": int(passed.sum()), "rejetees": int((~passed).sum()),
                                                        "raisons": Counter(str(x) for x, ok in zip(rap.reasons, passed) if not ok)}, indent=1, ensure_ascii=False) + "\n")
        meta = {"unite": nom, "type": kind, "systeme": sk, "n": int(r.n_frames), "lag_ns": u.get("lag"), "mode": "equilibrium",
                "sampler": sampler_desc or "defaut de EmuSession.propagate (edm, step_scale 1,0 : sampler_from_preset(mode='equilibrium'))",
                "precision": prec, "batch_size": bs, "seed": seed, "msa": PROF[sk], "temps_s": round(dt, 1), "s_par_structure": round(dt / n, 3),
                "pic_alloue_Gio": round(pic_a, 2), "pic_reserve_Gio": round(pic_r, 2), "slurm_job": os.environ.get("SLURM_JOB_ID"),
                "gpu": torch.cuda.get_device_name(0), "spec": u, "n_residus": int(r.n_residues), "rapport_indices": cvr["rapport"]}
        (d / "meta.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False, default=str) + "\n")
        print("UNITE", json.dumps({k: meta[k] for k in ("unite", "n", "temps_s", "s_par_structure", "pic_alloue_Gio", "pic_reserve_Gio")}, ensure_ascii=False),
              f"| physiques {int(passed.sum())}/{n} | d1 med {np.median(cvr['d1']):.3f} nm | d1<2,5 : {int((cvr['d1'] < 2.5).sum())}", flush=True)
        return "ok"
    except Exception:
        print("ERREUR", nom, traceback.format_exc(), flush=True)
        return "echec"


bilan = Counter()
for i, u in enumerate(PLAN):
    bilan[faire(u, i)] += 1
(RACINE).mkdir(parents=True, exist_ok=True)
(RACINE / "bilan.json").write_text(json.dumps({"bilan": dict(bilan), "duree_s": round(time.time() - T0, 1)}, indent=1) + "\n")
print(f"FIN : {dict(bilan)} · {time.time()-T0:.0f} s", flush=True)
sys.exit(1 if bilan["echec"] else 0)
