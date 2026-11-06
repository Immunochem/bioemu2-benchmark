#!/usr/bin/env python3
"""plancher.py — CPU. Pour chaque systeme PRET des strates : (1) le PLANCHER DU HASARD (sans modele) ; (2) la VERIFICATION des leurres exactement tels que le job GPU les construira.
  F1 = placement aleatoire (rotation uniforme + direction uniforme) a d_COM = native, sans choc  ·  F2 = placement aleatoire amene au contact (distance atomique min 3,5 A)
  memes metriques que pour le modele (fnat, iRMSD, lRMSD, DockQ, classe CAPRI) contre le meme cristal.
Aucun rejet silencieux : leurres refuses (choc impossible apres N essais) et leurres dont la verification echoue sont comptes par condition."""
import json, multiprocessing as mp, os, sys, time, traceback
from pathlib import Path
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "1"
sys.dont_write_bytecode = True
R = Path("/gpfs/phherman/bioemu2_bench/redocking")
sys.path.insert(0, str(R / "scripts"))
import numpy as np

N_PLANCHER = 600


def systemes_prets(strates):
    out = []
    for nom in strates:
        m = json.loads((R / f"strate_{nom}/manifeste.json").read_text())
        S = m.get("systemes", m)
        for pid, e in sorted(S.items()):
            if e["etat"] == "PRET":
                out.append((nom, pid, str(R / f"strate_{nom}/{pid}/x0.pdb")))
    return out


def travailler(arg):
    nom, pid, x0 = arg
    t0 = time.time()
    try:
        import redock_lib as RL, redock_cond as RC, redock_bioemu as RB
        info = RB.systeme_depuis_x0(x0)
        Rf = RL.Reference(info["ref"])
        pl = RL.plancher(Rf, N_PLANCHER, N_PLANCHER, RL.seed_de(pid, "plancher"))
        frames, audit = RC.construire_leurres(Rf, pid, RC.COND_POS)
        rej = [a for a in audit if not a["accepte"]]
        res = {"strate": nom, "id": pid, "x0": x0, "tokens": info["tokens"], "ponts": info["ponts"], "n_ponts_inter_chaine": info["n_ponts_inter_chaine"],
               "d_nat_A": round(Rf.d_nat, 3), "n_contacts_natifs": Rf.n_contacts_nat, "taille_chaines": [len(r) for r in Rf.res], "chaine_mobile": Rf.lig,
               "plancher": pl, "leurres_audit": audit, "leurres_refuses": len(rej), "duree_s": round(time.time() - t0, 1)}
        d = R / "resultats/plancher"
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{nom}__{pid}.json").write_text(json.dumps(res) + "\n")
        row = {"strate": nom, "id": pid, "tokens": info["tokens"], "d_nat": res["d_nat_A"], "n_nat": Rf.n_contacts_nat, "refuses": len(rej), "n_leurres": len(audit),
               "inter_chaine": info["n_ponts_inter_chaine"]}
        for k in ("F1", "F2"):
            m = pl[k]
            if m is None:
                row[k] = None; continue
            f, dq = np.array(m["fnat"]), np.array(m["dockq"])
            cl = np.array(m["classe_capri"])
            row[k] = {"n": len(f), "fnat_p50": float(np.median(f)), "fnat_p99": float(np.percentile(f, 99)), "fnat_max": float(f.max()), "dockq_p99": float(np.percentile(dq, 99)),
                      "dockq_max": float(dq.max()), "frac_acceptable_ou_mieux": float((cl != "incorrect").mean()), "frac_medium_ou_mieux": float(np.isin(cl, ["medium", "high"]).mean())}
        return row
    except Exception as e:
        return {"strate": nom, "id": pid, "ERREUR": f"{type(e).__name__}: {str(e)[:300]}", "trace": traceback.format_exc()[-800:]}


if __name__ == "__main__":
    strates = sys.argv[1:] or ["A_corpus", "A_apparie", "B"]
    taches = systemes_prets(strates)
    print(f"{len(taches)} systemes PRETS dans {strates}", flush=True)
    n = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    with mp.Pool(n) as p:
        rows = []
        for r in p.imap_unordered(travailler, taches):
            rows.append(r)
            if "ERREUR" in r:
                print("ERREUR", r["strate"], r["id"], r["ERREUR"], flush=True); continue
            print(f"{r['strate']:9s} {r['id']} tok={r['tokens']:3d} d_nat={r['d_nat']:6.1f} n_nat={r['n_nat']:3d} leurres_refuses={r['refuses']}/{r['n_leurres']} ponts_inter={r['inter_chaine']} | "
                  + " | ".join(f"{k}: fnat p50={r[k]['fnat_p50']:.3f} p99={r[k]['fnat_p99']:.3f} max={r[k]['fnat_max']:.3f} DockQ max={r[k]['dockq_max']:.3f} ≥acceptable={r[k]['frac_acceptable_ou_mieux']:.3f}"
                                for k in ("F1", "F2") if r[k]), flush=True)
    (R / "resultats/plancher_resume.json").write_text(json.dumps(rows) + "\n")
    err = [r for r in rows if "ERREUR" in r]
    ok = [r for r in rows if "ERREUR" not in r]
    print(f"BILAN_PLANCHER : {len(ok)} systemes mesures, {len(err)} en erreur (REJETS = {len(err)}) ; leurres refuses au total {sum(r['refuses'] for r in ok)}/{sum(r['n_leurres'] for r in ok)}", flush=True)
    sys.exit(1 if err else 0)
