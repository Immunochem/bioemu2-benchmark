"""ponts_analyse.py — agrege, depuis le manifeste, les ponts disulfure par systeme / pont / usage / lag, et ce que check_physicality en dit."""
import json
from collections import defaultdict

VHH = {"5E5M", "4POU", "5SV3", "4Y7M", "5HGG", "4M3K", "6CWG", "6DBG", "3RJQ", "5JMO", "5VNW"}


def unites(L):
    """Une entree par unite (la derniere ok), usages U1/U2/U3/U3D."""
    u = {}
    for d in L:
        if d.get("type") == "unite" and d.get("statut") == "ok" and d.get("usage") in ("U1", "U2", "U3", "U3D"):
            u[d["unite"]] = d
    return u


def par_pont(u):
    """{(systeme, pont): {usage_regime: dict}}"""
    out = defaultdict(dict)
    for k, d in u.items():
        ds = d.get("disulfures") or {}
        for pid, v in (ds.get("par_pont") or {}).items():
            out[(d["systeme"], pid)][d["usage"] + ":" + d["regime"]] = {**v, "n": d["n"], "ligne": d["_ligne"], "retenus": d["retenus"]}
    return out


def par_systeme_regime(u):
    """Par (systeme, usage, regime) : n structures, n avec un pont rompu, n retenues, n rompues ET retenues."""
    rows = []
    for k, d in sorted(u.items()):
        ds = d.get("disulfures") or {}
        rows.append({"systeme": d["systeme"], "usage": d["usage"], "regime": d["regime"], "lag_ns": d.get("lag_ns"), "n": d["n"],
                     "n_ponts": ds.get("n_ponts", 0), "n_struct_pont_rompu": ds.get("n_structures_avec_un_pont_rompu"),
                     "n_retenues": d["retenus"], "n_rompues_et_retenues": ds.get("n_structures_rompues_et_retenues_par_check_physicality"),
                     "ligne": d["_ligne"], "vhh": d["systeme"] in VHH,
                     "raisons": d.get("raisons_check_physicality_des_rejetees", {})})
    return rows


def totaux(rows, filtre):
    sel = [r for r in rows if filtre(r) and r["n_ponts"] > 0]
    n = sum(r["n"] for r in sel)
    rompu = sum(r["n_struct_pont_rompu"] or 0 for r in sel)
    ret = sum(r["n_retenues"] for r in sel)
    rr = sum(r["n_rompues_et_retenues"] or 0 for r in sel)
    return {"n_unites": len(sel), "n_systemes": len({r["systeme"] for r in sel}), "n_structures": n, "n_avec_pont_rompu": rompu,
            "n_retenues_par_check_physicality": ret, "n_rompues_et_retenues": rr}


def charger_disulfures(ssh):
    """Tous les disulfures.json de sorties/ en UN appel ssh : {unite: json}."""
    prog = ("import glob,json,os\n"
            "o={}\n"
            "for p in glob.glob('/gpfs/phherman/bioemu2_bench/generation/sorties/*/disulfures.json'):\n"
            "    o[os.path.basename(os.path.dirname(p))]=json.load(open(p))\n"
            "print(json.dumps(o))\n")
    import chargerman
    return json.loads(chargerman.ssh_retry("python3 -", entree=prog))


def par_structure(dis, ref, unite, pdb):
    """Pour chaque structure d une unite : (rompue_3A, hors_plage_MD, retenue). Plage MD = [p0,5 ; p99,5] des images poolees des 3 repliques (reference_v2.json)."""
    d = dis[unite]
    ponts = {p["id"]: p for p in ref["systemes"][pdb]["ponts"]}
    n = len(d["retenue_par_check_physicality"])
    rompue, hors = [False] * n, [False] * n
    for pid, dist in d["distances_SG_SG_A"].items():
        lo, hi = ponts[pid]["md_pool_p0_5_A"], ponts[pid]["md_pool_p99_5_A"]
        for k, x in enumerate(dist):
            if x > d["repere_pont_rompu_A"]:
                rompue[k] = True
            if x < lo or x > hi:
                hors[k] = True
    return rompue, hors, d["retenue_par_check_physicality"]
