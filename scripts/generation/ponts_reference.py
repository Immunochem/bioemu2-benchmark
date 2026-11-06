#!/usr/bin/env python3
"""ponts_reference.py — les ponts disulfure de chacun des 28 systemes et leur distance SG–SG DANS LA MD (la reference).

DEFINITION D UN PONT : une liaison SG–SG de la liste `!NBOND` du `topology.psf` de la MD (c est la MD qui les DECLARE).
Controle croise : la meme paire, lue sur l IMAGE 0 du xtc, doit valoir la distance de l IMAGE 0 relue dans le x0 CURE —
c est ce qui prouve que la correspondance PSF -> x0 est la bonne (🔴 ORDINALE, jamais par resSeq : 5HGG, 4M3K, 3RJQ portent des
codes d insertion que mdtraj ecrase).
CORRESPONDANCE : segment PSF `PRO<c>` -> chaine d origine <c> -> chaine du x0 (table_chaines du JSON compagnon) ; le residu est
le n-ieme du segment (n = ordinal 1..L) ; le JSON compagnon (`map_residus`) doit donner la meme resSeq et un nom CYS/CYX.
SORTIE : generation/ponts/reference.json — par systeme : ponts (chaine+ordinal des deux cysteines, resSeq d origine), distance de
reference par replique (n images, image 0, mediane, min, p5, p95, max), mediane des medianes et plage des medianes par replique,
compte de cysteines libres, et le compte de systemes SANS aucun pont (ils ne comptent ni pour ni contre : le denominateur).
Un echec de lecture est une ENTREE (statut=ecarte + raison), jamais un saut.
"""
import json, os, sys, time, traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L
import mdtraj as md, numpy as np

import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument("--sortie", default="reference.json")
_a = _ap.parse_args()
C = L.BENCH / "corpus/systemes"
OUT = L.G / "ponts"
OUT.mkdir(exist_ok=True)


def lire_psf(chemin):
    """(atomes, liaisons) ; atomes = [(segid, resid_str, resname, name)], liaisons = [(i, j)] indices 0-based."""
    atomes, liaisons = [], []
    mode, attendu = None, 0
    with open(chemin) as f:
        for ligne in f:
            if "!NATOM" in ligne:
                mode, attendu = "atom", int(ligne.split()[0]); continue
            if "!NBOND" in ligne:
                mode, attendu = "bond", int(ligne.split()[0]); continue
            if "!" in ligne and mode:
                mode = None; continue
            if mode == "atom" and ligne.strip():
                t = ligne.split()
                atomes.append((t[1], t[2], t[3], t[4]))
            elif mode == "bond" and ligne.strip():
                v = [int(x) for x in ligne.split()]
                liaisons += [(v[k] - 1, v[k + 1] - 1) for k in range(0, len(v), 2)]
    if len(atomes) == 0:
        raise ValueError("aucun atome lu dans le PSF")
    return atomes, liaisons


def residus_par_segment(atomes):
    """{segid: [(resid_str, resname)] dans l ordre d apparition} + atome -> (segid, rang 1-based)."""
    seg, rang_atome, vu = {}, [], {}
    for (s, r, n, a) in atomes:
        cle = (s, r)
        if cle not in vu:
            vu[cle] = len(seg.setdefault(s, [])) + 1
            seg[s].append((r, n))
        rang_atome.append((s, vu[cle]))
    return seg, rang_atome


def systeme(pdb):
    rec = {"pdb": pdb, "statut": "ok"}
    j1 = json.loads((L.CUR / f"{pdb}__img0_rep1.json").read_text())
    table = j1["table_chaines"]
    mapr = j1["map_residus"]
    par_chaine = {}
    for m in mapr:
        par_chaine.setdefault(m["chaine"], {})[m["ordinal"]] = m
    longueurs = {c: len(v) for c, v in par_chaine.items()}
    x0 = md.load(str(L.chemin_x0(pdb, 1)))
    x0_ch = {chr(ord("A") + i): ch for i, ch in enumerate(x0.topology.chains)}
    ponts_ref = None
    rec["repliques"] = {}
    for rep in (1, 2, 3):
        d = C / pdb / f"rep{rep}"
        atomes, liaisons = lire_psf(d / "topology.psf")
        seg, rang_atome = residus_par_segment(atomes)
        prot = {s: v for s, v in seg.items() if s.startswith("PRO")}
        # segment PRO<c> -> chaine du x0
        corr = {}
        for s in prot:
            c_orig = s[3:]
            if c_orig not in table:
                raise ValueError(f"segment {s} : chaine d origine {c_orig!r} absente de table_chaines {table}")
            corr[s] = table[c_orig]
        for s, c in corr.items():
            if len(prot[s]) != longueurs.get(c):
                raise ValueError(f"segment {s} -> chaine {c} : {len(prot[s])} residus dans le PSF, {longueurs.get(c)} dans le x0")
        sg = [i for i, (s, r, n, a) in enumerate(atomes) if a == "SG" and s in prot]
        ponts = []
        sgset = set(sg)
        for (i, j) in liaisons:
            if i in sgset and j in sgset:
                (si, ri), (sj, rj) = rang_atome[i], rang_atome[j]
                a, b = sorted([(corr[si], ri, i), (corr[sj], rj, j)])
                ponts.append({"chaine_i": a[0], "ordinal_i": a[1], "chaine_j": b[0], "ordinal_j": b[1], "atome_md_i": a[2], "atome_md_j": b[2]})
        ponts.sort(key=lambda p: (p["chaine_i"], p["ordinal_i"], p["chaine_j"], p["ordinal_j"]))
        # coherence avec le JSON compagnon : resSeq d origine + nom CYS/CYX
        for p in ponts:
            for k in ("i", "j"):
                m = par_chaine[p[f"chaine_{k}"]][p[f"ordinal_{k}"]]
                if m["nom_orig"] not in ("CYS", "CYX") or m["nom"] != "CYS":
                    raise ValueError(f"pont {p}: le residu ({p[f'chaine_{k}']},{p[f'ordinal_{k}']}) est {m['nom_orig']}/{m['nom']} dans map_residus")
                p[f"resSeq_orig_{k}"] = m["resSeq_orig"]
                resid_psf = atomes[p[f"atome_md_{k}"]][1]
                if str(m["resSeq_orig"]) != resid_psf:
                    p[f"alerte_resid_psf_{k}"] = f"PSF {resid_psf} != map_residus {m['resSeq_orig']} (code d insertion ? l ordinal fait foi)"
        cys_tot = sum(1 for m in mapr if m["nom"] == "CYS")
        info = {"n_ponts": len(ponts), "n_cys_total": cys_tot, "n_cys_libres": cys_tot - 2 * len(ponts), "ponts": ponts}
        if ponts_ref is None:
            ponts_ref = ponts
        elif [(p["chaine_i"], p["ordinal_i"], p["chaine_j"], p["ordinal_j"]) for p in ponts] != [(p["chaine_i"], p["ordinal_i"], p["chaine_j"], p["ordinal_j"]) for p in ponts_ref]:
            raise ValueError(f"rep{rep} : ponts differents de rep1")
        # distances dans la MD : seulement les SG des ponts
        stats = {}
        brutes = {}
        if ponts:
            idx = sorted({p["atome_md_i"] for p in ponts} | {p["atome_md_j"] for p in ponts})
            top = md.load(str(d / "structure.pdb")).topology
            for i in idx:
                at = top.atom(i)
                if at.name != "SG" or at.residue.name not in ("CYS", "CYX"):
                    raise ValueError(f"atome MD {i} = {at.name} {at.residue.name} dans structure.pdb, SG attendu (ordre PSF != ordre PDB ?)")
            t = md.load(str(d / "trajectory.xtc"), top=str(d / "structure.pdb"), atom_indices=idx)
            pos = {a: k for k, a in enumerate(idx)}
            for p in ponts:
                dist = np.linalg.norm(t.xyz[:, pos[p["atome_md_i"]]] - t.xyz[:, pos[p["atome_md_j"]]], axis=1) * 10.0
                pid = f"{p['chaine_i']}{p['ordinal_i']}-{p['chaine_j']}{p['ordinal_j']}"
                brutes[pid] = dist
                stats[pid] = {"n_images": int(len(dist)), "image0_A": round(float(dist[0]), 3), "mediane_A": round(float(np.median(dist)), 3),
                              "min_A": round(float(dist.min()), 3), "p5_A": round(float(np.percentile(dist, 5)), 3),
                              "p95_A": round(float(np.percentile(dist, 95)), 3), "max_A": round(float(dist.max()), 3)}
        rec["repliques"][f"rep{rep}"] = {**{k: v for k, v in info.items() if k != "ponts"}, "distances_md": stats}
        rec.setdefault("_brutes", {})[rep] = brutes
    # pont par pont : agregat sur les 3 repliques + controle x0 / image 0
    rec["ponts"] = []
    rec["n_ponts"] = len(ponts_ref)
    rec["n_cys_libres"] = rec["repliques"]["rep1"]["n_cys_libres"]
    for p in ponts_ref:
        pid = f"{p['chaine_i']}{p['ordinal_i']}-{p['chaine_j']}{p['ordinal_j']}"
        meds = [rec["repliques"][f"rep{r}"]["distances_md"][pid]["mediane_A"] for r in (1, 2, 3)]
        # x0 curé : SG–SG sur l image 0 de la replique 1
        def sg(chaine, ordinal):
            res = list(x0_ch[chaine].residues)[ordinal - 1]
            return [a.index for a in res.atoms if a.name == "SG"][0], res.name
        (ia, na), (ib, nb) = sg(p["chaine_i"], p["ordinal_i"]), sg(p["chaine_j"], p["ordinal_j"])
        d_x0 = float(np.linalg.norm(x0.xyz[0, ia] - x0.xyz[0, ib]) * 10.0)
        d_md0 = rec["repliques"]["rep1"]["distances_md"][pid]["image0_A"]
        rec["ponts"].append({
            "id": pid, **{k: v for k, v in p.items() if not k.startswith("atome_md")},
            "residus_x0": [na, nb], "x0_img0_rep1_A": round(d_x0, 3), "md_img0_rep1_A": d_md0,
            "controle_x0_vs_md_image0_A": round(abs(d_x0 - d_md0), 4),
            "mediane_des_medianes_par_replique_A": round(float(np.median(meds)), 3),
            "plage_des_medianes_par_replique_A": [round(min(meds), 3), round(max(meds), 3)],
            "plage_toutes_images_A": [min(rec["repliques"][f"rep{r}"]["distances_md"][pid]["min_A"] for r in (1, 2, 3)),
                                      max(rec["repliques"][f"rep{r}"]["distances_md"][pid]["max_A"] for r in (1, 2, 3))],
        })
        if abs(d_x0 - d_md0) > 0.01:
            rec["statut"] = "ecarte"
            rec["raison"] = f"pont {pid} : SG–SG du x0 ({d_x0:.3f} A) != image 0 de la MD ({d_md0:.3f} A) : la correspondance ordinale est fausse"
    # quantiles poolés (3 repliques) par pont, et NUL par systeme : fraction des images de la MD qui ont >= 1 pont hors de [p0,5 ; p99,5] poolé
    brutes = rec.pop("_brutes")
    if rec["ponts"]:
        bornes, out_any = {}, None
        for pt in rec["ponts"]:
            pid = pt["id"]
            pool = np.concatenate([brutes[r][pid] for r in (1, 2, 3)])
            lo, hi = float(np.percentile(pool, 0.5)), float(np.percentile(pool, 99.5))
            pt["md_pool_p0_5_A"], pt["md_pool_p99_5_A"] = round(lo, 3), round(hi, 3)
            pt["md_pool_n_images"] = int(len(pool))
            bornes[pid] = (lo, hi)
            o = (pool < lo) | (pool > hi)
            out_any = o if out_any is None else (out_any | o)
        rec["nul_md_fraction_images_avec_un_pont_hors_p0_5_p99_5"] = round(float(out_any.mean()), 4)
        rec["nul_md_n_images"] = int(len(out_any))
        rec["nul_md_note"] = ("fraction des images de la MD (3 repliques poolees) dont >= 1 pont sort de [p0,5 ; p99,5] POOLE de son pont : "
                              "c est le taux que ce critere donne a la MD elle-meme, a lire comme le zero du critere (~ 1-0,99^k pour k ponts independants)")
    return rec


tout, rejets = {}, {}
t0 = time.time()
for pdb in sorted(L.corpus_28()):
    try:
        r = systeme(pdb)
    except Exception as e:
        r = {"pdb": pdb, "statut": "ecarte", "raison": f"{type(e).__name__}: {str(e)[:400]}", "trace": traceback.format_exc()[-800:]}
    tout[pdb] = r
    if r["statut"] != "ok":
        rejets[pdb] = r["raison"]
    print(f"PONTS {pdb} {r['statut']} n_ponts={r.get('n_ponts')} libres={r.get('n_cys_libres')} "
          f"{[(p['id'], p['md_img0_rep1_A'], p['mediane_des_medianes_par_replique_A'], p['controle_x0_vs_md_image0_A']) for p in r.get('ponts', [])]} "
          f"{r.get('raison', '')}", flush=True)
sans = sorted(p for p, r in tout.items() if r["statut"] == "ok" and r["n_ponts"] == 0)
bilan = {"n_systemes": len(tout), "n_ok": sum(1 for r in tout.values() if r["statut"] == "ok"), "n_rejets": len(rejets), "rejets": rejets,
         "systemes_sans_pont": sans, "n_sans_pont": len(sans),
         "n_avec_pont": sum(1 for r in tout.values() if r["statut"] == "ok" and r["n_ponts"] > 0),
         "n_ponts_total": sum(r.get("n_ponts", 0) for r in tout.values() if r["statut"] == "ok"), "duree_s": round(time.time() - t0, 1)}
(OUT / _a.sortie).write_text(json.dumps({"bilan": bilan, "systemes": tout}, indent=1, ensure_ascii=False) + "\n")
print("PONTS_BILAN", json.dumps(bilan, ensure_ascii=False), f"(rejets = {len(rejets)})", flush=True)
sys.exit(0 if not rejets else 1)
