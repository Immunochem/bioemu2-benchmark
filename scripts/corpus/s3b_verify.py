"""J3b — verification PAR LE CONTENU de chaque replique telechargee (pas par la taille).

Par replique : images, atomes, chaines proteiques, residus proteiques, residus non proteiques (eau/ions),
axe `time` (c'est un numero d'image), ruptures Calpha-Calpha > 4,5 A sur les 1000 images, distance entre
les deux chaines, paires de residus a 5,0 A (atomes lourds) a l'image 0 puis toutes les 50 images.
Par systeme : les trois repliques different-elles entre elles (ecart brut en nm) ; et pour les systemes deja
sur staging, `production1.dcd` est-il la replique 1 (ecart ~0,0005 nm) et non la 2 ou la 3 (0,7-2,8 nm) ?
Une anomalie est une DONNEE ecrite dans le resultat.
"""
import json
import os
import sys
import time

import mdtraj as md
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import (CACHE, CONTACT_A, ROOT, STAGING_INPUTS, jdump, jload, now_iso, sha256_file)  # noqa: E402

STRIDE_PAIRES = 50


def paires(xyz_A, idx0, idx1, res_of):
    t0, t1 = cKDTree(xyz_A[idx0]), cKDTree(xyz_A[idx1])
    s = set()
    for i, lst in enumerate(t0.query_ball_tree(t1, CONTACT_A + 1e-6)):
        for j in lst:
            s.add((res_of[idx0[i]], res_of[idx1[j]]))
    return len(s)


def verifier_replique(d, api):
    out = {"dossier": d, "date": now_iso(), "anomalies": []}
    top_pdb = d + "/structure.pdb"
    pdb = md.load(top_pdb)
    t = md.load(d + "/trajectory.xtc", top=top_pdb)
    top = t.topology
    out["n_images"] = t.n_frames
    out["n_atomes"] = t.n_atoms
    out["n_atomes_pdb"] = pdb.n_atoms
    out["time_premier"] = float(t.time[0])
    out["time_dernier"] = float(t.time[-1])
    out["time_pas"] = float(t.time[1] - t.time[0]) if t.n_frames > 1 else None
    out["boite_nm_image0"] = [round(float(x), 3) for x in t.unitcell_lengths[0]] if t.unitcell_lengths is not None else None
    chains = list(top.chains)
    prot = [c for c in chains if sum(1 for r in c.residues if r.is_protein) >= 0.9 * len(list(c.residues))]
    nonprot = sum(1 for r in top.residues if not r.is_protein)
    water = sum(1 for r in top.residues if r.is_water)
    out.update({"n_chaines": len(chains), "n_chaines_proteiques": len(prot),
                "n_residus_proteine": sum(1 for r in top.residues if r.is_protein),
                "n_residus_non_proteine": nonprot, "n_eau": water,
                "residus_par_chaine": [len(list(c.residues)) for c in prot]})
    # ----- anomalies par rapport a l'API
    if t.n_frames != api["mdFrames"]:
        out["anomalies"].append("images %d != API %s" % (t.n_frames, api["mdFrames"]))
    if len(prot) != 2:
        out["anomalies"].append("chaines proteiques %d != 2" % len(prot))
    if out["n_residus_proteine"] != api["PROTRES"]:
        out["anomalies"].append("residus proteiques %d != API PROTRES %s" % (out["n_residus_proteine"], api["PROTRES"]))
    if out["n_atomes"] != api["mdAtoms"]:
        out["anomalies"].append("atomes %d != API mdAtoms %s" % (out["n_atomes"], api["mdAtoms"]))
    if nonprot:
        out["anomalies"].append("%d residus non proteiques (eau/ions/ligands)" % nonprot)
    if t.n_atoms != pdb.n_atoms:
        out["anomalies"].append("atomes xtc %d != structure.pdb %d" % (t.n_atoms, pdb.n_atoms))
    # ----- ecart image 0 du xtc / structure.pdb (meme replique : ~0,0005 nm ; 1B6C rep2 = anomalie connue L0a)
    if t.n_atoms == pdb.n_atoms:
        out["ecart_max_xtc0_vs_structurepdb_nm"] = round(float(np.abs(t.xyz[0] - pdb.xyz[0]).max()), 4)
    # ----- ruptures CA-CA, distance entre chaines
    elem = [(a.element.symbol if a.element is not None else None) for a in top.atoms]
    heavy = np.array([e not in ("H", "D", None) for e in elem])
    res_of = np.array([a.residue.index for a in top.atoms])
    rupt = []
    cms = []
    for c in prot:
        ca = [a.index for a in c.atoms if a.name == "CA" and a.residue.is_protein]
        dd = np.linalg.norm(np.diff(t.xyz[:, ca, :], axis=1), axis=2) * 10.0   # A, (images, n-1)
        rupt.append({"dmax_A": round(float(dd.max()), 3), "images_avec_rupture_gt_4.5A": int((dd.max(axis=1) > 4.5).sum())})
        cms.append(t.xyz[:, ca, :].mean(axis=1))
    out["ca_ruptures_par_chaine"] = rupt
    if len(cms) == 2:
        dc = np.linalg.norm(cms[0] - cms[1], axis=1) * 10.0
        out["distance_centres_chaines_A"] = {"min": round(float(dc.min()), 2), "max": round(float(dc.max()), 2),
                                             "image0": round(float(dc[0]), 2), "derniere": round(float(dc[-1]), 2)}
        idx = [[a.index for a in c.atoms if heavy[a.index]] for c in prot]
        frames = list(range(0, t.n_frames, STRIDE_PAIRES))
        pp = [paires(t.xyz[f] * 10.0, idx[0], idx[1], res_of) for f in frames]
        out["paires_5A"] = {"image0": pp[0], "min": int(min(pp)), "mediane": float(np.median(pp)), "max": int(max(pp)),
                            "images_testees": len(frames), "stride": STRIDE_PAIRES}
    return out, t


def main():
    sel = jload(ROOT + os.environ.get("CORPUS_SELECTION", "/selection.json"))["selection"]
    inv = {}
    projs = {q["accession"]: q for q in jload("/gpfs/phherman/bioemu2_bench/L0a_work/all_projects.json")}
    SUF = os.environ.get("CORPUS_SUFFIX", "")
    rep_out = open(CACHE + "/verif_replicas%s.jsonl" % SUF, "w", encoding="utf-8")
    systemes = {}
    t0 = time.time()
    for s in sel:
        pdb, acc = s["pdb"], s["accession"]
        api = projs[acc]["metadata"]
        base = "%s/systemes/%s" % (ROOT, pdb)
        reps = sorted(int(x[3:]) for x in os.listdir(base) if x.startswith("rep"))
        T = {}
        SY = {"pdb": pdb, "accession": acc, "repliques": reps, "date": now_iso()}
        for k in reps:
            d = "%s/rep%d" % (base, k)
            try:
                r, t = verifier_replique(d, api)
                T[k] = t
            except Exception as ex:  # noqa: BLE001
                r = {"dossier": d, "date": now_iso(), "anomalies": ["verification impossible: %s: %s" % (type(ex).__name__, ex)]}
            r.update({"pdb": pdb, "accession": acc, "replique": k, "conforme": not r["anomalies"]})
            rep_out.write(json.dumps(r, ensure_ascii=False) + "\n")
            rep_out.flush()
            print("VERIF %s rep%d images=%s atomes=%s chaines_prot=%s residus=%s non_prot=%s ruptures=%s paires5A=%s ecart_xtc0_pdb=%s anomalies=%s" % (
                pdb, k, r.get("n_images"), r.get("n_atomes"), r.get("n_chaines_proteiques"), r.get("n_residus_proteine"),
                r.get("n_residus_non_proteine"), [x["images_avec_rupture_gt_4.5A"] for x in r.get("ca_ruptures_par_chaine", [])],
                r.get("paires_5A"), r.get("ecart_max_xtc0_vs_structurepdb_nm"), r["anomalies"] or "-"), flush=True)
        # ---- les repliques different-elles ?
        ent = {}
        for i in reps:
            for j in reps:
                if i < j and i in T and j in T and T[i].n_atoms == T[j].n_atoms:
                    ent["%d-%d" % (i, j)] = {
                        "image0_nm": round(float(np.abs(T[i].xyz[0] - T[j].xyz[0]).max()), 3),
                        "image500_nm": round(float(np.abs(T[i].xyz[500] - T[j].xyz[500]).max()), 3)}
        SY["ecart_entre_repliques_nm"] = ent
        SY["repliques_distinctes"] = bool(ent) and all(v["image0_nm"] > 0.1 and v["image500_nm"] > 0.1 for v in ent.values())
        print("ENTRE_REPLIQUES %s %s distinctes=%s" % (pdb, json.dumps(ent), SY["repliques_distinctes"]), flush=True)
        # ---- fichiers deja sur staging : production1.dcd = replique 1 ?
        dcd = "%s/%s/production1.dcd" % (STAGING_INPUTS, pdb)
        inp = "%s/%s/%s_input.pdb" % (STAGING_INPUTS, pdb, pdb)
        if os.path.exists(dcd) and os.path.exists(inp):
            SY["staging"] = {"dcd": dcd, "dcd_octets": os.path.getsize(dcd), "dcd_sha256": sha256_file(dcd),
                             "input_pdb": inp, "input_octets": os.path.getsize(inp), "input_sha256": sha256_file(inp)}
            cmp = {}
            for k in reps:
                if k not in T:
                    continue
                ec = []
                for fr in (0, 500, 999):
                    x = md.load_frame(dcd, fr, top=inp)
                    if x.n_atoms != T[k].n_atoms:
                        ec.append(None)
                    else:
                        ec.append(round(float(np.abs(x.xyz[0] - T[k].xyz[fr]).max()), 4))
                cmp["rep%d" % k] = ec
            SY["staging"]["dcd_vs_xtc_nm_images_0_500_999"] = cmp
            r1 = cmp.get("rep1")
            SY["staging"]["dcd_est_replique_1"] = bool(r1 and all(v is not None and v < 0.002 for v in r1)
                                                       and all(all(v is None or v > 0.1 for v in cmp[x]) for x in cmp if x != "rep1"))
            print("DCD %s %s dcd_est_replique_1=%s" % (pdb, json.dumps(cmp), SY["staging"]["dcd_est_replique_1"]), flush=True)
        systemes[pdb] = SY
    rep_out.close()
    jdump({"date": now_iso(), "systemes": systemes}, CACHE + "/verif_systemes%s.json" % SUF)
    print("FINI s3b secondes=%.1f %s" % (time.time() - t0, now_iso()), flush=True)


if __name__ == "__main__":
    main()
