#!/usr/bin/env python3
"""rmsd_ferme.py — d1 < 2,5 nm suffit-il a dire « etat ferme » ?  Non : d1 est UNE distance. On regarde si une structure a d1 < 2,5 nm RESSEMBLE aux 221 structures
fermees de l'auteur (118 exposed-closed + 103 buried-closed, tirees de sa MD, 'trajfit-cluster<N>_frame<M>.pdb').

Mesure : RMSD des Ca (appariement ORDINAL, residus 1-162 : le plus petit commun denominateur 5LZM/papier), apres superposition optimale (md.rmsd).
Etalon (vrai negatif de la distance) : pour chaque reference, son RMSD a la plus proche reference d'un AUTRE cluster de l'auteur = ce qu'un etat ferme authentique
mais d'un autre bassin affiche. Une structure « fermee » de BioEmu2 qui reste plus loin que ce 95e centile n'est pas reconnue comme son etat ferme.
Controle de coherence : le RMSD moyen au cristal 5LZM (ouvert) des references fermees dit a quelle distance se trouve l'etat ouvert.
Usage : rmsd_ferme.py <racine unites> <unite> [<unite> ...]"""
import glob, json, os, re, sys
import numpy as np, mdtraj as md
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C
T = "/gpfs/phherman/bioemu2_bench/t4l"
Z = f"{T}/donnees/zenodo/extrait/T4L-bechmarking/Closed-state"
N = 162


def ca162(traj):
    top = traj.topology
    res = list(list(top.chains)[0].residues)
    idx = [a.index for r in res[:N] for a in r.atoms if a.name == "CA"]
    assert len(idx) == N, len(idx)
    return traj.atom_slice(idx)


racine, unites = sys.argv[1], sys.argv[2:]
refs, clus, noms = [], [], []
for lab in ("Exposed-closed", "Buried-closed"):
    for f in sorted(glob.glob(f"{Z}/{lab}/*.pdb")):
        t = md.load(f)
        C.indices(t.topology)                                   # garde ordinal
        refs.append(ca162(t)); clus.append(int(re.search(r"cluster(\d+)_", f).group(1))); noms.append(f"{lab}/{os.path.basename(f)}")
R = md.join(refs)
clus = np.array(clus)
print(f"references fermees : {R.n_frames}, clusters {sorted(set(clus.tolist()))}")
# etalon : plus proche reference d'un AUTRE cluster
etalon = []
for i in range(R.n_frames):
    d = md.rmsd(R, R, frame=i) * 10.0
    d[clus == clus[i]] = np.inf
    etalon.append(d.min())
etalon = np.array(etalon)
q = np.percentile(etalon, [50, 95, 100])
cris = ca162(md.load(f"{T}/x0/x0_5lzm.pdb"))
d_cris = md.rmsd(R, cris, frame=0) * 10.0
out = {"n_refs": int(R.n_frames), "etalon_autre_cluster_A_med_95_max": [round(float(x), 2) for x in q],
       "refs_vers_cristal_5LZM_A_min_med_max": [round(float(d_cris.min()), 2), round(float(np.median(d_cris)), 2), round(float(d_cris.max()), 2)], "unites": {}}
print(json.dumps({k: v for k, v in out.items() if k != "unites"}))
for u in unites:
    d = f"{T}/{racine}/{u}"
    tr = md.load_xtc(f"{d}/brut/samples.xtc", top=f"{d}/brut/topology.cif")
    z = np.load(f"{d}/cv.npz")
    A = ca162(tr)
    mn = np.array([ (md.rmsd(R, A, frame=i) * 10.0).min() for i in range(A.n_frames)])
    d_open = md.rmsd(A, cris, frame=0) * 10.0
    f = z["d1"] < 2.5
    e = {"n": int(A.n_frames), "n_ferme_d1": int(f.sum()),
         "rmsd_min_aux_refs_fermees_A_quantiles_0_25_50_75_100_toutes": [round(float(x), 2) for x in np.percentile(mn, [0, 25, 50, 75, 100])],
         "rmsd_au_cristal_A_mediane_toutes": round(float(np.median(d_open)), 2)}
    # NB : le critere « sous le 95e centile de l'etalon » (11,6 A) est SATISFAIT PAR L'ECHEC : les structures OUVERTES non guidees (6,1 A) le passent toutes.
    # On le rapporte quand meme, a cote de son vrai negatif, et on juge sur la DISTRIBUTION du RMSD min et sur le seuil median de l'etalon.
    ouv = ~f
    e["vrai_negatif_ouvertes_sous_etalon_95"] = int((mn[ouv] <= q[1]).sum()) if ouv.any() else None
    if f.any():
        e["fermees"] = {"rmsd_min_aux_refs_A_min_med_max": [round(float(mn[f].min()), 2), round(float(np.median(mn[f])), 2), round(float(mn[f].max()), 2)],
                        "n_sous_etalon_95_(NON DISCRIMINANT)": int((mn[f] <= q[1]).sum()), "etalon_95_A": round(float(q[1]), 2),
                        "n_sous_etalon_mediane": int((mn[f] <= q[0]).sum()), "etalon_mediane_A": round(float(q[0]), 2),
                        "rmsd_au_cristal_A_min_med_max": [round(float(d_open[f].min()), 2), round(float(np.median(d_open[f])), 2), round(float(d_open[f].max()), 2)],
                        "d2_nm_min_med_max": [round(float(z["d2"][f].min()), 3), round(float(np.median(z["d2"][f])), 3), round(float(z["d2"][f].max()), 3)],
                        "p_median": round(float(np.median(z["p"][f])), 3), "exposed_closed": int((z["p"][f] > 0).sum()), "buried_closed": int((z["p"][f] < 0).sum()),
                        "physiques": int(z["physique"][f].sum())}
    fp = f & z["physique"].astype(bool)
    if fp.any():
        e["fermees_ET_physiques"] = {"n": int(fp.sum()), "rmsd_min_aux_refs_A_min_med_max": [round(float(mn[fp].min()), 2), round(float(np.median(mn[fp])), 2), round(float(mn[fp].max()), 2)],
                                     "n_sous_etalon_mediane": int((mn[fp] <= q[0]).sum()), "rmsd_au_cristal_A_min_med_max": [round(float(d_open[fp].min()), 2), round(float(np.median(d_open[fp])), 2), round(float(d_open[fp].max()), 2)],
                                     "d1_nm_min_med_max": [round(float(z["d1"][fp].min()), 3), round(float(np.median(z["d1"][fp])), 3), round(float(z["d1"][fp].max()), 3)],
                                     "d2_nm_min_med_max": [round(float(z["d2"][fp].min()), 3), round(float(np.median(z["d2"][fp])), 3), round(float(z["d2"][fp].max()), 3)],
                                     "p_median": round(float(np.median(z["p"][fp])), 3), "exposed_closed": int((z["p"][fp] > 0).sum()), "buried_closed": int((z["p"][fp] < 0).sum())}
    if ouv.any():
        e["ouvertes_rmsd_min_aux_refs_A_min_med_max"] = [round(float(mn[ouv].min()), 2), round(float(np.median(mn[ouv])), 2), round(float(mn[ouv].max()), 2)]
    out["unites"][u] = e
    print(u, json.dumps(e))
json.dump(out, open(f"{T}/resultats/rmsd_ferme_{racine}.json", "w"), indent=1)
