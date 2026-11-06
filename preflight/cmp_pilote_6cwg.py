import sys, os
sys.dont_write_bytecode = True
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import numpy as np, mdtraj as md
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/prep_fix")
import curer_entree as ce
for nom_pilote, nom_mien in (("x0_md0.pdb", "6CWG__img0_rep1.pdb"), ("x0_cristal.pdb", "6CWG__input.pdb")):
    a = md.load(f"/gpfs/phherman/bioemu2_bench/6CWG/x0/{nom_pilote}")
    b = md.load(f"/gpfs/phherman/bioemu2_bench/prep_fix/curation/x0/{nom_mien}")
    ka, kb = ce.prep6.cles(a.topology), ce.prep6.cles(b.topology)
    communs = sorted(set(ka) & set(kb))
    d = np.array([np.linalg.norm(a.xyz[0][ka[k]] - b.xyz[0][kb[k]]) * 10 for k in communs])
    print(f"{nom_pilote} (pilote) vs {nom_mien} (curation) : atomes pilote {len(ka)} · curation {len(kb)} · communs {len(communs)} · seulement pilote {len(set(ka) - set(kb))} · seulement curation {len(set(kb) - set(ka))} · ecart max {d.max():.4f} A · moyen {d.mean():.5f} A")
