import json, sys, glob
import numpy as np
R = "/gpfs/phherman/bioemu2_bench/redocking/sorties"
for d in sorted(glob.glob(R + "/*/*/meta.json")):
    meta = json.load(open(d)); base = d.rsplit("/", 1)[0]
    if meta.get("statut") != "ok": continue
    c = json.load(open(base + "/conditions.json")); m = json.load(open(base + "/lag0050/mesures.json"))["mesures"]
    print(f"== {meta['id']} {meta['strate']} tok={meta['tokens']} d_nat={meta['d_nat_A']} n_nat={meta['n_contacts_natifs']} s/struct={meta['lags']['50']['s_par_structure']} conc={meta['lags']['50']['jobs_concurrents_gpu']['min']}-{meta['lags']['50']['jobs_concurrents_gpu']['max']}")
    by = {}
    for a in c["audit"]:
        if a["accepte"]: by.setdefault(a["cond"], []).append(a["indice_image"])
    for k, ii in by.items():
        g = lambda key: np.median([m[key][i] for i in ii])
        print(f"   {k:7s} fnat {g('fnat'):.3f} iRMSD {g('irmsd_A'):5.2f} lRMSD {g('lrmsd_A'):6.2f} DockQ {g('dockq'):.3f} d_COM {g('d_com_A'):5.1f} (x0 {np.median([ [a for a in c['audit'] if a.get('indice_image')==i][0]['d_com_A'] for i in ii]):5.1f}) dmin {min(m['dmin_inter_A'][i] for i in ii):.2f} classes {[m['classe_capri'][i][:3] for i in ii]}")
