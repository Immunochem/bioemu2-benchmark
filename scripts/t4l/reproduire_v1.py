#!/usr/bin/env python3
"""reproduire_v1.py — LA PREMIERE MESURE : refaire, avec MON code (cv_t4l.py), les chiffres que l'auteur publie sur SES donnees.

Controles, dans l'ordre (chacun a un critere ecrit avant de lire le resultat) :
  C1  d1, d2 recalcules depuis les .xtc de Zenodo  == ses distances.npy (GitHub)          critere : ecart max < 1e-3 nm
  C2  p recalcule (cycle de Phe4)                  == ses locked.pkl (GitHub)             critere : ecart max < 1e-3
  C3  vrais positifs : 221 structures « Closed-state » de l'auteur (118 exposed-closed, 103 buried-closed)
        critere : d1 < 2,5 nm pour toutes ; p > 0 pour Exposed-closed ; p < 0 pour Buried-closed
  C4  son camembert « Gen-AI » (100 % exposed/open) : refait EXACTEMENT comme sa cellule [48] (lock = colonne 1 de distances.npy),
        puis refait avec le VRAI p (locked.pkl), puis avec mon p.
  C5  BioEmu v1 : d1 (ouvert/ferme) ; p est-il DEFINI ? (cycle de Phe4 present dans topology.pdb ?)
Un fichier rejete devient une ENTREE de `rejets` (jamais un saut silencieux) et le compte est imprime a chaque fois.
"""
import glob, io, json, os, pickle, sys
import numpy as np
import mdtraj as md

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_t4l as C

T = "/gpfs/phherman/bioemu2_bench/t4l"
Z = f"{T}/donnees/zenodo/extrait/T4L-bechmarking"
G = f"{T}/donnees/github/extrait/Gen-AI-t4l-bechmarking-main/GenAI-ensemble"
METHODES = {  # nom : (dossier Zenodo, xtc, pdb, dossier Analysis du depot)
    "AlphaFlow":      ("AlphaFlow", "combined.xtc", "prot-masses.pdb", f"{G}/AlphaFlow/Analysis"),
    "ESMFlow":        ("ESMFlow", "combined.xtc", "prot-masses.pdb", f"{G}/ESMFlow/Analysis"),
    "AF-cluster":     ("AF-cluster", "combined.xtc", "prot-masses.pdb", f"{G}/AF-cluster/Analysis"),
    "ConforFold":     ("ConforFold", "combined.xtc", "prot-masses.pdb", f"{G}/ConforFold/Analysis"),
    "ConfRover":      ("ConfRover", "T4L.xtc", "T4L.pdb", f"{G}/Confrover/Analysis"),
    "rMSA-AF2":       ("rMSA-AF2", "combined.xtc", "prot-masses.pdb", f"{G}/rmsa/pdb/Analysis"),
    "rMSA-AF2-recyc": ("rMSA-AF2-recyc", "combined.xtc", "prot-masses.pdb", f"{G}/rmsa/pdb/recycles"),
    "BioEmu":         ("BioEmu", "samples.xtc", "topology.pdb", f"{G}/BioEmu/Analysis"),
}
CINQ = ["AlphaFlow", "ESMFlow", "AF-cluster", "ConforFold", "ConfRover"]   # ce que concatene sa cellule [48]


class Restreint(pickle.Unpickler):
    OK = {("numpy.core.multiarray", "_reconstruct"), ("numpy._core.multiarray", "_reconstruct"), ("numpy", "ndarray"),
          ("numpy", "dtype"), ("numpy.core.multiarray", "scalar"), ("numpy._core.multiarray", "scalar"), ("builtins", "list"),
          ("_codecs", "encode")}   # protocole 2 : octets d'un ndarray encodes en latin1 ; encode() ne peut rien executer
    def find_class(self, module, name):
        if (module, name) in self.OK:
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f"REFUSE {module}.{name}")


def lire_pkl(chemin):
    obj = Restreint(io.BytesIO(open(chemin, "rb").read())).load()
    return np.concatenate([np.asarray(o).ravel() for o in obj]) if isinstance(obj, list) else np.asarray(obj).ravel()


def pct(c, n):
    return round(100.0 * c / n, 2) if n else None


def main():
    rejets, res, calc = [], {}, {}
    print("=== C1/C2 : mes d1, d2, p contre ses tableaux, methode par methode ===")
    for nom, (dz, xtc, pdb, dan) in METHODES.items():
        try:
            tr = md.load(f"{Z}/Genai-models/{dz}/{xtc}", top=f"{Z}/Genai-models/{dz}/{pdb}")
            r = C.cv(tr)
        except Exception as e:
            rejets.append({"methode": nom, "etape": "cv", "raison": repr(e)[:300]})
            print(f"REJET {nom} : {e!r}"[:300]); continue
        calc[nom] = r
        his = np.load(f"{dan}/distances.npy", allow_pickle=False)
        e = {"n_frames_xtc": tr.n_frames, "n_lignes_distances_npy": int(his.shape[0]), "n_atomes": tr.n_atoms,
             "resSeq_ordinaux": r["rapport"]["resSeq_ordinaux"], "cycle_phe4": r["rapport"]["atomes_cycle_phe4"],
             "p_defini": r["rapport"]["p_defini"], "diffs_sequence_vs_papier": r["rapport"]["diffs_sequence_vs_papier"]}
        if his.shape[0] == tr.n_frames:
            e["ecart_max_d1_nm"] = float(np.abs(r["d1"] - his[:, 0]).max())
            e["ecart_max_d2_nm"] = float(np.abs(r["d2"] - his[:, 1]).max())
        else:
            e["ecart_max_d1_nm"] = e["ecart_max_d2_nm"] = None
            rejets.append({"methode": nom, "etape": "C1", "raison": f"{tr.n_frames} images xtc contre {his.shape[0]} lignes distances.npy"})
        e["d1_min_med_max"] = [float(r["d1"].min()), float(np.median(r["d1"])), float(r["d1"].max())]
        e["d1_ferme_lt_2_5"] = int((r["d1"] < 2.5).sum())
        pk = f"{dan}/locked.pkl"
        if os.path.isfile(pk):
            ph = lire_pkl(pk)
            e["locked_pkl_n"] = int(ph.size)
            if ph.size == tr.n_frames and r["rapport"]["p_defini"]:
                e["ecart_max_p"] = float(np.abs(r["p"] - ph).max())
            else:
                e["ecart_max_p"] = None
            e["p_a_lui"] = {"min": float(ph.min()), "med": float(np.median(ph)), "max": float(ph.max()),
                            "p_gt_0": int((ph > 0).sum()), "p_lt_0": int((ph < 0).sum()), "p_eq_0": int((ph == 0).sum())}
        else:
            e["locked_pkl_n"] = None
        e["p_a_moi"] = ({"p_gt_0": int((r["p"] > 0).sum()), "p_lt_0": int((r["p"] < 0).sum()), "p_eq_0": int((r["p"] == 0).sum())}
                        if r["rapport"]["p_defini"] else "NON DEFINI : cycle de Phe4 absent du topology (p = NaN pour les %d images)" % tr.n_frames)
        res[nom] = e
        print(f"{nom:15s} n={tr.n_frames:5d} atomes={tr.n_atoms:5d} | d1 max|ecart|={e['ecart_max_d1_nm']} d2={e['ecart_max_d2_nm']} "
              f"| p defini={e['p_defini']} max|ecart p|={e.get('ecart_max_p')} | d1<2,5 : {e['d1_ferme_lt_2_5']}/{tr.n_frames}")

    print("\n=== C3 : vrais positifs, les 221 structures « Closed-state » de l'auteur ===")
    c3 = {}
    for lab, attendu_p in (("Exposed-closed", +1), ("Buried-closed", -1)):
        fich = sorted(glob.glob(f"{Z}/Closed-state/{lab}/*.pdb"))
        d1l, pl, rej = [], [], 0
        for f in fich:
            try:
                r = C.cv(md.load(f))
                d1l.append(float(r["d1"][0])); pl.append(float(r["p"][0]))
            except Exception as ex:
                rej += 1; rejets.append({"fichier": f, "etape": "C3", "raison": repr(ex)[:200]})
        d1a, pa = np.array(d1l), np.array(pl)
        c3[lab] = {"n_fichiers": len(fich), "n_rejets": rej, "n_calcules": len(d1l),
                   "d1_lt_2_5": int((d1a < 2.5).sum()), "d1_min_med_max": [float(d1a.min()), float(np.median(d1a)), float(d1a.max())] if len(d1a) else None,
                   "p_du_bon_signe": int(((pa * attendu_p) > 0).sum()), "p_min_med_max": [float(pa.min()), float(np.median(pa)), float(pa.max())] if len(pa) else None}
        print(lab, json.dumps(c3[lab]))

    print("\n=== C4 : le camembert « Gen-AI » (100 % exposed/open) ===")
    c4 = {}
    # (a) EXACTEMENT sa cellule [48] : all_genai = concat(distances.npy des CINQ) ; dist = col 0 ; lock = col 1 (= d2 !)
    allg = np.concatenate([np.load(f"{METHODES[m][3]}/distances.npy", allow_pickle=False) for m in CINQ], axis=0)
    dist, lock = allg[:, 0], allg[:, 1]
    cnt = [int(((lock > 0) & (dist < 2.5)).sum()), int(((lock > 0) & (dist > 2.5)).sum()), int(((lock < 0) & (dist < 2.5)).sum()), int(((lock < 0) & (dist > 2.5)).sum())]
    c4["a_sa_cellule_48_lock_egal_colonne_1"] = {"n": int(allg.shape[0]), "comptes_[exp/closed, exp/open, bur/closed, bur/open]": cnt,
                                                  "colonne_1_min": float(lock.min()), "colonne_1_est_toujours_positive": bool((lock > 0).all())}
    # (b) avec SON vrai p (locked.pkl)
    P = {m: lire_pkl(f"{METHODES[m][3]}/locked.pkl") for m in CINQ + ["rMSA-AF2", "rMSA-AF2-recyc"] if os.path.isfile(f"{METHODES[m][3]}/locked.pkl")}
    D = {m: np.load(f"{METHODES[m][3]}/distances.npy", allow_pickle=False) for m in P}
    for nom_jeu, liste in (("b_vrai_p_les_cinq", CINQ), ("b_vrai_p_sept_methodes", CINQ + ["rMSA-AF2", "rMSA-AF2-recyc"])):
        d1 = np.concatenate([D[m][:, 0] for m in liste if m in P]); pp = np.concatenate([P[m] for m in liste if m in P])
        pop = C.populations(d1, pp)
        c4[nom_jeu] = {"n": int(d1.size), **pop}
    # (c) avec mon p
    d1 = np.concatenate([calc[m]["d1"] for m in CINQ + ["rMSA-AF2", "rMSA-AF2-recyc"] if m in calc and calc[m]["rapport"]["p_defini"]])
    pp = np.concatenate([calc[m]["p"] for m in CINQ + ["rMSA-AF2", "rMSA-AF2-recyc"] if m in calc and calc[m]["rapport"]["p_defini"]])
    c4["c_mon_p_sept_methodes"] = {"n": int(d1.size), **C.populations(d1, pp)}
    print(json.dumps(c4, indent=1))

    print("\n=== C5 : BioEmu v1 ===")
    b = calc.get("BioEmu")
    if b is not None:
        res["BioEmu"]["resume"] = {"ouvert_d1_gt_2_5": int((b["d1"] > 2.5).sum()), "ferme_d1_lt_2_5": int((b["d1"] < 2.5).sum()),
                                   "p_defini": b["rapport"]["p_defini"], "resSeq_ordinaux": b["rapport"]["resSeq_ordinaux"]}
        print(json.dumps(res["BioEmu"]["resume"]))
    sortie = {"rejets": rejets, "n_rejets": len(rejets), "methodes": res, "C3": c3, "C4": c4}
    os.makedirs(f"{T}/resultats", exist_ok=True)
    json.dump(sortie, open(f"{T}/resultats/reproduction_v1.json", "w"), indent=1, ensure_ascii=False)
    print(f"\nREJETS : {len(rejets)}")
    for r in rejets:
        print(" ", r)
    print(f"ecrit : {T}/resultats/reproduction_v1.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
