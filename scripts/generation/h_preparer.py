#!/usr/bin/env python3
"""h_preparer.py — entrees de l etage hydrogenes (CPU) : une image = un fichier `hydrogenes/in/<unite>/img_NNN.pdb` (atomes LOURDS, un MODEL),
plus une image de reference `in/ref_<PDB>/img_000.pdb` (le x0 cure) par systeme. Format et ecriture DERIVES de 6CWG/scripts/reconstruire_H.py
(`cles`, `ecrire_modele`) — ecart declare : les entrees viennent de generation/sorties/ (U3), les noms d atomes du x0 CURE servent de modele.
🔴 Appariement des atomes ORDINAL : cle = (index de chaine, ordinal 1-based du residu dans la chaine, nom d atome). Jamais resSeq.
Un echantillon dont l ensemble d atomes lourds != celui du x0 est un REJET (entree au manifeste), jamais un appariement partiel.
"""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L
import mdtraj as md

import os
HI = L.G / os.environ.get("H_RACINE", "hydrogenes")
IN, LISTES = HI / "in", HI / "listes"
for d in (IN, LISTES):
    d.mkdir(parents=True, exist_ok=True)
USAGE = sys.argv[1] if len(sys.argv) > 1 else "U3"


def cles(top):
    out = {}
    for ch in top.chains:
        for o, r in enumerate(ch.residues, start=1):
            for at in r.atoms:
                if at.element.symbol != "H":
                    out[(ch.index, o, at.name)] = at.index
    return out


def ecrire_modele(top, xyz, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".tmp")
    md.Trajectory(xyz=xyz[None], topology=top).save_pdb(str(tmp))
    atomes = [l for l in tmp.read_text().splitlines(True) if l.startswith("ATOM")]
    tmp.unlink()
    path.write_text("MODEL        0\n" + "".join(atomes) + "ENDMDL\n")


bilan, rejets, tsv = {"systemes": 0, "images": 0, "unites": 0}, [], []
t0 = time.time()
for pdb in sorted(L.corpus_28()):
    try:
        x0 = md.load(str(L.chemin_x0(pdb, 1)))
        L.verifier_x0(pdb, 1)
        kref = cles(x0.topology)
        ordre = sorted(kref, key=lambda c: kref[c])
        i_ref = [kref[c] for c in ordre]
        ecrire_modele(x0.topology, x0.xyz[0][i_ref], IN / f"ref_{pdb}" / "img_000.pdb")
        liste = [str(IN / f"ref_{pdb}" / "img_000.pdb")]
        unites = sorted(d for d in (L.G / "sorties").glob(f"{pdb}__{USAGE}__*") if d.is_dir() and not d.name.endswith(".tmp") and "perime" not in d.name)
        if len(unites) == 0:
            raise L.Rejet(f"aucune unite {USAGE} pour {pdb}")
        for ud in unites:
            meta = json.loads((ud / "meta.json").read_text())
            if meta.get("statut") != "ok":
                raise L.Rejet(f"{ud.name} : meta.json statut {meta.get('statut')}")
            tr = md.load_xtc(str(ud / "samples.xtc"), top=str(ud / "topology.cif"))
            ke = cles(tr.topology)
            if set(ke) != set(kref):
                raise L.Rejet(f"{ud.name} : {len(set(ke) ^ set(kref))} atomes lourds en desaccord ordinal avec le x0 (ref {len(kref)}, echantillon {len(ke)})")
            i_ech = [ke[c] for c in ordre]
            if tr.n_frames != meta["n"]:
                raise L.Rejet(f"{ud.name} : {tr.n_frames} images dans le xtc, n = {meta['n']}")
            for i in range(tr.n_frames):
                p = IN / ud.name / f"img_{i:03d}.pdb"
                if not (p.is_file() and p.stat().st_size > 0):
                    ecrire_modele(x0.topology, tr.xyz[i][i_ech], p)
                liste.append(str(p))
                bilan["images"] += 1
            bilan["unites"] += 1
        (LISTES / f"{pdb}.txt").write_text("\n".join(liste) + "\n")
        tsv.append(f"{pdb}\t{LISTES / (pdb + '.txt')}\t{L.BENCH / 'corpus/systemes' / pdb / 'rep1/topology.psf'}")
        bilan["systemes"] += 1
        print(f"H_PREP {pdb} : {len(unites)} unites, {len(liste)} images (reference comprise), {len(kref)} atomes lourds", flush=True)
    except Exception as e:
        rejets.append({"systeme": pdb, "raison": f"{type(e).__name__}: {str(e)[:300]}"})
        L.ajouter_manifeste({"type": "rejet", "unite": f"{pdb}__{USAGE}__etage_H_preparation", "systeme": pdb, "statut": "rejet_preparation_H",
                             "motif": f"{type(e).__name__}: {str(e)[:300]}", "jobid": os.environ.get("SLURM_JOB_ID")})
        print(f"H_PREP REJET {pdb} : {e}", flush=True)
(HI / "taches.tsv").write_text("\n".join(tsv) + "\n")
bilan.update(rejets=rejets, n_rejets=len(rejets), duree_s=round(time.time() - t0, 1), usage=USAGE)
(HI / "preparation.json").write_text(json.dumps(bilan, indent=1, ensure_ascii=False) + "\n")
print("H_PREP_BILAN", json.dumps({k: v for k, v in bilan.items() if k != "rejets"}), f"(rejets = {len(rejets)})")
sys.exit(0 if not rejets else 1)
