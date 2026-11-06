#!/usr/bin/env python3
"""h_piles.py — assemble les piles `[reference H + 25 echantillons H]` (format de assembler_pile.py / reconstruire_H.assembler : `MODEL %4d`, lignes ATOM,
`ENDMDL`, `END`) et REJOUE les gardes AVANT DynaPIN, independamment de l outil :
  42 restitution : pour CHAQUE image, h_core.ecart_lourds_texte(entree, R3) doit rendre 0.0 (appariement ORDINAL) ; 40 0 hydrogene ; 41 topologie variable
  (garde_H.verifier, copie inchangee de hydrogenes/orig) ; 43 nombre d images R3 != nombre d entrees.
Un systeme/unite en echec est une ENTREE du manifeste (type pile_H, statut echec), jamais un saut ; il n a pas de pile, donc pas de DynaPIN.
A lancer avec le python de l env evomd (h_core importe openmm ; aucun Context n est cree)."""
import hashlib, json, os, sys, time
from pathlib import Path
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/hydrogenes/scripts")
sys.path.insert(0, "/gpfs/phherman/bioemu2_bench/hydrogenes/orig")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import h_core as H
import garde_H
import gen_lib as L

US = os.environ.get("H_USAGE", "U3")
HI = L.G / os.environ.get("H_RACINE", "hydrogenes")
OUT, PILES = HI / "out", HI / "piles"
PILES.mkdir(exist_ok=True)
JOB = os.environ.get("SLURM_JOB_ID")
taches = [l.split("\t")[0] for l in (HI / "taches.tsv").read_text().splitlines() if l.strip()]
if len(sys.argv) > 1:
    taches = [t for t in taches if t in sys.argv[1:]]          # reprise : seulement ces systemes
n_ok = n_ko = 0
t0 = time.time()


def atomes(p):
    return [l for l in open(p) if l.startswith("ATOM")]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


for pdb in taches:
    r3 = OUT / pdb / "H_t10_i2000" / "R3"
    ref_r3, ref_in = r3 / f"ref_{pdb}" / "img_000.pdb", HI / "in" / f"ref_{pdb}" / "img_000.pdb"
    unites = sorted(d.name for d in (HI / "in").glob(f"{pdb}__{US}__*") if d.is_dir())
    for u in unites:
        pile = PILES / u
        ligne = {"type": "pile_H", "unite": u, "systeme": pdb, "usage": US, "regime": "pile_H", "jobid": JOB, "pile": str(pile)}
        try:
            entrees = sorted((HI / "in" / u).glob("img_*.pdb"))
            sorties = sorted((r3 / u).glob("img_*.pdb")) if (r3 / u).is_dir() else []
            if not ref_r3.is_file() or len(sorties) != len(entrees) or len(entrees) == 0:
                raise RuntimeError(f"43 : {len(sorties)} images R3 pour {len(entrees)} entrees (reference R3 presente : {ref_r3.is_file()})")
            ecart = 0.0
            for pin, pout in [(ref_in, ref_r3)] + list(zip(entrees, sorties)):
                if pin.name != pout.name:
                    raise RuntimeError(f"43 : {pin.name} apparie a {pout.name}")
                e, nc, sans = H.ecart_lourds_texte(str(pin), str(pout))
                ecart = max(ecart, e)
                if e != 0.0:
                    raise RuntimeError(f"42 : lourds non restitues a 0,000 A dans {pout} (ecart {e} A, {sans} atomes sans correspondant)")
            pile.mkdir(parents=True, exist_ok=True)
            with open(pile / "pile.pdb", "w") as f:
                for k, p in enumerate([ref_r3] + sorties):
                    f.write("MODEL     %4d\n" % k)
                    f.write("".join(atomes(p)))
                    f.write("ENDMDL\n")
                f.write("END\n")
            g = garde_H.verifier(str(pile / "pile.pdb"))
            (pile / "garde_H.txt").write_text(f"code {g['code']} · {g['n_modeles']} modeles · H {g['H_min']}-{g['H_max']} · {g['atomes_par_modele']} atomes · {g['differences']}\n")
            ligne.update(code_garde_H=g["code"], n_modeles=g["n_modeles"], H_min=g["H_min"], H_max=g["H_max"], atomes_par_modele=g["atomes_par_modele"],
                         ecart_lourds_restitues_max_A=ecart, sha256_pile=sha(pile / "pile.pdb"))
            if g["code"] != 0:
                raise RuntimeError(f"{g['code']} : garde_H {g['differences']}")
            ligne["statut"] = "ok"
            n_ok += 1
        except Exception as e:
            ligne.update(statut="echec", erreur=str(e)[:300])
            (pile / "dynapin.code").unlink(missing_ok=True) if pile.is_dir() else None
            if pile.is_dir():
                (pile / "pile_a_ne_pas_juger.txt").write_text(str(e)[:300] + "\n")
                (pile / "pile.pdb").rename(pile / "pile.pdb.refusee") if (pile / "pile.pdb").is_file() else None
            n_ko += 1
        L.ajouter_manifeste(ligne)
        print(f"H_PILE {u} {ligne['statut']} {ligne.get('erreur', '')}", flush=True)
print(f"H_PILES_FIN ok={n_ok} echec={n_ko} (rejets = {n_ko}) {time.time() - t0:.0f} s")
sys.exit(0 if n_ko == 0 else 1)
