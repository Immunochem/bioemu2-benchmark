#!/usr/bin/env python3
"""h_reordonner.py PDB [PDB ...] — liste d entrees reordonnee pour les systemes dont la 1re image de certains workers fait echouer l outil.
POURQUOI (MESURE, 1AY7 : « No template found for residue 6 (CYS) … 1 S atom too many », 7 workers sur 8) : l outil cree UN Context par worker, sur la topologie de la
PREMIERE image qu il traite ; or OpenMM (Topology.createDisulfideBonds, seuil 0,3 nm) infere une liaison S–S entre deux SG a < 3,0 A, quelle que soit la variante
imposee par le PSF. Une image dont deux cysteines libres sont a < 3,0 A (cas A) ou dont un pont declare est a > 3,0 A (cas B) ne peut donc pas OUVRIR un worker ; elle passe
sans dommage quand elle vient apres (le Context reutilise garde la topologie de la 1re image, donc l etat du PSF).
On ne touche NI a l outil NI aux gardes : on met en tete la reference puis 7 images NON touchees (indices 0..7 de la liste = 1res images des 8 shards), le reste suit dans l ordre.
Sortie : generation/hydrogenes/listes_v2/<PDB>.txt et taches_reprise.tsv ; la liste d origine est conservee. Aucune image n est retiree : le compte est verifie."""
import json, os, sys
from pathlib import Path
G = Path("/gpfs/phherman/bioemu2_bench/generation")
HI = G / os.environ.get("H_RACINE", "hydrogenes")
diag = json.loads((HI / "diag_ss.json").read_text())["systemes"]
(HI / "listes_v2").mkdir(exist_ok=True)
tsv_orig = {l.split("\t")[0]: l.split("\t") for l in (HI / "taches.tsv").read_text().splitlines() if l.strip()}
lignes = []
systemes = sys.argv[1:]
if systemes == ["TOUS"]:
    systemes = [l.split("\t")[0] for l in (HI / "taches.tsv").read_text().splitlines() if l.strip()]
for pdb in systemes:
    liste = [l.strip() for l in (HI / "listes" / f"{pdb}.txt").read_text().splitlines() if l.strip()]
    touchees = set()
    for u, v in diag[pdb]["unites"].items():
        for f in v["images_touchees"]:
            touchees.add(str(HI / "in" / u / f"img_{f:03d}.pdb"))
    bons = [p for p in liste[1:] if p not in touchees]
    tete = [liste[0]] + bons[:7]
    k = len(tete)          # nombre de workers qui peuvent DEMARRER sur une image propre (reference + images propres, au plus 8) ; si k < 8 l outil tourne avec k workers sur ce systeme
    neuf = tete + [p for p in liste if p not in set(tete)]
    assert sorted(neuf) == sorted(liste) and len(neuf) == len(liste), "le reordonnancement a perdu ou duplique une image"
    (HI / "listes_v2" / f"{pdb}.txt").write_text("\n".join(neuf) + "\n")
    lignes.append(f"{pdb}\t{HI / 'listes_v2' / (pdb + '.txt')}\t{tsv_orig[pdb][2]}\t{k}")
    print(f"REORDONNE {pdb} : {len(liste)} images, {len(touchees)} touchees, {k} premieres propres -> {k} workers")
(HI / "taches_reprise.tsv").write_text("\n".join(lignes) + "\n")
