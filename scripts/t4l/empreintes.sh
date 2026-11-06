#!/bin/bash
# empreintes.sh : liste markdown des fichiers cles avec taille et sha256 (16 premiers caracteres). Les sha se calculent ICI, sur staging.
cd /gpfs/phherman/bioemu2_bench/t4l || exit 1
echo "| fichier | octets | sha256 (16) |"; echo "|---|---|---|"
for f in scripts/cv_t4l.py scripts/reproduire_v1.py scripts/proxy_cb.py scripts/inspecter_5lzm.py scripts/bioemu2_t4l.py scripts/bioemu2_t4l_steer.py scripts/verifier_enligne.py scripts/analyser_t4l.py scripts/rmsd_ferme.py scripts/rmsd_cristal.py scripts/ks_unites.py scripts/refs_fermees_cv.py scripts/tables_rapport.py scripts/analyse_tout.sh scripts/j1_t4l.sbatch scripts/j2_t4l.sbatch scripts/j3_t4l.sbatch \
  plans/smoke.json plans/prod_a.json plans/prod_b.json plans/steer.json \
  x0/5LZM_proteine.pdb x0/x0_5lzm.pdb x0/rapport_curation_5lzm.json msa/profondeur.json \
  donnees/zenodo/T4L-bechmarking.zip donnees/github/repo.tar.gz donnees/pdb5lzm/5LZM.pdb donnees/notebooks/src_github_analysis.txt donnees/notebooks/src_zenodo_analysis.txt \
  resultats/reproduction_v1.json resultats/proxy_cb.json resultats/5lzm_inspection.json resultats/refs_fermees_cv.json resultats/analyse_unites.json resultats/rmsd_ferme_unites.json resultats/rmsd_cristal.json resultats/ks_unites.json resultats/verif_enligne_unites.json resultats/tables_rapport.md \
  resultats/analyse_unites_steer.json resultats/rmsd_ferme_unites_steer.json resultats/tables_rapport_unites_steer.md \
  figs/fig_comparaison_unites.png figs/fig_comparaison_unites_steer.png \
  unites/eq164_n500/cv.npz unites/eq164_n500/brut/samples.xtc unites/eq164_n500/meta.json; do
  if [ -f "$f" ]; then printf "| \`%s\` | %s | \`%s\` |\n" "$f" "$(stat -c %s "$f")" "$(sha256sum "$f" | cut -c1-16)"; else printf "| \`%s\` | ABSENT | — |\n" "$f"; fi
done
