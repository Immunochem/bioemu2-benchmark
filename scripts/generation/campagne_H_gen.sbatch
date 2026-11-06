#!/bin/bash
#SBATCH --job-name=gen_H
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# DERIVE de hydrogenes/scripts/campagne_H.sbatch (gabarit « non soumis ») — ECARTS DECLARES : (1) sortie du journal et racine sous generation/ (jamais d ecriture
# sous hydrogenes/) ; (2) G13 : le job refuse de demarrer s il y a plus de PREVUS calculs GPU, et echantillonne la concurrence ; (3) temps et code par systeme ecrits
# pour que h_bilan.py ajoute UNE ligne `etage_H` par systeme au manifeste ; (4) apres la boucle, h_bilan.py. La boucle, l outil et ses options sont ceux du gabarit.
# Arguments : PREVUS [WORKERS] [TACHES.tsv]. Le --time est passe en ligne de commande.
# --time 02:00:00 : debit MESURE par l agent hydrogenes, 211,6 s pour 277 images de 6CWG a 8 workers (0,764 s/image), 1,87 h estimees pour 28 x 277 images ;
# ici 28 x 126 = 3 528 images -> 1,87 h x 126/277 = 0,85 h x 1,35 = 1,15 h -> arrondi a l heure = 2 h. --reprise : une image deja reconstruite ET verifiee est sautee.
set -uo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
Hd=/gpfs/phherman/bioemu2_bench/hydrogenes; B=/gpfs/phherman/bioemu2_bench; GEN=$B/generation; HR="${H_RACINE:-hydrogenes}"
PY=/home/phherman/micromamba/envs/evomd/bin/python
PYB=/gpfs/phherman/bioemu2/venv/bin/python
PREVUS="${1:?PREVUS}"; W="${2:-8}"
TACHES="${3:-$GEN/$HR/taches.tsv}"; RACINE=$GEN/$HR/out
RUN=$GEN/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID; mkdir -p "$RUN" "$RACINE"
if [ -z "${CUDA_VISIBLE_DEVICES:-}" ]; then echo "ARRET 33 : CUDA_VISIBLE_DEVICES vide : aucune slice visible"; exit 33; fi
if touch "$RUN/.essai"; then rm -f "$RUN/.essai"; else echo "ARRET 13 : ecriture impossible sous $RUN"; exit 13; fi
TOT=$(python3 $B/L0b/outils/concurrence.py maintenant --champ total); CQ=$?
[ "$CQ" -eq 0 ] || { echo "ARRET 22 : echantillon de concurrence impossible"; exit 22; }
if [ "$TOT" -gt "$PREVUS" ]; then echo "ARRET 13 : $TOT calcul(s) GPU en cours (moi compris), $PREVUS prevu(s)."; exit 13; fi
echo "=== campagne_H $SLURM_JOB_ID $(hostname) $(date -Is) ; CUDA_VISIBLE_DEVICES=[$CUDA_VISIBLE_DEVICES] ; $(wc -l < "$TACHES") systemes ; $W workers ; G13 au demarrage : $TOT (prevu $PREVUS)"
nvidia-smi -L
sha256sum $Hd/scripts/reconstruire_H_gpu.py $Hd/scripts/h_core.py $Hd/scripts/executer_variante.py $Hd/orig/garde_H.py | cut -c1-16,65- > "$RUN/sha_outil.txt"
cat "$RUN/sha_outil.txt"
python3 $B/L0b/outils/concurrence.py echantillonner --sortie "$RUN/concurrence.jsonl" --intervalle 15 --mode staging &
SAMP=$!
BAD=0
: > "$RUN/systemes_temps.tsv"
while IFS=$'\t' read -r SYS LISTE PSF WK; do
  [ -n "$SYS" ] || continue
  WKS="${WK:-$W}"; [ "$WKS" -le "$W" ] || WKS=$W
  T0=$(date +%s)
  $PY $Hd/scripts/reconstruire_H_gpu.py --entrees "$LISTE" --sortie "$RACINE/$SYS" --psf "$PSF" --workers "$WKS" --reprise > "$RACINE/$SYS.log" 2>&1; C=$?
  T1=$(date +%s)
  printf '%s\t%s\t%s\t%s\n' "$SYS" "$T0" "$T1" "$C" >> "$RUN/systemes_temps.tsv"
  echo "$SYS code=$C workers=$WKS $((T1-T0)) s $(grep -E 'ARRET|ECHEC' "$RACINE/$SYS.log" | head -1 | cut -c1-160)"
  [ "$C" -eq 0 ] || BAD=$((BAD+1))
done < "$TACHES"
kill -TERM "$SAMP"; wait "$SAMP"
echo "=== $BAD systeme(s) en echec ; $(date -Is)"
python3 $B/L0b/outils/concurrence.py resume --fichier "$RUN/concurrence.jsonl" > "$RUN/concurrence_resume.json"; CR=$?
echo "=== G13 apres coup : $(python3 -c "import json;r=json.load(open('$RUN/concurrence_resume.json'));print('max', r.get('max'), 'fractions', r.get('fraction_du_temps_par_valeur'))") (code resume $CR)"
"$PYB" $GEN/scripts/h_bilan.py "$RACINE" "$RUN"; CB=$?
echo "=== h_bilan -> $CB"
exit "$BAD"
