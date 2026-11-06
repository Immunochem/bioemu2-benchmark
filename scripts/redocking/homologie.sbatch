#!/bin/bash
#SBATCH --job-name=rd_homol
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/redocking/logs/%x_%j.log
# Reseau seul (API de recherche RCSB, une requete a la fois, pause 1 s) : 30 systemes B x 2 chaines x 3 seuils = 180 recherches de sequence + pagination. --time : plafond DECLARE, non mesure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/redocking/scripts/commun_redock.sh
RUN=$R/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_redock_cpu "redock__homologie_B" "$RUN"
"$PYB" -B "$R/scripts/homologie_B.py" > "$RUN/homol.stdout" 2> "$RUN/homol.stderr"; CODE=$?
echo "=== homologie_B -> $CODE · $(date -Is)"
grep -E "^BILAN_HOMOLOGIE|^ERREUR" "$RUN/homol.stdout"
if [ "$CODE" -ne 0 ]; then tail -15 "$RUN/homol.stderr"; fi
exit "$CODE"
