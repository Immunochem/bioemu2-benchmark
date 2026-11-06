#!/bin/bash
#SBATCH --job-name=rd_neg
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/redocking/logs/%x_%j.log
# CPU. Construit les paires du TEMOIN NEGATIF (chaines de deux complexes differents, sans parente de sequence) et les passe par la curation. Variables : NEG (nom de strate), N, SOURCES (strates sources).
# --time 20 min : plafond DECLARE, non mesure (16 paires x curation ~1 s + alignements locaux de quelques centaines de residus, estimation ~2 min).
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/redocking/scripts/commun_redock.sh
RUN=$R/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_redock_cpu "redock__negatifs" "$RUN"
"$PYB" -B "$R/scripts/construire_negatifs.py" ${NEG:?} ${N:?} ${SOURCES:?} > "$RUN/neg.stdout" 2> "$RUN/neg.stderr"; CODE=$?
echo "=== construire_negatifs -> $CODE · $(date -Is)"
cut -c1-260 "$RUN/neg.stdout"
if [ "$CODE" -ne 0 ]; then tail -20 "$RUN/neg.stderr"; exit "$CODE"; fi
[ -s "$R/strate_${NEG}/manifeste.json" ] || { echo "ARTEFACT ABSENT"; exit 30; }
exit 0
