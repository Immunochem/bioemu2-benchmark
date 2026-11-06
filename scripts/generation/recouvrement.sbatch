#!/bin/bash
#SBATCH --job-name=gen_recouvrement
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# RMSD Calpha croise U2 x U3-lag1000 (CPU). --time 01:00:00 : debit MESURE voisin — relecture de 144 dossiers de sorties en < 1 min (job 11121) ; 28 x (50 x 50) RMSD de 100-800 Calpha ~ secondes -> plancher d une heure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__recouvrement_U2_U3" "$RUN"
"$PYB" "$G/scripts/recouvrement.py" > "$RUN/rec.stdout" 2> "$RUN/rec.stderr"; CODE=$?
echo "=== recouvrement -> $CODE"; tail -4 "$RUN/rec.stdout" | cut -c1-300
if [ "$CODE" -ne 0 ]; then tail -20 "$RUN/rec.stderr"; fi
exit "$CODE"
