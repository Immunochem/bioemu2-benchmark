#!/bin/bash
#SBATCH --job-name=gen_ponts2
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Memes ponts que reference.json (NON ecrase : les workers de production le lisent) + quantiles poolés de la MD et NUL du critere
# « hors plage MD » -> ponts/reference_v2.json. --time 01:00:00 : debit MESURE du meme script, job 11091 = 10,6 s (duree_s du bilan) x 1,35 -> 1 h.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__ponts_reference_v2" "$RUN"
"$PYB" "$G/scripts/ponts_reference.py" --sortie reference_v2.json > "$RUN/ponts.stdout" 2> "$RUN/ponts.stderr"; CODE=$?
echo "=== ponts_reference v2 -> $CODE"; tail -3 "$RUN/ponts.stdout" | cut -c1-400
if [ "$CODE" -ne 0 ]; then tail -25 "$RUN/ponts.stderr"; fi
exit "$CODE"
