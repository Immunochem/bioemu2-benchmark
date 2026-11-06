#!/bin/bash
#SBATCH --job-name=gen_ablanc
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# A BLANC de generer.py (CPU, aucun GPU, aucune requete) : systeme + MSA depuis le cache + x0 + step_scale lu par mode.
# --time 01:00:00 : le preflight (11062) a fait la meme chose pour 5 systemes en 13 s (sacct) x 1,35 -> arrondi a l heure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__ablanc_generer" "$RUN"
printf '%s\n' "$@" > "$RUN/taches.txt"
"$PYB" "$G/scripts/generer.py" --a-blanc --vague ablanc --usages U3 U2 U1 --taches "$RUN/taches.txt" --batch 4 --etiquette ablanc > "$RUN/ablanc.stdout" 2> "$RUN/ablanc.stderr"; CODE=$?
echo "=== generer --a-blanc -> $CODE"; cat "$RUN/ablanc.stdout"
if [ "$CODE" -ne 0 ]; then tail -25 "$RUN/ablanc.stderr"; fi
exit "$CODE"
