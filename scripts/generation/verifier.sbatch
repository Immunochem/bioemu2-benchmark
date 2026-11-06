#!/bin/bash
#SBATCH --job-name=gen_verifier
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Relecture des ARTEFACTS (sha256, images du xtc, physicality.json) contre meta.json et le manifeste.
# --time 01:00:00 : debit MESURE sur le meme geste — le smoke a ecrit 2 x 4 dossiers ; relire ~200 petits dossiers (1 Mo) ~ qq min ; 1 h = plancher de la regle (x 1,35, arrondi a l heure).
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__verifier" "$RUN"
"$PYB" "$G/scripts/verifier_sorties.py" > "$RUN/verif.stdout" 2> "$RUN/verif.stderr"; CODE=$?
echo "=== verifier_sorties -> $CODE"; cat "$RUN/verif.stdout" | cut -c1-600
if [ "$CODE" -ne 0 ]; then tail -20 "$RUN/verif.stderr"; fi
exit "$CODE"
