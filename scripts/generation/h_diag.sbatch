#!/bin/bash
#SBATCH --job-name=gen_H_diag
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Diagnostic CPU : images U3 dont un pont declare est rompu (> 3,0 A) ou dont deux cysteines libres sont a < 3,0 A — cas ou OpenMM infere une liaison S–S contraire au PSF.
# --time 01:00:00 : debit MESURE sur le geste voisin — 28 systemes x 5 xtc relus par ponts/rapport (job de verification 11121 : 144 dossiers en < 1 min) x 1,35 -> arrondi a l heure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__H_diag_ss" "$RUN"
"$PYB" "$G/scripts/h_diag_ss.py" > "$RUN/diag.stdout" 2> "$RUN/diag.stderr"; CODE=$?
echo "=== h_diag_ss -> $CODE"; cut -c1-300 "$RUN/diag.stdout"
if [ "$CODE" -ne 0 ]; then tail -20 "$RUN/diag.stderr"; fi
exit "$CODE"
