#!/bin/bash
#SBATCH --job-name=gen_ponts
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Ponts disulfure de reference (PSF + xtc de la MD) pour les 28 systemes. CPU seul, aucun --gres (partition gpu comme les jobs CPU du
# pilote : un job GPU en attente dans `shared` bloque cette partition, mesure du 2026-11-05).
# --time 01:00:00 : debit MESURE sur le lot comparable le plus proche — curation_audit (lecture de toutes les entrees des 28 systemes)
# = 3 min 00 s (sacct 11041) x 1,35 -> arrondi a l heure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__ponts_reference" "$RUN"
"$PYB" "$G/scripts/ponts_reference.py" > "$RUN/ponts.stdout" 2> "$RUN/ponts.stderr"; CODE=$?
echo "=== ponts_reference -> $CODE"; cut -c1-330 "$RUN/ponts.stdout"
if [ "$CODE" -ne 0 ]; then tail -25 "$RUN/ponts.stderr"; fi
exit "$CODE"
