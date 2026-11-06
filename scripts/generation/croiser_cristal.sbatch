#!/bin/bash
#SBATCH --job-name=gen_cristal
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Croisement cristal (mmCIF RCSB struct_conn) / PSF de la MD / cysteines libres proches, 28 systemes. CPU + reseau sortant (28 telechargements en serie, pause 1 s).
# --time 01:00:00 : debit MESURE — le corpus a telecharge 210 fichiers DynaBench en 0,3 s chacun (fetch_structures.jsonl) ; les 28 mmCIF RCSB pesent ~0,1-1 Mo ; lecture des
# xtc 3 repliques x 14 systemes a pont ~ 10 s (job ponts 11091, 10,6 s) -> << 1 h ; plancher de la regle (x 1,35 arrondi a l heure).
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__croisement_cristal" "$RUN"
"$PYB" "$G/scripts/croiser_cristal.py" > "$RUN/cristal.stdout" 2> "$RUN/cristal.stderr"; CODE=$?
echo "=== croiser_cristal -> $CODE"; cut -c1-420 "$RUN/cristal.stdout"
if [ "$CODE" -ne 0 ]; then tail -20 "$RUN/cristal.stderr"; fi
exit "$CODE"
