#!/bin/bash
#SBATCH --job-name=gen_msa
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# MSA ColabFold des systemes absents du cache, en SERIE (une requete a la fois, pause 20 s, chien de garde 900 s).
# --time 03:00:00 : débit MESURE sur la 1re unite du lot MSA precedent (msa_j1, 1B6C : 405,8 s/systeme pour 2 requetes + pauses)
# x 18 systemes x 1,35 = 9 861 s = 2,74 h -> arrondi a l heure = 3 h. Le lot J2 a coute en realite 3 092 s pour 10 systemes
# (309 s/systeme, de 68 s a 655 s) : la formule sur-dimensionne, ce qui est le sens sur. REPRENABLE : le cache est ecrit
# atomiquement, un TIMEOUT perd au plus la requete en vol.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__msa_serie" "$RUN"
echo "=== msa_serie · $(date -Is)"
"$PYB" "$G/scripts/msa_serie.py" --pause 20 --timeout 900 --reprise-pause 120 > "$RUN/msa.stdout" 2> "$RUN/msa.stderr"; CODE=$?
echo "=== msa_serie -> $CODE · $(date -Is)"
cat "$RUN/msa.stdout"
if [ "$CODE" -ne 0 ]; then echo "--- stderr (30 dernieres lignes)"; tail -30 "$RUN/msa.stderr"; fi
exit "$CODE"
