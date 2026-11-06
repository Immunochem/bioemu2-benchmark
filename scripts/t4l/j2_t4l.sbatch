#!/bin/bash
#SBATCH --job-name=t4l_j2
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/t4l/logs/%x_%j.log
# J2 T4L (production). Variables : PLAN (nom du JSON sous t4l/plans), ATTENDU (nombre de cv.npz attendus en tout), PREVUS (calculs GPU prevus, moi compris).
# --time : PASSE A LA SOUMISSION, = s/structure MESURE par le smoke J1 (job 11159) x structures x 1,35, + 1 min de chargement. Le 01:00:00 ci-dessus n est qu un placeholder.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/6CWG/scripts/commun_gpu.sh
T=/gpfs/phherman/bioemu2_bench/t4l
RUN=$T/runs/j2_${PLAN:?}_$SLURM_JOB_ID
garde_gpu "t4l__j2_${PLAN}" "${PREVUS:?}" "$RUN"
export BIOEMU2_MSA_CACHE=$T/msa_cache
echo "=== production $PLAN · $(date -Is)"
"$PYB" "$B0/outils/contourner_nvml.py" --mode nvml -- "$T/scripts/bioemu2_t4l.py" --plan "$T/plans/$PLAN.json" --racine unites > "$RUN/prod.stdout" 2> "$RUN/prod.stderr"; CODE=$?
echo "=== bioemu2_t4l -> $CODE · $(date -Is)"
fin_gpu "$RUN" "$PREVUS"
grep -E "^UNITE|^FIN|^ERREUR|^SAUTE|^NON_RECLAME" "$RUN/prod.stdout" | cut -c1-280
N=$(ls $T/unites/*/cv.npz | awk 'END{print NR}')
echo "=== artefacts : $N cv.npz sur ${ATTENDU:?} attendus"
[ "$N" -ge "$ATTENDU" ] || { echo "ARTEFACT INCOMPLET"; tail -25 "$RUN/prod.stderr"; exit 30; }
exit "$CODE"
