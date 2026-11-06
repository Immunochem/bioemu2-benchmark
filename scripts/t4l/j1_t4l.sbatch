#!/bin/bash
#SBATCH --job-name=t4l_j1
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:50:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/t4l/logs/%x_%j.log
# J1 T4L : MSA ColabFold SERIE de s164 et s162 (cache sous t4l/msa_cache) + smoke GPU de 7 unites (112 structures) : mesure s/structure, VRAM, et
# fait traverser le chemin GPU entier aux trois types d unites (equilibrium, propagate, steer).
# --time 00:50:00 : PLAFOND DECLARE (HYPOTHESE). Aucun debit mesure a 164 residus sur cette pile ; borne haute : 112 structures x 11,94 s (mesure 6CWG, plus gros systeme) x 1,35 = 30 min
# + MSA (2 requetes) + chargement du modele. Le --time de la production sera calcule sur le s/structure MESURE ici.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/6CWG/scripts/commun_gpu.sh
T=/gpfs/phherman/bioemu2_bench/t4l
RUN=$T/runs/j1_$SLURM_JOB_ID
garde_gpu "t4l__j1_smoke" "${PREVUS:?}" "$RUN"
export BIOEMU2_MSA_CACHE=$T/msa_cache; mkdir -p "$BIOEMU2_MSA_CACHE"
echo "=== smoke · $(date -Is) · MSA cache $BIOEMU2_MSA_CACHE"
"$PYB" "$B0/outils/contourner_nvml.py" --mode nvml -- "$T/scripts/bioemu2_t4l.py" --plan "$T/plans/smoke.json" --racine unites_smoke > "$RUN/smoke.stdout" 2> "$RUN/smoke.stderr"; CODE=$?
echo "=== bioemu2_t4l -> $CODE · $(date -Is)"
fin_gpu "$RUN" "$PREVUS"
grep -E "^MSA|^x0 5LZM|^UNITE|^FIN|^ERREUR|^SAUTE|^NON_RECLAME" "$RUN/smoke.stdout" | cut -c1-260
N=$(ls $T/unites_smoke/*/cv.npz | awk 'END{print NR}')
echo "=== artefacts : $N cv.npz sur 7 attendus"
[ "$N" -eq 7 ] || { echo "ARTEFACT INCOMPLET"; tail -25 "$RUN/smoke.stderr"; exit 30; }
exit "$CODE"
