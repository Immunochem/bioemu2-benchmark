#!/bin/bash
#SBATCH --job-name=rd_gpu
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/redocking/logs/%x_%j.log
# REDOCKING BioEmu2 — UNE slice GPU. Variables : TACHES (json), PREVUS (jobs GPU prevus, moi compris : G13), ETIQ (dossier de sortie stable : reprise), LAGS (defaut 50), BATCH (defaut 4),
# IDS (optionnel : sous-ensemble d identifiants), KREP (images par condition, defaut 3), BUDGET (secondes : une tache n est lancee que si 1,35 x son estimation tient). Le --time ci-dessus est un PLACEHOLDER : la soumission le remplace (--time = debit MESURE x 1,35 arrondi a l heure).
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/redocking/scripts/commun_redock.sh
: "${TACHES:?}" "${PREVUS:?}" "${ETIQ:?}"
RUN=$R/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_redock_gpu "redock__gpu_${ETIQ}" "$PREVUS" "$RUN"
OUT=$R/sorties/$ETIQ
mkdir -p "$OUT"
echo "=== redock_gpu · etiquette $ETIQ · taches $TACHES · lags ${LAGS:-50} · batch ${BATCH:-4} · budget ${BUDGET:-<aucun>} s · $(date -Is)"
"$PYB" "$B0/outils/contourner_nvml.py" --mode nvml -- "$R/scripts/redock_gpu.py" --taches "$TACHES" --sortie "$OUT" --lags ${LAGS:-50} --batch ${BATCH:-4} ${BUDGET:+--budget-s $BUDGET} ${IDS:+--ids $IDS} ${KREP:+--K $KREP} > "$RUN/gpu.stdout" 2> "$RUN/gpu.stderr"; CODE=$?
echo "=== redock_gpu -> $CODE · $(date -Is)"
fin_gpu "$RUN" "$PREVUS"
grep -E "^UNITE|^FIN|^REJET|^OOM|^ECHEC|^SAUTE|^NON LANCE" "$RUN/gpu.stdout"
N=$(ls "$OUT"/*/lag*/samples.xtc 2>/dev/null | awk 'END{print NR}')
echo "=== artefacts : $N samples.xtc sous $OUT"
if [ "$CODE" -ne 0 ]; then echo "--- stderr (30 dernieres lignes)"; tail -30 "$RUN/gpu.stderr"; fi
[ "$N" -ge 1 ] || { echo "ARTEFACTS ABSENTS"; exit 30; }
exit "$CODE"
