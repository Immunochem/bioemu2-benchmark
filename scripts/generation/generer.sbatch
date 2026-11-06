#!/bin/bash
#SBATCH --job-name=gen_prod
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Production U1/U2/U3. Arguments : VAGUE "USAGES" BATCH PREVUS TACHES_FILE [BUDGET_S] [COUT_JSON] [BATCH_JSON] ["ARGS SUPPLEMENTAIRES"]
#   PREVUS = nombre de workers GPU prevus, moi compris (G13). Le --time est PASSE EN LIGNE DE COMMANDE (placeholder ici) :
#   il vient du debit MESURE x 1,35, arrondi a l heure — voir la decision de la vague.
set -uo pipefail
VAGUE="$1"; USAGES="$2"; BATCH="$3"; PREVUS="$4"; TACHES="$5"; BUDGET="${6:-0}"; COUT="${7:-}"; BJSON="${8:-}"; SUPP="${9:-}"
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen "gen__${VAGUE}__${USAGES// /+}" "$PREVUS" "$RUN"
EXTRA=""
[ -n "$COUT" ] && EXTRA="$EXTRA --cout-json $COUT"
[ -n "$BJSON" ] && EXTRA="$EXTRA --batch-json $BJSON"
echo "=== generer · vague $VAGUE · usages $USAGES · batch $BATCH · taches $TACHES · budget ${BUDGET}s · prevus $PREVUS · $(date -Is)"
"$PYB" "$B0/outils/contourner_nvml.py" --mode nvml -- "$G/scripts/generer.py" --vague "$VAGUE" --usages $USAGES --taches "$TACHES" --batch "$BATCH" --budget-s "$BUDGET" --etiquette "c${PREVUS}" $EXTRA $SUPP > "$RUN/generer.stdout" 2> "$RUN/generer.stderr"; CODE=$?
echo "=== generer -> $CODE · $(date -Is)"
fin_gpu "$RUN" "$PREVUS"
"$PYB" "$G/scripts/gen_lib.py" fin_job --run "$RUN" --code "$CODE" --unite "gen__${VAGUE}"
grep -E "^MODELE|^RECLAME|^TERMINE|^NON_RECLAME|^UNITE|^OOM|^ECHEC|^REJET|^SAUTE|^EMBEDDING|^FIN" "$RUN/generer.stdout" | tail -n 120
if [ "$CODE" -ne 0 ]; then echo "--- stderr (25 dernieres lignes)"; tail -25 "$RUN/generer.stderr"; fi
exit "$CODE"
