#!/bin/bash
#SBATCH --job-name=gen_cap
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Courbe de capacite + echelle de batch + concurrence. Arguments : ETIQ PREVUS "SYSTEMES" "BATCHES" [BARRIERE_N BARRIERE_DIR]
#   ETIQ   = regime de concurrence PREVU (solo | c4), ecrit dans chaque ligne du manifeste
#   PREVUS = nombre de jobs GPU prevus, moi compris (G13 : le job REFUSE de demarrer s il y en a davantage)
# Le --time est PASSE EN LIGNE DE COMMANDE (placeholder ici) : voir la decision du lot.
set -uo pipefail
ETIQ="$1"; PREVUS="$2"; SYSTEMES="$3"; BATCHES="$4"; BARN="${5:-0}"; BARD="${6:-}"
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen "gen__capacite__${ETIQ}" "$PREVUS" "$RUN"
echo "=== capacite · $SYSTEMES · batches $BATCHES · etiquette $ETIQ · barriere $BARN · $(date -Is)"
"$PYB" "$B0/outils/contourner_nvml.py" --mode nvml -- "$G/scripts/capacite.py" $SYSTEMES --batches $BATCHES --etiquette "$ETIQ" --barriere-n "$BARN" --barriere-dir "$BARD" > "$RUN/capacite.stdout" 2> "$RUN/capacite.stderr"; CODE=$?
echo "=== capacite -> $CODE · $(date -Is)"
fin_gpu "$RUN" "$PREVUS"
"$PYB" "$G/scripts/gen_lib.py" fin_job --run "$RUN" --code "$CODE" --unite "gen__capacite__${ETIQ}"
grep -E "^PREP|^REJET|^MODELE|^BARRIERE|premier appel|bs[0-9]+ |^FIN" "$RUN/capacite.stdout"
if [ ! -s "$G/capacite/job${SLURM_JOB_ID}__${ETIQ}.json" ]; then echo "ARTEFACT ABSENT : capacite/job${SLURM_JOB_ID}__${ETIQ}.json"; tail -25 "$RUN/capacite.stderr"; exit 30; fi
exit "$CODE"
