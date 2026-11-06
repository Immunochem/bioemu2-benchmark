#!/bin/bash
#SBATCH --job-name=rd_plancher
#SBATCH --partition=shared
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/redocking/logs/%x_%j.log
# CPU. (1) CONTROLE POSITIF DE LA MESURE (CONTROLE=1) : mes fnat/iRMSD/lRMSD/DockQ contre DynaPIN sur les piles du pilote 6CWG. (2) PLANCHER DU HASARD + verification des leurres sur les strates $STRATES.
# Le --time ci-dessus est un PLACEHOLDER : la soumission le remplace (--time = debit MESURE x 1,35 arrondi a l heure ; le 1er lot A_corpus MESURE le debit pour les suivants).
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/redocking/scripts/commun_redock.sh
RUN=$R/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_redock_cpu "redock__plancher" "$RUN"
if [ "${CONTROLE:-0}" = "1" ]; then
  "$PYB" -B "$R/scripts/controle_mesure.py" > "$RUN/controle.stdout" 2> "$RUN/controle.stderr"; CC=$?
  echo "=== controle_mesure -> $CC"; tail -8 "$RUN/controle.stdout"
  [ "$CC" -eq 0 ] || { echo "🔴 la mesure ne retrouve pas DynaPIN : on s arrete ici"; tail -15 "$RUN/controle.stderr"; exit "$CC"; }
fi
T0=$(date +%s)
"$PYB" -B "$R/scripts/plancher.py" ${STRATES:?} > "$RUN/plancher.stdout" 2> "$RUN/plancher.stderr"; CP=$?
echo "=== plancher -> $CP · $(( $(date +%s) - T0 )) s · strates ${STRATES}"
grep -E "^BILAN_PLANCHER|^ERREUR" "$RUN/plancher.stdout"
cut -c1-330 "$RUN/plancher.stdout" | head -80
[ "$CP" -eq 0 ] || { tail -20 "$RUN/plancher.stderr"; exit "$CP"; }
N=$(ls "$R"/resultats/plancher/*.json 2>/dev/null | awk 'END{print NR}')
echo "=== artefacts : $N fichiers plancher"
[ "$N" -ge 1 ] || { echo "ARTEFACT ABSENT"; exit 30; }
exit 0
