#!/bin/bash
#SBATCH --job-name=rd_prep
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/redocking/logs/%x_%j.log
# Preparation CPU des x0 cristallins des strates A et B (extraction -> curation). CPU, aucune ecriture hors redocking/.
# --time 90 min : AUCUNE mesure de debit pour cette etape. Seule mesure : 22 cristaux de A_corpus traites en 8 s (job 11145, extraction + curation, hors telechargement).
# Estimation : ~300 telechargements x (pause 1 s + ~1 s de lecture) + curation ~1 s = ~15 min, x 1,35 puis marge pour l API RCSB (metadonnees A_apparie) : plafond DECLARE, pas mesure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/redocking/scripts/commun_redock.sh
RUN=$R/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_redock_cpu "redock__prep_AB" "$RUN"
"$PYB" -B "$R/scripts/preparer_AB.py" > "$RUN/prep.stdout" 2> "$RUN/prep.stderr"; CODE=$?
echo "=== preparer_AB -> $CODE · $(date -Is)"
grep -E "^BILAN_|^FIN|^B : |^A_apparie : " "$RUN/prep.stdout"
if [ "$CODE" -ne 0 ]; then echo "--- stderr (30 dernieres lignes)"; tail -30 "$RUN/prep.stderr"; exit "$CODE"; fi
[ -s "$R/strate_A_corpus/manifeste.json" ] && [ -s "$R/strate_B/manifeste.json" ] && [ -s "$R/strate_A_apparie/manifeste.json" ] || { echo "ARTEFACT ABSENT"; exit 30; }
exit 0
