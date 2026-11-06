#!/bin/bash
#SBATCH --job-name=rd_msa
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/redocking/logs/%x_%j.log
# MSA ColabFold des systemes du redocking absents du cache, en SERIE (une requete a la fois, pause 20 s, chien de garde 900 s). Variable : ATTENTE (json), MAX (optionnel : N premiers a chercher).
# Le --time ci-dessus est un PLACEHOLDER : la soumission le remplace. DEBIT MESURE (generation/msa/journal_requetes.jsonl, 18 systemes) : requete non appariee mediane 42,6 s (max 172,8),
# requete appariee mediane 28,5 s (max 78,9) + 2 pauses de 20 s = ~110 s par systeme ; --time = N a chercher x 110 s x 1,35 arrondi a l heure. REPRENABLE : cache ecrit atomiquement.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/redocking/scripts/commun_redock.sh
: "${ATTENTE:?}"
RUN=$R/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_redock_cpu "redock__msa_serie" "$RUN"
# un seul job MSA a la fois sur cette IP : un autre job dont le nom contient « msa » (le mien exclu) -> ARRET 41
AUTRES=$(squeue -h -u "$USER" -o "%i %j %T" | awk -v me="$SLURM_JOB_ID" '$1!=me && tolower($2) ~ /msa/ && $3=="RUNNING"{print}')
if [ -n "$AUTRES" ]; then echo "ARRET 41 : un autre job MSA tourne ($AUTRES)"; exit 41; fi
# amorcage du cache par COPIE (jamais ecrit dans les originaux) : meme cle sha256 = meme requete
for SRC in /gpfs/phherman/bioemu2_bench/generation/msa_cache /gpfs/phherman/bioemu2_bench/msa/cache /gpfs/phherman/bioemu2_bench/6CWG/msa_cache /gpfs/phherman/bioemu2_bench/vram_1030/msa_cache; do
  [ -d "$SRC" ] && cp -n "$SRC"/*.json "$BIOEMU2_MSA_CACHE"/ 2>> "$RUN/amorcage.err"
done
echo "=== cache amorce : $(ls "$BIOEMU2_MSA_CACHE" | awk 'END{print NR}') fichiers"
"$PYB" -B "$R/scripts/redock_msa.py" --attente "$ATTENTE" ${MAX:+--max $MAX} > "$RUN/msa.stdout" 2> "$RUN/msa.stderr"; CODE=$?
echo "=== redock_msa -> $CODE · $(date -Is)"
grep -E "^CLASSIFICATION|^CACHE|^MSA|^REPRISE|^MSA_BILAN" "$RUN/msa.stdout" | cut -c1-300
if [ "$CODE" -ne 0 ]; then echo "--- stderr (30 dernieres lignes)"; tail -30 "$RUN/msa.stderr"; fi
exit "$CODE"
