#!/bin/bash
#SBATCH --job-name=gen_preflight
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# PREFLIGHT (CPU, aucun GPU, aucune requete reseau) : amorce le cache MSA, CLASSE les 28 systemes (dans le cache / a chercher),
# verifie les sha256 des x0 curés, construit System + x0 + MSA pour les systemes de la vague 1 (« a blanc »).
# --time 01:00:00 : débit MESURE sur un lot comparable — curation_suite (tests + PSF + tokens, 28 systemes) = 14 s (sacct 11047) ;
# x 1,35 -> arrondi a l heure = 1 h. Le cache MSA a lire pese 115 Mo (du -sh msa/cache).
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__preflight" "$RUN"
echo "=== amorcage du cache MSA · $(date -Is)"
python3 "$G/scripts/semer_cache.py" > "$RUN/semer.stdout" 2> "$RUN/semer.stderr"; C1=$?
echo "=== semer_cache -> $C1"; cat "$RUN/semer.stdout"; [ "$C1" -eq 0 ] || { tail -20 "$RUN/semer.stderr"; exit 30; }
echo "=== classification (sans reseau) · $(date -Is)"
"$PYB" "$G/scripts/msa_serie.py" --classifier-seulement > "$RUN/classif.stdout" 2> "$RUN/classif.stderr"; C2=$?
echo "=== classification -> $C2"; cat "$RUN/classif.stdout"; [ "$C2" -eq 0 ] || { tail -20 "$RUN/classif.stderr"; exit 31; }
echo "=== a blanc capacite · $(date -Is)"
"$PYB" "$G/scripts/capacite.py" --a-blanc --etiquette preflight "$@" > "$RUN/ablanc.stdout" 2> "$RUN/ablanc.stderr"; C3=$?
echo "=== a blanc -> $C3"; cat "$RUN/ablanc.stdout"; [ "$C3" -eq 0 ] || { tail -20 "$RUN/ablanc.stderr"; exit 32; }
exit 0
