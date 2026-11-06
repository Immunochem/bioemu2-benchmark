#!/bin/bash
#SBATCH --job-name=gen_H_dynapin
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Etage H, piles + DynaPIN + agregation (CPU, aucun GPU). Arguments : [PARALLELE=6] [SYSTEMES... pour une reprise].
# --time 02:00:00 : debit MESURE par l agent hydrogenes, job 11100 : 20 piles de 6CWG (26 modeles) en 249 s a 6 en parallele = 12,5 s/pile ; 140 piles de 28 systemes
# de 126 a 777 jetons (moyenne ~0,93 x 6CWG) = 1 630 s = 27 min x 1,35 = 37 min -> arrondi a l heure = 1 h ; DOUBLE a 2 h parce que le cout de DynaPIN en fonction
# de la taille n est mesure que sur 6CWG (390 jetons) : les piles de 777 jetons peuvent couter beaucoup plus que la moyenne. Reprenable : dynapin.code == 0 => « deja fait ».
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
P="${1:-6}"; shift
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__H_piles_dynapin" "$RUN"
Hd=/gpfs/phherman/bioemu2_bench/hydrogenes
echo "=== piles + gardes 40/41/42/43 · $(date -Is)"
/home/phherman/micromamba/envs/evomd/bin/python "$G/scripts/h_piles.py" "$@" > "$RUN/piles.stdout" 2> "$RUN/piles.stderr"; C1=$?
echo "=== h_piles -> $C1"; tail -3 "$RUN/piles.stdout" | cut -c1-300
grep -E "echec" "$RUN/piles.stdout" | head -10 | cut -c1-300
echo "=== DynaPIN ($P en parallele, departs etales) · $(date -Is) · charge $(cut -d' ' -f1-3 /proc/loadavg)"
bash $Hd/scripts/juger_piles.sh $G/${H_RACINE:-hydrogenes}/piles "$P" > "$RUN/dynapin.stdout" 2> "$RUN/dynapin.stderr"; C2=$?
echo "=== juger_piles -> $C2 · $(date -Is)"
awk '{print $1, $2}' "$RUN/dynapin.stdout" | sort | uniq -c | head
echo "=== agregation · $(date -Is)"
"$PYB" "$G/scripts/h_agreger.py" > "$RUN/agreger.stdout" 2> "$RUN/agreger.stderr"; C3=$?
echo "=== h_agreger -> $C3"; cut -c1-500 "$RUN/agreger.stdout"
if [ "$C3" -ne 0 ]; then tail -10 "$RUN/agreger.stderr"; fi
[ "$C1" -eq 0 ] && [ "$C2" -eq 0 ] && [ "$C3" -eq 0 ] && exit 0
exit 1
