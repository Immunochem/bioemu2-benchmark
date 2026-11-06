#!/bin/bash
#SBATCH --job-name=gen_MD_dynapin
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Reference MD pour la comparaison U2 / U3-lag1000 : DynaPIN v2b sur les 84 piles MD (28 systemes x 3 repliques, stride 5, avec H du champ de force). CPU seul.
# Arguments : [PARALLELE=8]. --time 04:00:00 : debit MESURE par le pilote 6CWG (J4-MD : 3 x 200 images de 6CWG en 9 min 55 s) ; 28 systemes ~ 26 « 6CWG-equivalents »
# (somme des jetons / 390) x 3 repliques x 9,9 min = 12,9 h de conteneur / 8 en parallele = 1,6 h x 1,35 = 2,2 h -> 3 h ; 4 h parce que le cout n est mesure que sur 6CWG.
# Reprenable : dynapin.code == 0 => « deja fait ». Departs etales (juger_piles.sh : 0-8 s aleatoires, une reprise), jamais 12 conteneurs dans la meme seconde.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
P="${1:-8}"
export H_RACINE=md_ref H_USAGE=MD
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__MD_reference_dynapin" "$RUN"
Hd=/gpfs/phherman/bioemu2_bench/hydrogenes
echo "=== piles MD · $(date -Is)"
"$PYB" "$G/scripts/md_piles.py" > "$RUN/mdpiles.stdout" 2> "$RUN/mdpiles.stderr"; C1=$?
echo "=== md_piles -> $C1"; tail -2 "$RUN/mdpiles.stdout" | cut -c1-300
echo "=== DynaPIN ($P en parallele) · $(date -Is) · charge $(cut -d' ' -f1-3 /proc/loadavg)"
bash $Hd/scripts/juger_piles.sh $G/md_ref/piles "$P" > "$RUN/dynapin.stdout" 2> "$RUN/dynapin.stderr"; C2=$?
echo "=== juger_piles -> $C2 · $(date -Is)"; awk '{print $1, $2}' "$RUN/dynapin.stdout" | sort | uniq -c | head
"$PYB" "$G/scripts/h_agreger.py" > "$RUN/agreger.stdout" 2> "$RUN/agreger.stderr"; C3=$?
echo "=== h_agreger -> $C3"; cut -c1-400 "$RUN/agreger.stdout"
[ "$C1" -eq 0 ] && [ "$C2" -eq 0 ] && [ "$C3" -eq 0 ] && exit 0
exit 1
