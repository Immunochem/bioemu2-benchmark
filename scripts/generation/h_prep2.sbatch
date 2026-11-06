#!/bin/bash
#SBATCH --job-name=gen_H_prep2
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Etage H, preparation GENERIQUE (CPU, aucun GPU). Arguments : USAGE RACINE   (ex. U2 hydrogenes_U2).
# Memes etapes que h_prep.sbatch + diagnostic des ponts (cas A/B) + reordonnancement de TOUTES les listes (reference + 7 images propres en tete).
# --time 01:00:00 : debit MESURE — h_prep (job 11122) : 3 500 images + croisement PSF en 34 s ; h_diag (11127) < 1 min ; U2 = 28 unites x 25 images = 700 images (5 fois moins) x 1,35 -> plancher d une heure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
USAGE="${1:?USAGE}"; export H_USAGE="$USAGE" H_RACINE="${2:?RACINE}"
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__H_preparation_$USAGE" "$RUN"
"$PYB" "$G/scripts/h_preparer.py" "$USAGE" > "$RUN/prep.stdout" 2> "$RUN/prep.stderr"; C1=$?
echo "=== h_preparer -> $C1"; tail -2 "$RUN/prep.stdout" | cut -c1-300
[ "$C1" -eq 0 ] || { tail -20 "$RUN/prep.stderr"; exit 31; }
/home/phherman/micromamba/envs/evomd/bin/python "$G/scripts/h_psf_croise.py" > "$RUN/psf.stdout" 2> "$RUN/psf.stderr"; C2=$?
echo "=== croisement PSF -> $C2"; tail -1 "$RUN/psf.stdout" | cut -c1-300
[ "$C2" -eq 0 ] || { tail -20 "$RUN/psf.stderr"; exit 32; }
"$PYB" "$G/scripts/h_diag_ss.py" > "$RUN/diag.stdout" 2> "$RUN/diag.stderr"; C3=$?
echo "=== diagnostic ponts -> $C3"; tail -1 "$RUN/diag.stdout" | cut -c1-300
[ "$C3" -eq 0 ] || { tail -20 "$RUN/diag.stderr"; exit 33; }
"$PYB" "$G/scripts/h_reordonner.py" TOUS > "$RUN/reord.stdout" 2> "$RUN/reord.stderr"; C4=$?
echo "=== reordonnancement -> $C4"; tail -1 "$RUN/reord.stdout" | cut -c1-200
[ "$C4" -eq 0 ] || { tail -20 "$RUN/reord.stderr"; exit 34; }
exit 0
