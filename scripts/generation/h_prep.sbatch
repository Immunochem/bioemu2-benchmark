#!/bin/bash
#SBATCH --job-name=gen_H_prep
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Etage H, preparation (CPU, aucun GPU) : 140 unites U3 -> 3 500 images PDB d atomes lourds + 28 references, puis croisement PSF / tableau de protonation / ponts.
# --time 01:00:00 : debit MESURE sur le geste le plus proche — le job ponts (11091) a relu 28 x 3 xtc + PSF en 10,6 s ; ecrire 3 528 PDB de 3 000-6 300 atomes
# (~0,3 s/image par save_pdb mdtraj, hypothese) ~ 20 min au pire x 1,35 -> arrondi a l heure.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__H_preparation" "$RUN"
echo "=== h_preparer · $(date -Is)"
"$PYB" "$G/scripts/h_preparer.py" U3 > "$RUN/prep.stdout" 2> "$RUN/prep.stderr"; C1=$?
echo "=== h_preparer -> $C1"; tail -4 "$RUN/prep.stdout" | cut -c1-300
[ "$C1" -eq 0 ] || { tail -20 "$RUN/prep.stderr"; exit 31; }
echo "=== croisement PSF · $(date -Is)"
/home/phherman/micromamba/envs/evomd/bin/python "$G/scripts/h_psf_croise.py" > "$RUN/psf.stdout" 2> "$RUN/psf.stderr"; C2=$?
echo "=== croisement PSF -> $C2"; cut -c1-260 "$RUN/psf.stdout"
[ "$C2" -eq 0 ] || { tail -20 "$RUN/psf.stderr"; exit 32; }
exit 0
