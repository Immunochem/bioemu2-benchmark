#!/bin/bash
#SBATCH --job-name=gen_MD_reprise
#SBATCH --partition=gpu
#SBATCH --qos=owner
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/phherman/bioemu2_bench/generation/logs/%x_%j.log
# Reprise des 9 piles MD (3RJQ, 4M3K, 5HGG x 3 repliques) mortes dans DSSP (« FileNotFoundError ... .dssp ») : ces trois systemes portent des CODES D INSERTION, mdtraj ecrit des resSeq dupliques.
# Remede : numerotation ORDINALE 1..n par chaine (contrat d entree, regle 6) pour CES piles. Les 75 autres gardent la numerotation d auteur : leurs metriques (iRMSD, fnat, RMSD) n en dependent pas, leurs
# ETIQUETTES de contacts (ALA114) different de celles des piles BioEmu2 (ordinales) — a traduire par x0/<PDB>__img0_rep1.json:map_residus si on compare des contacts residu par residu.
# --time 02:00:00 : debit MESURE sur ce lot — 65 piles ok en ~80 min a 8 en parallele (journal 11140) soit ~10 min/pile ; 9 piles (3 systemes de 370-442 jetons) a 8 en parallele ~ 1 vague de ~20 min x 1,35 -> 1 h -> 2 h par prudence.
set -uo pipefail
. /gpfs/phherman/bioemu2_bench/generation/scripts/commun_gen.sh
export H_RACINE=md_ref H_USAGE=MD
RUN=$G/runs/${SLURM_JOB_NAME}_$SLURM_JOB_ID
garde_gen_cpu "gen__MD_reference_reprise" "$RUN"
Hd=/gpfs/phherman/bioemu2_bench/hydrogenes
for s in 3RJQ 4M3K 5HGG; do for r in 1 2 3; do rm -rf "$G/md_ref/piles/${s}__MD__rep${r}"; done; done
"$PYB" "$G/scripts/md_piles.py" --renumeroter 3RJQ 4M3K 5HGG > "$RUN/mdpiles.stdout" 2> "$RUN/mdpiles.stderr"; C1=$?
echo "=== md_piles (reprise) -> $C1"; tail -2 "$RUN/mdpiles.stdout" | cut -c1-300
bash $Hd/scripts/juger_piles.sh $G/md_ref/piles 9 > "$RUN/dynapin.stdout" 2> "$RUN/dynapin.stderr"; C2=$?
echo "=== juger_piles -> $C2"; awk '{print $1, $2}' "$RUN/dynapin.stdout" | sort | uniq -c | head
"$PYB" "$G/scripts/h_agreger.py" > "$RUN/agreger.stdout" 2> "$RUN/agreger.stderr"; C3=$?
echo "=== h_agreger -> $C3"; cut -c1-400 "$RUN/agreger.stdout"
[ "$C1" -eq 0 ] && [ "$C2" -eq 0 ] && [ "$C3" -eq 0 ] && exit 0
exit 1
