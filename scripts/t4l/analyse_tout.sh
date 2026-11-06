#!/bin/bash
# analyse_tout.sh <racine> : verification en ligne/hors ligne, analyse, RMSD au cristal, RMSD aux references fermees, tableaux. CPU seul. Codes de retour captures AVANT tout tube.
T=/gpfs/phherman/bioemu2_bench/t4l; PY=/gpfs/phherman/bioemu2/venv/bin/python
export CUDA_VISIBLE_DEVICES= PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
R=${1:?racine}; cd $T
UN=$(ls $R | grep -v bilan.json | tr '\n' ' ')
echo "unites : $UN"
$PY -I scripts/verifier_enligne.py $R > logs/verif_$R.out 2> logs/verif_$R.err; echo "verifier_enligne -> $?"
$PY -I scripts/analyser_t4l.py $R > logs/analyse_$R.out 2> logs/analyse_$R.err; echo "analyser -> $?"
$PY -I scripts/rmsd_cristal.py $R > logs/rmsd_cristal_$R.out 2> logs/rmsd_cristal_$R.err; echo "rmsd_cristal -> $?"
$PY -I scripts/rmsd_ferme.py $R $UN > logs/rmsd_ferme_$R.out 2> logs/rmsd_ferme_$R.err; echo "rmsd_ferme -> $?"
if [ "$R" = "unites" ]; then $PY -I scripts/tables_rapport.py > logs/tables.out 2> logs/tables.err; echo "tables -> $?"; fi
echo "FIN analyse_tout $R · $(date -Is)"
