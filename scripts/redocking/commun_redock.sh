# commun_redock.sh — harnais commun du REDOCKING BioEmu2, SOURCE. Il SOURCE le harnais de 6CWG (G2 sha du commit execute + code installe == upstream,
# G13 concurrence echantillonnee, contournement NVML, libnvrtc, HF hors ligne) au lieu de le reecrire, comme generation/scripts/commun_gen.sh. Deux ajouts seulement :
#   - BIOEMU2_MSA_CACHE pointe sous redocking/ (aucune ecriture sous 6CWG/ : le `mkdir -p` de garde_gpu vise un dossier qui existe deja) ;
#   - BIOEMU2_EMBEDS_CACHE est defait.
R=/gpfs/phherman/bioemu2_bench/redocking
. /gpfs/phherman/bioemu2_bench/6CWG/scripts/commun_gpu.sh     # pose B0 UP PYB SHA_UPSTREAM W ; definit garde_gpu / fin_gpu

garde_redock_gpu() {   # $1=unite $2=jobs GPU prevus (moi compris) $3=dossier de run
  garde_gpu "$1" "$2" "$3"
  export BIOEMU2_MSA_CACHE="$R/msa_cache"
  mkdir -p "$BIOEMU2_MSA_CACHE"
  unset BIOEMU2_EMBEDS_CACHE
  export REDOCK_RUN="$3"
  export PYTHONPATH="$R/scripts${PYTHONPATH:+:$PYTHONPATH}"
}

garde_redock_cpu() {   # $1=unite $2=dossier de run — jobs CPU (preparation, MSA, mesures) : G2 + code installe == upstream, ni G13 ni NVML
  local UNITE="$1" RUN="$2"
  mkdir -p "$RUN"
  if touch "$RUN/.essai"; then rm -f "$RUN/.essai"; else echo "ARRET 13 : ecriture impossible sous $RUN"; exit 13; fi
  echo "=== $UNITE · job ${SLURM_JOB_ID:-<hors slurm>} · $(hostname) · $(date -Is) · CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<non pose>}"
  . "$B0/outils/garde_commit.sh"
  garde_commit "$UP" "$SHA_UPSTREAM"; local G2=$?
  if [ "$G2" -ne 0 ]; then echo "=== G2 -> exit $G2"; exit "$G2"; fi
  local SITE
  if ! SITE=$("$PYB" -c "import importlib.util as u,os;print(os.path.dirname(u.find_spec('bioemu2').origin))"); then echo "ARRET 90 : bioemu2 introuvable"; exit 90; fi
  python3 "$B0/outils/arbre_sha.py" "$SITE" --comparer "$UP/src/bioemu2"; local CODE=$?
  if [ "$CODE" -ne 0 ]; then echo "=== code installe != upstream -> exit $CODE"; exit "$CODE"; fi
  export BIOEMU2_MSA_CACHE="$R/msa_cache"
  export TMPDIR="$RUN/.tmp"; mkdir -p "$TMPDIR" "$BIOEMU2_MSA_CACHE"
  export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
  export HF_HOME=/gpfs/phherman/bioemu2/hf_cache HF_HUB_OFFLINE=1
  export REDOCK_RUN="$RUN"
  export PYTHONPATH="$R/scripts${PYTHONPATH:+:$PYTHONPATH}"
}
