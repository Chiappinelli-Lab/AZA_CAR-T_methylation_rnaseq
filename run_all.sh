#!/usr/bin/env bash
# Pipeline driver: raw FASTQ/IDATs -> results/ -> figures/. Every stage is a SLURM job; later stages wait on
# earlier ones with --dependency=afterok. Logs: logs/<stage>_<jobid>.log.
#
#   ./run_all.sh             stage 1, then stages 2-5 once stage 1 has finished (the normal full run)
#   ./run_all.sh 1           stage 1 only: FASTQ -> SQuIRE gene counts (Snakemake; one to two days)
#   ./run_all.sh 2-5         stages 2-5 only, on the stage-1 counts (or RNASEQ_COUNTS=<matrix> ./run_all.sh 2-5)
#   ./run_all.sh 4-5         any contiguous range works; START_AT=<script> skips earlier scripts of the first stage
#
# Stages: 1 RNA-seq pipeline | 2 RNA-seq DE (deseq2) | 3 methylation (haydar_methylation; independent of 1-2)
#         4 integration (needs 2+3) | 5 figures (needs 4)
# Set up first: config/paths.yaml, and the conda envs in envs/ (see README.md).
source ~/.bashrc >/dev/null 2>&1 || true
set -euo pipefail
REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export REPO_ROOT
cd "$REPO_ROOT"
mkdir -p logs
RANGE=${1:-1-5}
FIRST=${RANGE%-*}; LAST=${RANGE#*-}
want() { (( $1 >= FIRST && $1 <= LAST )); }

cfg() {  # cfg <dotted.key> -> value from config/paths.yaml
  conda run -n data_analysis python -c "import yaml,sys,os; v=yaml.safe_load(open('config/paths.yaml'))
for k in sys.argv[1].split('.'): v=v[k]
print(os.path.expanduser(str(v)))" "$1"
}
ACCOUNT=$(cfg slurm.account); PARTITION=$(cfg slurm.partition)
SB=(sbatch --parsable --account="$ACCOUNT" --partition="$PARTITION" --export=ALL)

## ---- stage 1: Snakemake controller job ----
J1=""
if want 1; then
  WD=$REPO_ROOT/results/rnaseq
  mkdir -p "$WD/logs/snakemake" "$WD/slurm_profile"
  conda run -n data_analysis python workflow/01_rnaseq_pipeline/make_fastq_links.py \
    "$(cfg data.rnaseq_fastq_dir)" metadata/haydar_rnaseq_sample_manifest.tsv "$WD/raw_files/batch1"
  EMAIL=$(cfg slurm.email)
  EXTRA="${EMAIL:+--mail-type=FAIL --mail-user=$EMAIL}${SLURM_EXCLUDE:+ --exclude=$SLURM_EXCLUDE}"
  EXTRA="'${EXTRA# }'"
  sed -e "s|@BINDS@|$(cfg slurm.singularity_binds)|" -e "s|@JOIN@|$(cfg slurm.squire_join_binary)|" \
      -e "s|@ACCOUNT@|$ACCOUNT|" -e "s|@PARTITION@|$PARTITION|" -e "s|\"@SLURM_EXTRA@\"|\"$EXTRA\"|" \
      workflow/01_rnaseq_pipeline/slurm_profile/config.yaml.in > "$WD/slurm_profile/config.yaml"
  CONDA_PREFIX_SM=$(cfg slurm.snakemake_conda_prefix); APPTAINER=$(cfg slurm.apptainer_module)
  mkdir -p "$CONDA_PREFIX_SM"
  J1=$("${SB[@]}" --job-name=pub_stage1 --cpus-per-task=1 --mem=8G --time=7-00:00:00 \
        --output="logs/stage1_%j.log" --wrap "
source ~/.bashrc >/dev/null 2>&1 || true
conda activate snakemake.RNAseq.standard.with.SQuiRE   # built from workflow/01_rnaseq_pipeline/envs/snakemake.yaml
module load '$APPTAINER'      # before set -u: the cluster's module scripts aren't set -u safe
set -euo pipefail
cd '$WD'
snakemake --snakefile '$REPO_ROOT/workflow/01_rnaseq_pipeline/Snakefile' \
  --configfile '$REPO_ROOT/config/paths.yaml' '$REPO_ROOT/workflow/01_rnaseq_pipeline/config.yaml' \
  --config working_dir='$WD/' run_id=\$SLURM_JOB_ID \
  --profile slurm_profile --use-conda --conda-prefix '$CONDA_PREFIX_SM'")
  echo "stage 1 (Snakemake controller): job $J1"
fi

## ---- stages 2-5 ----
dep() { local d=""; for j in "$@"; do [[ -n "$j" ]] && d="${d:+$d:}$j"; done; [[ -n "$d" ]] && echo "--dependency=afterok:$d"; }
stage() {  # stage <n> <cpus> <mem> <time> <dependency jobs...>
  local n=$1 c=$2 m=$3 t=$4; shift 4
  "${SB[@]}" --job-name="pub_stage$n" --cpus-per-task="$c" --mem="$m" --time="$t" $(dep "$@") \
    --output="logs/stage${n}_%j.log" workflow/run_stage.sh "$n"
}
J2=""; J3=""; J4=""
want 2 && { J2=$(stage 2 4 32G 6:00:00 "$J1");  echo "stage 2 (RNA-seq DE): job $J2"; }
want 3 && { J3=$(stage 3 8 48G 8:00:00);        echo "stage 3 (methylation): job $J3"; }
want 4 && { J4=$(stage 4 2 32G 4:00:00 "$J2" "$J3"); echo "stage 4 (integration): job $J4"; }
want 5 && { J5=$(stage 5 2 32G 4:00:00 "$J4");  echo "stage 5 (figures): job $J5"; }
exit 0
