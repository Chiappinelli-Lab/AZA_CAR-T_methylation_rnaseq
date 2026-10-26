#!/usr/bin/env bash
# Runs one downstream stage (2-5) in order, each script in its conda env. Called by run_all.sh, inside a
# SLURM job; can also be run directly:  REPO_ROOT=<repo> workflow/run_stage.sh <2|3|4|5> [START_AT=<script>]
# Stage 1 (FASTQ -> counts, Snakemake) is run by run_all.sh itself.
source ~/.bashrc >/dev/null 2>&1 || true      # conda; ~/.bashrc isn't set -u safe
set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
export REPO_ROOT
STAGE=$1
W=$REPO_ROOT/workflow
STARTED=0; [[ -z "${START_AT:-}" ]] && STARTED=1
go() { [[ $STARTED == 1 ]] && return 0; [[ "$(basename "$1")" == "$START_AT" ]] && STARTED=1; [[ $STARTED == 1 ]]; }
run_r()  { go "$2" || return 0; echo "### [$1] $(basename "$2") ${*:3}  $(date '+%F %T')"
           (cd "$REPO_ROOT" && conda run -n "$1" --no-capture-output Rscript "$2" "${@:3}" 2>&1 | grep -v Fontconfig); }
run_py() { go "$1" || return 0; echo "### [data_analysis] $(basename "$1") ${*:2}  $(date '+%F %T')"
           (cd "$REPO_ROOT" && conda run -n data_analysis --no-capture-output python "$1" "${@:2}"); }

case $STAGE in
  2)  # RNA-seq differential expression (env deseq2)
      echo "### counts: ${RNASEQ_COUNTS:-<stage-1 output or paths.yaml override>}"
      for s in run_deseq2_gene_NS.R run_deseq2_gene_S1.R run_gsea_gene_NS.R; do run_r deseq2 "$W/02_rnaseq_de/$s"; done ;;
  3)  # methylation array (env haydar_methylation); 11 uses SLURM_CPUS_PER_TASK cores
      for s in 00_build_sample_sheet.R 00b_build_mm39_annotation.R 01_load_idats_qc.R 06_dmrcate.R \
               06b_gene_background.R 11_dmgsea.R 13_gsk_vs_aza.R; do
        run_r haydar_methylation "$W/03_methylation/$s"
      done ;;
  4)  # methylation x RNA: shared GO terms, gene-level slopes, candidate genes (needs stages 2 and 3)
      run_r deseq2             "$W/03_methylation/12_dmgsea_plots.R"
      run_r deseq2             "$W/04_integration/01_rna_donor_slopes.R"
      run_r haydar_methylation "$W/04_integration/02_methylation_donor_slopes.R"
      run_r deseq2             "$W/04_integration/03_correlate_and_plot.R"
      run_r deseq2             "$W/04_integration/10_candidate_gene_rna_data.R"
      run_r haydar_methylation "$W/04_integration/11_candidate_gene_meth_data.R"
      run_r haydar_methylation "$W/04_integration/12_candidate_gene_heatmaps.R" ;;
  5)  # paper figures: data prep (R) then plotting (Python, env data_analysis)
      F=$W/05_figures
      run_r deseq2             "$F/fig2_01_prep_rna.R"
      run_r haydar_methylation "$F/fig2_02_prep_meth.R"
      run_r deseq2             "$F/fig3_01_prep_rna.R"
      run_r haydar_methylation "$F/fig3_02_prep_meth.R"
      run_r haydar_methylation "$F/fig3_04_genomewide_fdr.R"
      run_r haydar_methylation "$F/figS_gsk_02_prep.R"
      run_py "$F/fig2_03_plot.py"                         # working files (figures/fig2_rna_methylation/)
      run_py "$F/fig2_03_plot.py" --font-bump 1 --tag v3  # the shared release (release_v3/)
      run_py "$F/fig3_03_plot.py"
      run_py "$F/fig3_03_plot.py" --absolute
      run_py "$F/fig3_03_plot.py" --bars
      run_py "$F/fig3_03_plot.py" --absolute --bars
      run_py "$F/figS_gsk_03_plot.py" ;;
  *)  echo "unknown stage: $STAGE" >&2; exit 1 ;;
esac
echo "### stage $STAGE done  $(date '+%F %T')"
