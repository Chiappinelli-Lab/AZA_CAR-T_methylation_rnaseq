# Every input and output path used by the R scripts. Sourced at the top of each script as
#   source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
# REPO_ROOT is set by run_all.sh; external inputs come from config/paths.yaml.

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (!nzchar(REPO_ROOT) || !file.exists(file.path(REPO_ROOT, "config", "paths.yaml")))
  stop("Set REPO_ROOT to the publication repo root (the folder holding config/paths.yaml)")
CFG <- yaml::read_yaml(file.path(REPO_ROOT, "config", "paths.yaml"))
.mkdir <- function(d) { dir.create(d, showWarnings = FALSE, recursive = TRUE); d }

## ---- external inputs ----
RNASEQ_MANIFEST  <- file.path(REPO_ROOT, "metadata", "haydar_rnaseq_sample_manifest.tsv")
METH_MANIFEST    <- file.path(REPO_ROOT, "metadata", "Methyl Array Sample Manifest EV.xlsx")
IDAT_DIR         <- CFG$data$methylation_idat_dir
GTF_VM32         <- CFG$references$gencode_gtf
CHAIN_MM10_MM39  <- CFG$references$chain_mm10_to_mm39
CCRE_MM10_BED    <- CFG$references$encode_ccre_mm10_bed
MM285_NATIVE_MM39 <- CFG$references$mm285_native_mm39_manifest

## ---- RNA-seq ----
# stage-1 output, unless RNASEQ_COUNTS (env) or rnaseq_counts_override (paths.yaml) points elsewhere
RNASEQ_COUNTS <- Sys.getenv("RNASEQ_COUNTS",
                            if (nzchar(CFG$rnaseq_counts_override %||% "")) CFG$rnaseq_counts_override else
                              file.path(REPO_ROOT, "results", "rnaseq", "output", "results", "counts", "squire_gene_counts.tsv"))
RNA_DE_DIR <- file.path(REPO_ROOT, "results", "rnaseq_de")          # = Haydar_RNAseq/analysis/deseq2/results
RNA_NS_DIR <- file.path(RNA_DE_DIR, "gene_NS")
RNA_S1_DIR <- file.path(RNA_DE_DIR, "gene_S1")

## ---- methylation ----
METH_DIR      <- file.path(REPO_ROOT, "results", "methylation")      # = Haydar_Methylation_Array/results
METH_SHEET    <- file.path(METH_DIR, "metadata", "sample_sheet.csv")
METH_QC_DIR   <- file.path(METH_DIR, "qc")
# folder with the deposited sdfs_QCDPBM.rds / betas_QCDPBM.rds; when set, 01 copies them instead of
# preprocessing the IDATs. METH_PREPROCESSED (env) or methylation_preprocessed_override (paths.yaml)
METH_PREPROCESSED <- Sys.getenv("METH_PREPROCESSED", CFG$methylation_preprocessed_override %||% "")
METH_DM_DIR   <- file.path(METH_DIR, "dm")
DMGSEA_DIR    <- file.path(METH_DIR, "dmgsea")
GSK_DIR       <- file.path(METH_DIR, "gsk_vs_aza")
MM39_REF_DIR  <- file.path(METH_DIR, "reference", "mm39")           # = Haydar_Methylation_Array/reference/mm39
BETAS_RDS     <- file.path(METH_QC_DIR, "betas_QCDPBM.rds")
DMRCATE_RDS   <- file.path(METH_DM_DIR, "dmrcate_aza_b7h3.rds")

## ---- integration and figures ----
INT_DIR  <- file.path(REPO_ROOT, "results", "integration")           # = integration/results
FIG_ROOT <- file.path(REPO_ROOT, "figures")                          # = integration/paper_figures
WORKFLOW <- file.path(REPO_ROOT, "workflow")
