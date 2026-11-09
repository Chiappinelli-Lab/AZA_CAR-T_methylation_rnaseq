# Load all 24 IDATs and run basic per-sample QC for the Haydar mouse
# methylation array project (Infinium Mouse Methylation BeadChip / MM285).
#
# Preprocessing choices:
#   - prep chain: QCDPBM (qualityMask -> infer Infinium-I channel ->
#     nonlinear dye-bias correction -> pOOBAH detection masking -> noob
#     background subtraction -> matchDesign). No strain inference (T)
#     added since all donors are already known C57BL/6J.
#   - qualityMask uses a CUSTOM mask (ref_issue + design_issue), not
#     sesame's "recommended" set, because "recommended" also includes
#     the `nonunique` mask (probes with mapQ==0, i.e. non-unique genomic
#     mapping) which would systematically remove repeat/TE-derived
#     probes.
#   - pOOBAH detection p-value threshold: 0.05 (sesame default).
#   - background subtraction: noob (sesame default), not scrub.
#   - min_beads: 1 (sesame default, minimal bead-count filtering).
#   - matchDesign (M) added after seeing real Infinium-I vs II probe-design
#     bias in sesameQC_plotBetaByDesign() on this data (left-shifted
#     Infinium-I Green methylated-mode peak). min_dbeta = 0.3 (sesame
#     default).
#   - No global between-sample quantile normalization of betas, and no
#     statistical chip/batch correction (e.g. ComBat): donor and Sentrix
#     chip are almost completely confounded (only EV07 spans both chips),
#     so a batch correction would just re-absorb donor variation. Donor
#     stays as the blocking/random effect in the DM model instead.
#   - matchDesign (M) is stochastic in sesame 1.28.1. set.seed(1) before each sample's prepSesame()
#     makes reruns identical. The paper's results come
#     from an earlier unseeded run, so a fresh run here does NOT reproduce them digit for digit: to
#     reproduce the paper, set METH_PREPROCESSED (env) or methylation_preprocessed_override
#     (config/paths.yaml) to the folder holding the deposited sdfs_QCDPBM.rds / betas_QCDPBM.rds,
#     and this script copies them in instead of preprocessing.
#   - Betas: openSesame(sdfs, func = getBetas) applies its default prep="QCDPB" again on top of the
#     QCDPBM SigDFs, so the betas used downstream are QCDPBM followed by QCDPB. This is how the
#     paper's betas were made (deterministic) and is kept as is.
#   - Sample-level pass/fail cutoffs are deliberately NOT applied here --
#     this script only computes and reports per-sample QC stats; exclusion
#     thresholds will be set after reviewing the actual distribution.

suppressMessages({
  library(sesame)
  library(dplyr)
  library(purrr)
  library(tibble)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
sample_sheet <- read.csv(METH_SHEET, stringsAsFactors = FALSE)

out_dir <- METH_QC_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (nzchar(METH_PREPROCESSED)) {
  qc_files <- c("sdfs_QCDPBM.rds", "betas_QCDPBM.rds", "per_sample_qc_stats.csv")
  src <- file.path(METH_PREPROCESSED, qc_files)
  stopifnot(all(file.exists(src[1:2])))
  file.copy(src[file.exists(src)], out_dir, overwrite = TRUE)
  writeLines(paste("Copied from", METH_PREPROCESSED, "on", Sys.time()), file.path(out_dir, "SOURCE.txt"))
  cat("Using deposited preprocessed SigDFs/betas from", METH_PREPROCESSED, "\n")
  quit(save = "no")
}

custom_mask_names <- c("ref_issue", "design_issue")
prep_args <- list(Q = list(mask_names = custom_mask_names))

cat("Loading", nrow(sample_sheet), "SigDFs and applying prep='QCDPBM' ...\n")
sdfs <- setNames(
  lapply(sample_sheet$Basename, function(bn) {
    sdf <- readIDATpair(bn)
    set.seed(1)  # matchDesign (M) is stochastic
    prepSesame(sdf, prep = "QCDPBM", prep_args = prep_args)
  }),
  sample_sheet$Code
)

cat("Computing per-sample QC stats ...\n")
qc_list <- lapply(sdfs, sesameQC_calcStats)
qc_df <- map_dfr(qc_list, ~ as.data.frame(as.list(sesameQC_getStats(.x))), .id = "Code")

qc_df <- left_join(sample_sheet %>% select(Code, Donor, CAR, Dose, Drug), qc_df, by = "Code")

write.csv(qc_df, file.path(out_dir, "per_sample_qc_stats.csv"), row.names = FALSE)

saveRDS(sdfs, file.path(out_dir, "sdfs_QCDPBM.rds"))

cat("Computing betas ...\n")
betas <- openSesame(sdfs, func = getBetas)
saveRDS(betas, file.path(out_dir, "betas_QCDPBM.rds"))

cat("\n=== Per-sample QC summary ===\n")
print(as.data.frame(qc_df), row.names = FALSE)

cat("\nWrote", nrow(qc_df), "rows to", file.path(out_dir, "per_sample_qc_stats.csv"), "\n")
cat("Saved processed SigDF list to", file.path(out_dir, "sdfs_QCDPBM.rds"), "\n")
cat("Saved betas to", file.path(out_dir, "betas_QCDPBM.rds"), "\n")
