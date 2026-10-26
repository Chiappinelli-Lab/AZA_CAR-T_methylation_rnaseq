#!/usr/bin/env Rscript
# Per-sample + per-dose methylation data for the collaborators' candidate gene sets,
# feeding the methylation panel of the per-pathway 3-panel heatmaps.
#
# Uses the same per-gene composite as 02_methylation_donor_slopes.R -- median
# beta over the gene's promoter (1.5kb upstream) + gene-body probes, dmGsea
# mapping (meth_gene_mapping.R) -- but keeps the per-sample matrix instead of
# collapsing straight to a per-donor slope, since the heatmap script needs both
# the full per-sample matrix and a per-dose (averaged across donors) trajectory.
# Switched 2026-09-24 from "mean region-average beta over all DMRs within 5kb";
# every candidate gene with a mapped probe now has data, DMR or not.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- file.path(INT_DIR, "candidate_genes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(WORKFLOW, "04_integration", "candidate_gene_sets.R"))
source(file.path(WORKFLOW, "04_integration", "meth_gene_mapping.R"))
all_genes <- sort(unique(unlist(candidate_gene_sets)))

## ---- 1. Probe -> gene map (promoter + body, dmGsea) ----
gmap <- load_meth_gene_probe_map()
found <- intersect(all_genes, unique(gmap$symbol))
missing <- setdiff(all_genes, found)
if (length(missing) > 0) cat("Candidate genes with NO promoter/body probe (no methylation data possible):\n  ", paste(missing, collapse = ", "), "\n")

## ---- 2. Load beta matrix, subset to AZA + B7H3 CAR arm (12 samples) ----
sample_sheet <- read.csv(METH_SHEET, stringsAsFactors = FALSE)
meta <- sample_sheet %>%
  filter(Drug == "5-Azacytidine", CAR == "B7H3") %>%
  transmute(Code, Donor = factor(Donor), dose_nM = as.numeric(sub("nM", "", Dose)))
stopifnot(nrow(meta) == 12)
setorder(setDT(meta), Donor, dose_nM)

betas <- readRDS(BETAS_RDS)
betas_sub <- betas[intersect(unique(gmap$Name), rownames(betas)), meta$Code]
rm(betas); gc()

## ---- 3. Genome-wide per-gene median beta (reference for the column-Z heatmaps) ----
# 12_candidate_gene_heatmaps.R Z-scores each candidate-gene sample column
# against this matrix's column mean/SD.
all_gene_beta <- gene_median_beta(betas_sub, gmap)
fwrite(as.data.table(all_gene_beta, keep.rownames = "symbol"), file.path(out_dir, "meth_beta_wide_all_genes.csv"))
cat(sprintf("Genome-wide methylation reference: %d genes with >=1 promoter/body probe\n", nrow(all_gene_beta)))

## ---- 4. Candidate genes ----
beta_dt <- as.data.table(all_gene_beta[found, , drop = FALSE], keep.rownames = "symbol")
fwrite(beta_dt, file.path(out_dir, "meth_beta_wide.csv"))
fwrite(meta[, .(sample = Code, donor = Donor, dose_nM)], file.path(out_dir, "meth_sample_meta.csv"))

n_probes <- gmap[symbol %in% found & Name %in% rownames(betas_sub), .(n_probes = .N), by = symbol]
fwrite(n_probes, file.path(out_dir, "meth_n_probes.csv"))

## ---- 5. Per-dose trajectory (mean beta across the 3 donors at that dose) ----
long <- melt(beta_dt, id.vars = "symbol", variable.name = "sample", value.name = "beta")
long <- merge(long, meta[, .(sample = Code, donor = Donor, dose_nM)], by = "sample")
traj <- long[, .(mean_beta = mean(beta, na.rm = TRUE), n_donors = sum(!is.na(beta))), by = .(symbol, dose_nM)]
setorder(traj, symbol, dose_nM)
fwrite(traj, file.path(out_dir, "meth_dose_trajectory.csv"))

cat(sprintf("Saved methylation beta wide matrix (%d genes x %d samples), sample meta, and dose trajectory to %s\n",
            nrow(beta_dt), ncol(beta_dt) - 1, out_dir))
cat("=== Done ===\n")
