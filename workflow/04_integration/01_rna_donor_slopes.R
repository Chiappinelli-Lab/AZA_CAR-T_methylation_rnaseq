#!/usr/bin/env Rscript
# Per-donor and average-across-donor AZA dose-response slopes, RNA-seq gene
# expression, NS (non-stimulated) samples. Feeds the RNA side of the
# RNA-seq <-> methylation per-gene correlation (see integration/scripts/03).
#
# Two donor-scope summaries are computed from the SAME slope-fitting method
# (vectorized closed-form OLS of log2(normalized count+1) ~ dose_nM, fit
# independently within each donor):
#   - rna_slope_EV07: the EV07-only slope (the one donor shared with the
#     methylation cohort).
#   - rna_slope_avg_donors: the mean of the 5 independent per-donor slopes
#     (EV07/17/18/19/20). EV19/EV20 are missing the 200nM NS sample, so their
#     slopes are fit from only 3 dose points (0/50/100nM) instead of 4.
# Neither uses DESeq2's own donor-blocked pooled dose_100nM coefficient
# (results_full.csv$log2FoldChange) -- that coefficient is still carried
# through here (for reference) and used in 03 to define the DE-gene set via
# the existing `sig` (padj<0.05) call, but is not itself the plotted value.

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- INT_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load DESeq2 object, get normalized counts ----
dds <- readRDS(file.path(RNA_NS_DIR, "dds.rds"))
nc <- counts(dds, normalized = TRUE)
log2nc <- log2(nc + 1)

cd <- as.data.frame(colData(dds))
stopifnot(identical(rownames(cd), colnames(log2nc)))
cat("Genes x samples:", nrow(log2nc), "x", ncol(log2nc), "\n")
cat("Donors:", paste(levels(cd$donor), collapse = ", "), "\n")

## ---- 2. Vectorized per-gene OLS slope, fit independently within each donor ----
# Closed-form simple linear regression slope: slope = sum((x-xbar)*(y-ybar)) / sum((x-xbar)^2)
# Applied to all genes at once via a single matrix-vector product per donor.
fit_donor_slope <- function(y_mat, x) {
  xc <- x - mean(x)
  as.vector(y_mat %*% xc) / sum(xc^2)
}

donors <- levels(cd$donor)
slope_mat <- matrix(NA_real_, nrow = nrow(log2nc), ncol = length(donors),
                     dimnames = list(rownames(log2nc), donors))

for (d in donors) {
  samp <- rownames(cd)[cd$donor == d]
  n_doses <- length(samp)
  cat(sprintf("Donor %s: %d dose points (%s)\n", d, n_doses,
              paste(sort(cd[samp, "dose_nM"]), collapse = ",")))
  slope_mat[, d] <- fit_donor_slope(log2nc[, samp, drop = FALSE], cd[samp, "dose_nM"])
}

rna_slope_EV07 <- slope_mat[, "EV07"]
rna_slope_avg_donors <- rowMeans(slope_mat)

## ---- 3. Attach gene symbols / existing DE significance call ----
res_full <- fread(file.path(RNA_NS_DIR, "results_full.csv"))
stopifnot(all(res_full$feature %in% rownames(slope_mat)))

dt <- data.table(
  feature = rownames(slope_mat),
  rna_slope_EV07 = rna_slope_EV07,
  rna_slope_avg_donors = rna_slope_avg_donors
)
dt <- merge(dt, res_full[, .(feature, ensembl_id, symbol, baseMean, log2FoldChange, padj, DE_sig = sig)],
            by = "feature", all.x = TRUE)

n_before <- nrow(dt)
dt <- dt[!is.na(symbol) & symbol != ""]
cat(sprintf("Dropped %d / %d rows with NA/empty symbol\n", n_before - nrow(dt), n_before))

dup_symbols <- dt[, .N, by = symbol][N > 1]
cat(sprintf("N symbols with >1 Ensembl ID: %d (%d rows) -- keeping max-baseMean row per symbol\n",
            nrow(dup_symbols), sum(dup_symbols$N)))
setorder(dt, symbol, -baseMean)
dt <- dt[, .SD[1], by = symbol]

cat(sprintf("Final RNA gene-slope table: %d unique gene symbols\n", nrow(dt)))

## ---- 4. Save ----
out_cols <- c("symbol", "ensembl_id", "baseMean", "log2FoldChange", "padj", "DE_sig",
              "rna_slope_EV07", "rna_slope_avg_donors")
fwrite(dt[, ..out_cols], file.path(out_dir, "rna_gene_slopes.csv"))
cat("Saved", file.path(out_dir, "rna_gene_slopes.csv"), "\n")
cat("=== Done ===\n")
