#!/usr/bin/env Rscript
# Per-sample + per-dose RNA-seq data for the  candidate gene sets
# feeding the RNA panel
# of the per-pathway 3-panel (expression / methylation / correlation)
# heatmaps built in 12_candidate_gene_heatmaps.R.
#
# Same log2(normalized count + 1), NS (non-stimulated), all 5-donor dds
# object as 01_rna_donor_slopes.R -- but here we keep the full per-sample
# matrix (not just the per-donor slope) plus a per-dose (averaged across
# donors) trajectory, since the heatmap script needs both.

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- file.path(INT_DIR, "candidate_genes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(WORKFLOW, "04_integration", "candidate_gene_sets.R"))
all_genes <- sort(unique(unlist(candidate_gene_sets)))
cat(sprintf("N candidate genes (union across %d gene sets): %d\n", length(candidate_gene_sets), length(all_genes)))

## ---- 1. Load DESeq2 object, normalized counts ----
dds <- readRDS(file.path(RNA_NS_DIR, "dds.rds"))
nc <- counts(dds, normalized = TRUE)
log2nc <- log2(nc + 1)
cd <- as.data.frame(colData(dds))
stopifnot(identical(rownames(cd), colnames(log2nc)))

## ---- 2. Map candidate symbols -> Ensembl feature rows ----
res_full <- fread(file.path(RNA_NS_DIR, "results_full.csv"))
map <- res_full[symbol %in% all_genes, .(feature, symbol, baseMean)]
n_before <- length(all_genes)
found <- unique(map$symbol)
missing <- setdiff(all_genes, found)
if (length(missing) > 0) cat("Candidate genes NOT found in RNA-seq results_full.csv symbol column:\n  ", paste(missing, collapse = ", "), "\n")

# Keep max-baseMean feature per symbol (dup Ensembl IDs mapping to one symbol)
setorder(map, symbol, -baseMean)
map <- map[, .SD[1], by = symbol]
stopifnot(all(map$feature %in% rownames(log2nc)))

## ---- 3. Per-sample wide matrix (genes x 18 NS samples) ----
expr_mat <- log2nc[map$feature, , drop = FALSE]
rownames(expr_mat) <- map$symbol
expr_dt <- as.data.table(expr_mat, keep.rownames = "symbol")
fwrite(expr_dt, file.path(out_dir, "rna_expr_wide.csv"))

# Genome-wide reference for the column-Z heatmaps: every
# gene in the dds (the 18,267 passing >=10 counts in >=3 samples), keyed by
# Ensembl feature ID; 12_candidate_gene_heatmaps.R Z-scores each sample column
# of the candidate genes against this matrix's column mean/SD.
fwrite(as.data.table(log2nc, keep.rownames = "feature"), file.path(out_dir, "rna_expr_wide_all_genes.csv"))

meta_dt <- data.table(sample = rownames(cd), donor = cd$donor, dose_nM = cd$dose_nM)
setorder(meta_dt, donor, dose_nM)
fwrite(meta_dt, file.path(out_dir, "rna_sample_meta.csv"))

## ---- 4. Per-dose trajectory (mean log2 expr across donors present at that dose) ----
long <- melt(expr_dt, id.vars = "symbol", variable.name = "sample", value.name = "log2expr")
long <- merge(long, meta_dt, by = "sample")
traj <- long[, .(mean_log2expr = mean(log2expr), n_donors = .N), by = .(symbol, dose_nM)]
setorder(traj, symbol, dose_nM)
fwrite(traj, file.path(out_dir, "rna_dose_trajectory.csv"))

cat(sprintf("Saved RNA expression wide matrix (%d genes x %d samples), sample meta, and dose trajectory to %s\n",
            nrow(expr_dt), ncol(expr_dt) - 1, out_dir))
cat("=== Done ===\n")
