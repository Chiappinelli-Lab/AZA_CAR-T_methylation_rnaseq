#!/usr/bin/env Rscript
# 3-panel (RNA expression | DNA methylation | RNA<->Meth correlation) heatmap
# per candidate gene set from the custom pathway list ("Candidate Genes Methyl Array
# RNA Seq EV.docx"). Inputs built by 10_candidate_gene_rna_data.R (deseq2
# env) and 11_candidate_gene_meth_data.R (haydar_methylation env).
#
# Panel 1 (RNA):  row-Z-scored, donor-corrected log2 normalized counts, NS
#   samples, all 5 donors x 4 doses (18 samples; EV19/EV20 missing the
#   200nM point), columns grouped by dose (donor nested within each dose
#   block). Donor correction via limma::removeBatchEffect(batch = donor,
#   design = ~dose_nM): removes the fitted donor baseline
#   while preserving the dose effect + real per-sample residual noise.
# Panel 2 (Meth): row-Z-scored, donor-corrected per-gene composite beta
#   (median beta over the gene's promoter 1.5kb + gene-body probes, dmGsea
#   mapping; was mean region beta over DMRs within 5kb before 2026-09-24),
#   AZA + B7H3 CAR arm, 3 donors x 4 doses (12 samples),
#   columns grouped by dose the same way. No column annotation (redundant
#   with panel 1, which shares the same dose grouping/column order logic).
# Panel 3 (Corr): Pearson r, per gene, between the RNA and methylation
#   AZA dose-response trajectories (mean-across-donors value at each of
#   0/50/100/200nM, n=4 points, computed from the DONOR-CORRECTED matrices
#   of panels 1/2 -- gene-level dose-response concordance
#   across the two independent cohorts,
#   NOT a per-sample/per-donor paired correlation.
#   p = two-sided cor.test() p-value (added 2026-10-09); with n=4 (df=2) it
#   needs |r| > ~0.95 for p<0.05, so it is saved for reference, not plotted.
#
# Genes with no RNA data (below DESeq2's >=10 counts in >=3 samples filter) or no
# methylation data (no promoter/body probe) are shown as blank/grey rows in the
# corresponding panel; row order is fixed to the list order (not
# clustered), so gene identity stays legible across all 3 panels. Row
# labels are bold if either flag applies, and annotated with * (DE gene,
# RNA padj<0.05 & |log2FC|>0.5, from
# rna_gene_slopes.csv$DE_sig) and/or † (>=1 promoter/body probe inside a
# significant DMR -- FDR-filtered + |meandiff|>0.10 -- from
# methylation_gene_slopes.csv$is_sig_DMR_gene).

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(data.table)
  library(dplyr)
  library(limma)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
results_dir <- INT_DIR
in_dir <- file.path(results_dir, "candidate_genes")
plot_dir <- file.path(results_dir, "plots/candidate_genesets")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(WORKFLOW, "04_integration", "candidate_gene_sets.R"))

## ---- 1. Load data ----
rna_wide  <- fread(file.path(in_dir, "rna_expr_wide.csv"))
rna_meta  <- fread(file.path(in_dir, "rna_sample_meta.csv"))

meth_wide <- fread(file.path(in_dir, "meth_beta_wide.csv"))
meth_meta <- fread(file.path(in_dir, "meth_sample_meta.csv"))

# Significance flags for gene-label annotation (computed on the full gene
# universe, not just the candidate list -- see 01/02_..._donor_slopes.R)
rna_sig  <- fread(file.path(results_dir, "rna_gene_slopes.csv"))[, .(symbol, DE_sig)]
meth_sig <- fread(file.path(results_dir, "methylation_gene_slopes.csv"))[, .(symbol, is_sig_DMR_gene)]

# Group columns by dose (donor nested within each dose block)
setorder(rna_meta, dose_nM, donor)
setorder(meth_meta, dose_nM, donor)

rna_mat_full  <- as.matrix(rna_wide[, -"symbol"]);  rownames(rna_mat_full)  <- rna_wide$symbol
meth_mat_full <- as.matrix(meth_wide[, -"symbol"]); rownames(meth_mat_full) <- meth_wide$symbol
rna_mat_full  <- rna_mat_full[,  rna_meta$sample,  drop = FALSE]
meth_mat_full <- meth_mat_full[, meth_meta$sample, drop = FALSE]

## ---- 1b. Donor correction (removeBatchEffect, dose preserved) ----
rna_mat_full <- removeBatchEffect(rna_mat_full, batch = factor(rna_meta$donor),
                                   design = model.matrix(~dose_nM, data = rna_meta))
meth_mat_full <- removeBatchEffect(meth_mat_full, batch = factor(meth_meta$donor),
                                    design = model.matrix(~dose_nM, data = meth_meta))

## ---- 2. Per-gene RNA<->Meth dose-trajectory correlation (all candidate genes) ----
# Trajectories = per-dose mean of the donor-corrected (not Z-scored) values.
dose_means <- function(mat, meta, value_name) {
  long <- as.data.table(mat, keep.rownames = "symbol")
  long <- melt(long, id.vars = "symbol", variable.name = "sample", value.name = "v", variable.factor = FALSE)
  long <- merge(long, meta[, .(sample, dose_nM)], by = "sample")
  setnames(long[, .(v = mean(v)), by = .(symbol, dose_nM)], "v", value_name)[]
}
rna_traj  <- dose_means(rna_mat_full,  rna_meta,  "mean_log2expr")
meth_traj <- dose_means(meth_mat_full, meth_meta, "mean_beta")
all_genes <- sort(unique(unlist(candidate_gene_sets)))
corr_dt <- rbindlist(lapply(all_genes, function(g) {
  rt <- rna_traj[symbol == g, .(dose_nM, mean_log2expr)]
  mt <- meth_traj[symbol == g, .(dose_nM, mean_beta)]
  if (nrow(rt) == 0 || nrow(mt) == 0) return(data.table(symbol = g, r = NA_real_, p = NA_real_, n = 0L))
  j <- merge(rt, mt, by = "dose_nM")
  if (nrow(j) < 3 || sd(j$mean_log2expr) == 0 || sd(j$mean_beta) == 0) {
    return(data.table(symbol = g, r = NA_real_, p = NA_real_, n = nrow(j)))
  }
  ct <- cor.test(j$mean_log2expr, j$mean_beta, method = "pearson")
  data.table(symbol = g, r = unname(ct$estimate), p = ct$p.value, n = nrow(j))
}))
fwrite(corr_dt, file.path(in_dir, "rna_meth_dose_trajectory_correlation.csv"))

## ---- 3. Colors (shared across all gene-set figures) ----
rna_donor_colors  <- c(EV07 = "#E87BA4", EV17 = "#2A78D6", EV18 = "#008300", EV19 = "#EDA100", EV20 = "#7B4FA6")
meth_donor_colors <- c(EV07 = "#E87BA4", EV08 = "#EDA100", ZY03 = "#008300")
dose_levels <- c(0, 50, 100, 200)
dose_ramp <- colorRamp2(c(1, 4), c("#EDEFEC", "#2A3A33"))
dose_colors <- setNames(dose_ramp(1:4), dose_levels)

z_ramp <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
r_ramp <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

rna_col_ann <- HeatmapAnnotation(
  Donor = rna_meta$donor, Dose = factor(rna_meta$dose_nM, levels = dose_levels),
  col = list(Donor = rna_donor_colors, Dose = dose_colors),
  annotation_name_side = "left", show_legend = c(Donor = TRUE, Dose = TRUE)
)
# Panel 2 (methylation) keeps its own Donor/Dose color bars -- its donors
# (EV07/EV08/ZY03) differ from panel 1's (EV07/17/18/19/20), so it keeps its
# own legend, but the "Donor"/"Dose" row-side text label is redundant next
# to panel 1's and is suppressed here.
meth_col_ann <- HeatmapAnnotation(
  Donor = meth_meta$donor, Dose = factor(meth_meta$dose_nM, levels = dose_levels),
  col = list(Donor = meth_donor_colors, Dose = dose_colors),
  show_annotation_name = FALSE
)

row_zscore <- function(m) t(scale(t(m)))
# Column Z: each sample column Z-scored against the
# GENOME-WIDE distribution of that sample -- mean/SD over all 18,267 expressed
# genes (RNA) or all genes with a promoter/body probe (methylation), after the same
# donor correction -- so colors say how high a gene sits relative to the
# genome, comparable across gene sets. (First version Z-scored within each
# ~10-gene set; switched at the user's request.) Written as
# *_3panel_heatmap_colz.{png,pdf}.
genome_col_stats <- function(file, key, meta) {
  w <- fread(file.path(in_dir, file))
  m <- as.matrix(w[, -key, with = FALSE])[, meta$sample, drop = FALSE]
  m <- removeBatchEffect(m, batch = factor(meta$donor), design = model.matrix(~dose_nM, data = meta))
  list(mu = colMeans(m, na.rm = TRUE), sd = apply(m, 2, sd, na.rm = TRUE), n = sum(rowSums(!is.na(m)) > 0))
}
rna_ref  <- genome_col_stats("rna_expr_wide_all_genes.csv",  "feature", rna_meta)
meth_ref <- genome_col_stats("meth_beta_wide_all_genes.csv", "symbol",  meth_meta)
cat(sprintf("Column-Z reference: RNA %d genes, methylation %d genes\n", rna_ref$n, meth_ref$n))
col_zscore <- function(m, ref) sweep(sweep(m, 2, ref$mu[colnames(m)], "-"), 2, ref$sd[colnames(m)], "/")
norm_modes <- list(
  row = list(fn = function(m, ref) row_zscore(m), suffix = "", legend = "z-score",
             subtitle = "Z-scored per gene (row): each gene vs. its own dose trajectory"),
  col = list(fn = col_zscore, suffix = "_colz", legend = "genome z",
             subtitle = "Z-scored per sample (column) against all measured genes: level relative to the genome"))

## ---- 4. Build + draw one 3-panel figure per gene set, per normalization ----
for (mode in names(norm_modes)) {
nm <- norm_modes[[mode]]
rna_leg <- paste0("RNA\n", nm$legend); meth_leg <- paste0("Meth.\n", nm$legend)
for (set_name in names(candidate_gene_sets)) {
  genes <- candidate_gene_sets[[set_name]]

  rna_sub <- matrix(NA_real_, nrow = length(genes), ncol = ncol(rna_mat_full),
                     dimnames = list(genes, colnames(rna_mat_full)))
  have_rna <- intersect(genes, rownames(rna_mat_full))
  rna_sub[have_rna, ] <- rna_mat_full[have_rna, , drop = FALSE]
  rna_z <- nm$fn(rna_sub, rna_ref)
  rna_z[is.nan(rna_z)] <- NA

  meth_sub <- matrix(NA_real_, nrow = length(genes), ncol = ncol(meth_mat_full),
                      dimnames = list(genes, colnames(meth_mat_full)))
  have_meth <- intersect(genes, rownames(meth_mat_full))
  meth_sub[have_meth, ] <- meth_mat_full[have_meth, , drop = FALSE]
  meth_z <- nm$fn(meth_sub, meth_ref)
  meth_z[is.nan(meth_z)] <- NA

  r_vec <- corr_dt$r[match(genes, corr_dt$symbol)]
  r_mat <- matrix(r_vec, ncol = 1, dimnames = list(genes, "r"))
  r_label <- ifelse(is.na(r_vec), "", sprintf("%.2f", r_vec))

  de_flag  <- rna_sig$DE_sig[match(genes, rna_sig$symbol)]
  de_flag[is.na(de_flag)] <- FALSE
  dmr_flag <- meth_sig$is_sig_DMR_gene[match(genes, meth_sig$symbol)]
  dmr_flag[is.na(dmr_flag)] <- FALSE
  gene_labels <- paste0(genes, ifelse(de_flag, "*", ""), ifelse(dmr_flag, "†", ""))

  pathway_label <- gsub("^[0-9]+_", "", set_name)
  pathway_label <- gsub("_", " ", pathway_label)

  ht_rna <- Heatmap(rna_z, name = rna_leg, col = z_ramp, na_col = "grey85",
                     top_annotation = rna_col_ann,
                     column_title = "Gene expression\n(NS, donor-corrected log2 norm. counts)",
                     column_title_gp = grid::gpar(fontsize = 10.5),
                     row_title = pathway_label, row_title_gp = grid::gpar(fontface = "bold", fontsize = 11),
                     row_labels = gene_labels,
                     cluster_rows = FALSE, cluster_columns = FALSE,
                     show_row_names = TRUE, row_names_side = "left", show_column_names = FALSE,
                     # Bold = DE gene or sig-DMR gene (2026-09-22); */† still say which.
                     row_names_gp = grid::gpar(fontsize = 10,
                                               fontface = ifelse(de_flag | dmr_flag, "bold", "plain")))

  ht_meth <- Heatmap(meth_z, name = meth_leg, col = z_ramp, na_col = "grey85",
                      top_annotation = meth_col_ann,
                      column_title = "DNA methylation\n(AZA+B7H3, donor-corrected median\npromoter + gene-body beta)",
                      column_title_gp = grid::gpar(fontsize = 10.5),
                      cluster_rows = FALSE, cluster_columns = FALSE,
                      show_row_names = FALSE, show_column_names = FALSE)

  ht_corr <- Heatmap(r_mat, name = "RNA-Meth\ndose corr.\n(Pearson r)", col = r_ramp, na_col = "grey85",
                      column_title = "Corr.",
                      column_title_gp = grid::gpar(fontsize = 10.5),
                      cluster_rows = FALSE, cluster_columns = FALSE,
                      show_row_names = FALSE, show_column_names = FALSE,
                      width = unit(1.6, "cm"),
                      cell_fun = function(j, i, x, y, w, h, fill) {
                        grid::grid.text(r_label[i], x, y, gp = grid::gpar(fontsize = 9))
                      })

  ht_list <- ht_rna + ht_meth + ht_corr

  footnote <- "Bold = DE and/or DMR gene.   * DE gene (RNA-seq padj < 0.05 & |log2FC| > 0.5)      † promoter/body probe in a significant DMR (FDR-filtered, |meandiff| > 0.10)"

  out_base <- file.path(plot_dir, paste0(set_name, "_3panel_heatmap", nm$suffix))
  png(paste0(out_base, ".png"), width = 14.5, height = 6.4, units = "in", res = 300)
  draw(ht_list, main_heatmap = rna_leg, merge_legend = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       ht_gap = unit(c(9, 9), "mm"), padding = unit(c(16, 2, 2, 2), "mm"),
       column_title = paste0(pathway_label, ": AZA dose-response, RNA-seq vs. methylation\n", nm$subtitle),
       column_title_gp = grid::gpar(fontface = "bold", fontsize = 13))
  grid::grid.text(footnote, x = unit(5, "mm"), y = unit(6, "mm"),
                   just = c("left", "bottom"), gp = grid::gpar(fontsize = 9, fontface = "italic"))
  dev.off()
  cairo_pdf(paste0(out_base, ".pdf"), width = 14.5, height = 6.4)
  draw(ht_list, main_heatmap = rna_leg, merge_legend = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       ht_gap = unit(c(9, 9), "mm"), padding = unit(c(16, 2, 2, 2), "mm"),
       column_title = paste0(pathway_label, ": AZA dose-response, RNA-seq vs. methylation\n", nm$subtitle),
       column_title_gp = grid::gpar(fontface = "bold", fontsize = 13))
  grid::grid.text(footnote, x = unit(5, "mm"), y = unit(6, "mm"),
                   just = c("left", "bottom"), gp = grid::gpar(fontsize = 9, fontface = "italic"))
  dev.off()
  cat("Saved", out_base, ".png/.pdf\n")
}
}
cat("=== Done ===\n")
