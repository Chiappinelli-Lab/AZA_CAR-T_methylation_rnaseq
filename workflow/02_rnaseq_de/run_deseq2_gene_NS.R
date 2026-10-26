#!/usr/bin/env Rscript
# DESeq2 gene-level AZA dose-response test in non-stimulated (NS) CAR-T cells.
# Design: counts ~ donor + dose_100nM (dose as a continuous covariate, nM/100).
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(pheatmap)
  library(org.Mm.eg.db)
})

set.seed(1)

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
manifest_path <- RNASEQ_MANIFEST
counts_path <- RNASEQ_COUNTS
cat("Counts:", counts_path, "\n")
out_dir <- RNA_NS_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Build colData for NS samples ----
manifest <- fread(manifest_path, sep = "\t", header = TRUE)
ns <- unique(manifest[treatment_2 == "Unstimulated",
                       .(external_event_id, external_participant_id, dose_1)])
stopifnot(nrow(ns) == 18, !anyDuplicated(ns$external_event_id))

coldata <- data.frame(
  sample = gsub("-", "_", ns$external_event_id),
  donor = factor(ns$external_participant_id),
  dose_nM = as.numeric(sub(" nM", "", ns$dose_1)),
  row.names = gsub("-", "_", ns$external_event_id)
)
coldata$dose_100nM <- coldata$dose_nM / 100
coldata <- coldata[order(coldata$donor, coldata$dose_nM), ]

## ---- 2. Load counts, subset to NS samples ----
dt <- fread(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
mat <- as.matrix(dt[, coldata$sample, with = FALSE])
rownames(mat) <- dt$feature
storage.mode(mat) <- "integer"

## ---- 3. Filter low-count genes ----
# Filter: >=10 raw counts in at least 3 samples (confirmed with the user
# 2026-09-22, replacing rowSums>=10). rowSums>=10 let through genes with a
# handful of reads in a single sample and zeros elsewhere; those produced
# spurious |log2FC| ~ 18-28 calls with small p-values (see
# results/shrinkage_ashr_vs_apeglm/ and plot_ashr_tail_gene_counts.R).
keep <- rowSums(mat >= 10) >= 3
message(sprintf("%d / %d genes pass >=10 counts in >=3 samples filter", sum(keep), nrow(mat)))
mat <- mat[keep, , drop = FALSE]

## ---- 4. Fit DESeq2 model ----
dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata,
                               design = ~ donor + dose_100nM)
dds <- DESeq(dds)
message("resultsNames(dds): ", paste(resultsNames(dds), collapse = ", "))

## ---- 5. Extract AZA dose-response effect ----
res <- results(dds, name = "dose_100nM", alpha = 0.05)
# apeglm shrinkage (confirmed with the user 2026-09-22, replacing ashr): ashr
# only sees the MLE + SE, so all-zero-in-one-group genes kept |LFC| ~ 18-28;
# apeglm uses the count likelihood and shrinks them to ~0.
res_shrunk <- lfcShrink(dds, coef = "dose_100nM", res = res, type = "apeglm")

## ---- 6. Annotate with gene symbols ----
res_df <- as.data.frame(res_shrunk)
res_df$feature <- rownames(res_df)
res_df$ensembl_id <- sub("\\..*$", "", res_df$feature)
res_df$symbol <- mapIds(org.Mm.eg.db, keys = res_df$ensembl_id,
                         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res_df$gene_name <- mapIds(org.Mm.eg.db, keys = res_df$ensembl_id,
                            column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
# Significance call: padj<0.05 AND |log2FC per 100nM AZA| > 0.5 (confirmed
# with the user 2026-09-14, replacing the earlier padj-only call -- an
# effect-size floor on top of the FDR test. Kept in per-100nM units rather
# than re-expressing as "fold-change over 0->200nM" because not every
# donor has a 200nM sample (EV19/EV20 max out at 100nM); 0.5 per 100nM is
# equivalent to a 2-fold change by 200nM only for the donors that reach it.)
LFC_CUTOFF <- 0.5
res_df$sig <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > LFC_CUTOFF
res_df <- res_df[order(res_df$padj), ]
res_df <- res_df[, c("feature", "ensembl_id", "symbol", "gene_name",
                      "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "sig")]

fwrite(res_df, file.path(out_dir, "results_full.csv"))
fwrite(res_df[res_df$sig, ], file.path(out_dir, "results_sig_padj0.05.csv"))
saveRDS(dds, file.path(out_dir, "dds.rds"))

## ---- 7. QC / diagnostic plots ----
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("donor", "dose_100nM"), returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = dose_100nM * 100, shape = donor)) +
  geom_point(size = 3) +
  scale_color_viridis_c(name = "Dose (nM)") +
  xlab(sprintf("PC1: %d%% variance", pct_var[1])) +
  ylab(sprintf("PC2: %d%% variance", pct_var[2])) +
  ggtitle("NS gene-level PCA") +
  theme_bw()
ggsave(file.path(out_dir, "PCA.png"), p_pca, width = 6, height = 5)

sample_dist <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dist)
ann_row <- coldata[rownames(sample_dist_mat), c("donor", "dose_nM"), drop = FALSE]
png(file.path(out_dir, "sample_distance_heatmap.png"), width = 800, height = 700)
pheatmap(sample_dist_mat, annotation_row = ann_row, main = "NS sample-to-sample distance (vst)")
dev.off()

png(file.path(out_dir, "dispersion_plot.png"), width = 600, height = 500)
plotDispEsts(dds, main = "NS gene-level dispersion estimates")
dev.off()

p_ma <- ggplot(res_df, aes(baseMean, log2FoldChange, color = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_x_log10() +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
  ggtitle("NS AZA dose-response: MA plot") +
  theme_bw()
ggsave(file.path(out_dir, "MA_plot.png"), p_ma, width = 6, height = 5)

# Horizontal dashed line: the raw p-value equivalent of the padj<0.05
# boundary (max raw pvalue among padj<0.05 genes, since BH is monotonic in
# p) -- the volcano's y-axis is raw -log10(p), not -log10(padj), so this is
# the honest way to show where that FDR cutoff actually falls on this axis.
p_cutoff_at_padj05 <- if (any(res_df$padj < 0.05, na.rm = TRUE)) {
  max(res_df$pvalue[which(res_df$padj < 0.05)], na.rm = TRUE)
} else NA_real_

# log2FoldChange is heavily right/left-skewed by a handful of near-zero-
# baseline genes reaching |log2FC| up to ~20, while 95% of genes fall
# within +/-0.06 -- a linear x-axis squashes essentially every point into
# a sliver near 0. scales::pseudo_log_trans (a ggplot2 dependency, no new
# package needed) is linear near zero -- so the +/-1 cutoff lines and the
# bulk of points stay readable -- and compresses log-like further out, so
# the few extreme outliers no longer dominate the axis range.
volcano_breaks <- c(-16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16)
p_volcano <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 2), breaks = volcano_breaks) +
  ggtitle("NS AZA dose-response: volcano") +
  xlab("log2FC per 100 nM AZA (pseudo-log scale)") +
  theme_bw()
if (!is.na(p_cutoff_at_padj05)) {
  p_volcano <- p_volcano +
    geom_hline(yintercept = -log10(p_cutoff_at_padj05), linetype = "dashed", color = "grey40", linewidth = 0.3) +
    labs(subtitle = sprintf("dashed lines: padj<0.05 (p<%.2g) and |log2FC|>%g", p_cutoff_at_padj05, LFC_CUTOFF))
}
ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 6, height = 5)

p_pval <- ggplot(res_df, aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.025), fill = "steelblue", color = "white") +
  ggtitle("NS AZA dose-response: raw p-value distribution") +
  theme_bw()
ggsave(file.path(out_dir, "pvalue_histogram.png"), p_pval, width = 6, height = 4)

## ---- 8. Summary ----
n_tested <- sum(!is.na(res_df$padj))
n_sig_05 <- sum(res_df$padj < 0.05, na.rm = TRUE)
n_sig_10 <- sum(res_df$padj < 0.10, na.rm = TRUE)
n_up_05 <- sum(res_df$sig & res_df$log2FoldChange > 0)
n_down_05 <- sum(res_df$sig & res_df$log2FoldChange < 0)

summary_lines <- c(
  "DESeq2 gene-level AZA dose-response, NS (non-stimulated) samples",
  "Design: ~ donor + dose_100nM (continuous, log2FC per 100 nM AZA)",
  sprintf("Samples: %d (donors: %s)", nrow(coldata), paste(levels(coldata$donor), collapse = ", ")),
  sprintf("Genes passing >=10 counts in >=3 samples filter: %d", nrow(mat)),
  sprintf("Genes with non-NA padj (tested): %d", n_tested),
  sprintf("Significant at padj<0.05: %d (up=%d, down=%d)", n_sig_05, n_up_05, n_down_05),
  sprintf("Significant at padj<0.10: %d", n_sig_10)
)
writeLines(summary_lines, file.path(out_dir, "summary.txt"))
message(paste(summary_lines, collapse = "\n"))
message("=== Done ===")
