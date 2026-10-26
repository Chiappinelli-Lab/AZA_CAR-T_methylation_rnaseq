#!/usr/bin/env Rscript
# DESeq2 gene-level AZA dose-response test in STIMULATED (S1: one 2:1 GL261-B7H3
# co-culture) CAR-T cells. Added 2026-09-24 for the collaborators' figure 3
# (IL-2 / IFN-g / TNF expression). Mirrors run_deseq2_gene_NS.R exactly:
# counts ~ donor + dose_100nM (dose continuous, nM/100), >=10 counts in >=3
# samples, apeglm shrinkage, sig = padj<0.05 & |log2FC| > 0.5.
# 12 samples: 3 donors (EV07, EV17, EV20) x 0/50/100/200 nM -- fewer donors than
# NS (5), so less power; EV18/EV19 have no stimulated samples.
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(org.Mm.eg.db)
})
set.seed(1)
source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
manifest_path <- RNASEQ_MANIFEST
counts_path <- RNASEQ_COUNTS
cat("Counts:", counts_path, "\n")
out_dir <- RNA_S1_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. colData for S1 samples ----
manifest <- fread(manifest_path, sep = "\t", header = TRUE)
s1 <- unique(manifest[treatment_2 == "GL261-B7H3 Stimulated",
                      .(external_event_id, external_participant_id, dose_1)])
stopifnot(nrow(s1) == 12, !anyDuplicated(s1$external_event_id))
coldata <- data.frame(
  sample = gsub("-", "_", s1$external_event_id),
  donor = factor(s1$external_participant_id),
  dose_nM = as.numeric(sub(" nM", "", s1$dose_1)),
  row.names = gsub("-", "_", s1$external_event_id)
)
coldata$dose_100nM <- coldata$dose_nM / 100
coldata <- coldata[order(coldata$donor, coldata$dose_nM), ]

## ---- 2. Counts, filter, fit ----
dt <- fread(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
mat <- as.matrix(dt[, coldata$sample, with = FALSE])
rownames(mat) <- dt$feature
storage.mode(mat) <- "integer"
keep <- rowSums(mat >= 10) >= 3
message(sprintf("%d / %d genes pass >=10 counts in >=3 samples filter", sum(keep), nrow(mat)))
mat <- mat[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~ donor + dose_100nM)
dds <- DESeq(dds)
res <- results(dds, name = "dose_100nM", alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "dose_100nM", res = res, type = "apeglm")

## ---- 3. Annotate + save ----
res_df <- as.data.frame(res_shrunk)
res_df$feature <- rownames(res_df)
res_df$ensembl_id <- sub("\\..*$", "", res_df$feature)
res_df$symbol <- mapIds(org.Mm.eg.db, keys = res_df$ensembl_id, column = "SYMBOL", keytype = "ENSEMBL",
                        multiVals = "first")
res_df$gene_name <- mapIds(org.Mm.eg.db, keys = res_df$ensembl_id, column = "GENENAME", keytype = "ENSEMBL",
                           multiVals = "first")
LFC_CUTOFF <- 0.5
res_df$sig <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > LFC_CUTOFF
res_df <- res_df[order(res_df$padj), c("feature", "ensembl_id", "symbol", "gene_name", "baseMean",
                                       "log2FoldChange", "lfcSE", "pvalue", "padj", "sig")]
fwrite(res_df, file.path(out_dir, "results_full.csv"))
fwrite(res_df[res_df$sig, ], file.path(out_dir, "results_sig_padj0.05.csv"))
saveRDS(dds, file.path(out_dir, "dds.rds"))
cat(sprintf("S1: %d genes tested, %d with padj<0.05, %d sig (padj<0.05 & |LFC|>%.1f)\n",
            sum(!is.na(res_df$padj)), sum(res_df$padj < 0.05, na.rm = TRUE), sum(res_df$sig), LFC_CUTOFF))
