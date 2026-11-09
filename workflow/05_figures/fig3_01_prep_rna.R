#!/usr/bin/env Rscript
# Figure 3 cytokine panel: RNA side. Env `deseq2`.
# Per-sample normalized counts for Il2 / Ifng / Tnf + their dose-model stats, for both conditions.
#
# Normalization:
#   DESeq2 median-of-ratios size factors estimated ONCE on all 30 samples together (18 unstimulated + 12
#   stimulated; same >=10 counts in >=3 samples gene filter), normalized count = raw / size factor.
# Stats are unchanged -- they come from the per-condition dose models (scaling doesn't affect them):
#   Same design in both: ~ donor + dose_100nM, apeglm log2FC per 100 nM.
suppressPackageStartupMessages({ library(data.table); library(DESeq2) })
source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
deseq_dir <- RNA_DE_DIR
out_dir <- file.path(FIG_ROOT, "fig3_cytokines", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
GENES <- c("Il2", "Ifng", "Tnf")
CONDS <- c(unstimulated = "Unstimulated", stimulated = "GL261-B7H3 Stimulated")

## ---- joint size factors on all 30 samples ----
manifest <- fread(RNASEQ_MANIFEST)
meta <- unique(manifest[treatment_2 %in% CONDS, .(sample = gsub("-", "_", external_event_id),
                                                  donor = external_participant_id,
                                                  dose_nM = as.numeric(sub(" nM", "", dose_1)),
                                                  condition = names(CONDS)[match(treatment_2, CONDS)])])
stopifnot(nrow(meta) == 30, sum(meta$condition == "stimulated") == 12)
dt <- fread(RNASEQ_COUNTS, check.names = FALSE)
mat <- as.matrix(dt[, meta$sample, with = FALSE]); rownames(mat) <- dt$feature
storage.mode(mat) <- "integer"
mat <- mat[rowSums(mat >= 10) >= 3, , drop = FALSE]
dds_all <- estimateSizeFactors(DESeqDataSetFromMatrix(mat, as.data.frame(meta), design = ~ 1))
sf <- sizeFactors(dds_all)
fwrite(data.table(meta, size_factor_joint = sf[meta$sample]), file.path(out_dir, "rna_size_factors_joint.csv"))
cat(sprintf("Joint size factors (%d genes, 30 samples): %.2f-%.2f\n", nrow(mat), min(sf), max(sf)))

## ---- per condition: joint-normalized counts + that condition's dose stats ----
for (label in names(CONDS)) {
  fit <- c(unstimulated = "gene_NS", stimulated = "gene_S1")[[label]]
  res <- fread(file.path(deseq_dir, fit, "results_full.csv"))[symbol %in% GENES]
  stopifnot(setequal(res$symbol, GENES))
  m <- meta[condition == label]
  long <- rbindlist(lapply(seq_len(nrow(res)), function(i) {
    raw <- mat[res$feature[i], m$sample]
    data.table(symbol = res$symbol[i], m[, .(sample, donor, dose_nM)], raw_count = raw,
               norm_count = raw / sf[m$sample])
  }))
  fwrite(long, file.path(out_dir, sprintf("rna_%s_long.csv", label)))
  fwrite(res[, .(symbol, baseMean, log2FC_per100nM = log2FoldChange, lfcSE, pvalue, padj)],
         file.path(out_dir, sprintf("rna_%s_stats.csv", label)))
  cat("==", label, "(mean joint-normalized counts)\n")
  print(long[, .(mean_norm = round(mean(norm_count), 1)), by = symbol])
}

## ---- type I IFNs -> rna_<cond>_long_ifn.csv ----
# Ifnb1 and the Ifna family pooled (sum of raw counts over Ifna1-16 + Ifnab; no single Ifna gene is expressed).
# None passes the >=10-in->=3 filter, so they are not in the DESeq2 fits (no dose stats); counts are taken from
# the unfiltered matrix and divided by the same joint size factors.
suppressPackageStartupMessages(library(org.Mm.eg.db))
ens <- sub("\\..*", "", dt$feature)
sym <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, unique(ens), "SYMBOL", "ENSEMBL"))[ens]
ifn_groups <- list(Ifnb1 = which(sym == "Ifnb1"), Ifna = which(grepl("^Ifna([0-9]+|b)$", sym)))
cat("Ifna family genes pooled:", paste(sort(sym[ifn_groups$Ifna]), collapse = ", "), "\n")
full <- as.matrix(dt[, meta$sample, with = FALSE])
for (label in names(CONDS)) {
  m <- meta[condition == label]
  long <- rbindlist(lapply(names(ifn_groups), function(g) {
    raw <- colSums(full[ifn_groups[[g]], m$sample, drop = FALSE])
    data.table(symbol = g, m[, .(sample, donor, dose_nM)], raw_count = raw, norm_count = raw / sf[m$sample])
  }))
  fwrite(long, file.path(out_dir, sprintf("rna_%s_long_ifn.csv", label)))
  cat("==", label, "type I IFN (total raw reads, samples > 0)\n")
  print(long[, .(total_reads = sum(raw_count), n_nonzero = sum(raw_count > 0), n = .N), by = symbol])
}
