#!/usr/bin/env Rscript
# Paper figure (RNA-seq x methylation integration), methylation-side and
# gene-set plot-ready data. Run in env `haydar_methylation` after
# fig2_01_prep_rna.R. Added 2026-09-24.
#
# Writes to paper_figures/fig2_rna_methylation/data/:
#   dmr_heatmap_matrix.csv / dmr_heatmap_samples.csv -- panel B: the 507
#     filtered DMRs, region-average donor-adjusted beta, computed exactly as the
#     donor-adjusted heatmap in Haydar_Methylation_Array/scripts/08 (complete-
#     case probes, removeBatchEffect(batch = Donor, design = ~AZA), mean over
#     each DMR's member probes).
#   pathway_termsim_wang.csv / pathway_tree_hclust.csv -- panel E: Wang GO
#     semantic similarity of the 26 shared GO:BP terms, and the hclust tree
#     (ward.D on 1 - similarity, 5 clusters) -- the same clustering enrichplot
#     treeplot() used for meth_rna_shared_GO_BP_tree.png (12_dmgsea_plots.R).
#   geneset_{rna,meth}_corrected.csv -- panel G: every candidate gene, donor-
#     corrected (removeBatchEffect, design ~dose_nM), from
#     integration/results/candidate_genes/ (11_candidate_gene_meth_data.R:
#     median promoter + gene-body probe beta).
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(GenomicRanges)
  library(sesameData)
  library(GOSemSim)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- file.path(FIG_ROOT, "fig2_rna_methylation", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Panel B: 507 DMRs, donor-adjusted region beta ----
res <- readRDS(DMRCATE_RDS)
dmrs_filt <- res$dmrs[abs(mcols(res$dmrs)$meandiff) > 0.10]
ss <- fread(METH_SHEET)
meta <- ss[Drug == "5-Azacytidine" & CAR == "B7H3", .(Code, Donor, AZA = as.numeric(sub("nM", "", Dose)))]
stopifnot(nrow(meta) == 12)
betas <- readRDS(BETAS_RDS)[, meta$Code]
complete <- rownames(betas)[complete.cases(betas)]
betas_adj <- removeBatchEffect(betas[complete, ], batch = meta$Donor, design = model.matrix(~AZA, data = meta))

source(file.path(WORKFLOW, "03_methylation", "mm39_annotation.R"))  # mm39 / GENCODE vM32 annotation (since 2026-10-07)

mm285_gr <- mm285_gr_mm39()
mm285_gr <- mm285_gr[names(mm285_gr) %in% complete]
ov <- findOverlaps(dmrs_filt, mm285_gr)
region <- rowsum(betas_adj[names(mm285_gr)[subjectHits(ov)], ], group = queryHits(ov)) /
  as.vector(table(queryHits(ov)))
idx <- as.integer(rownames(region))
d <- dmrs_filt[idx]
dmr_dt <- data.table(dmr = sprintf("%s:%d-%d", seqnames(d), start(d), end(d)),
                     meandiff = mcols(d)$meandiff, region)
fwrite(dmr_dt, file.path(out_dir, "dmr_heatmap_matrix.csv"))
fwrite(meta[, .(sample = Code, donor = Donor, dose_nM = AZA)], file.path(out_dir, "dmr_heatmap_samples.csv"))
cat(sprintf("Panel B: %d of %d DMRs with complete-case probes x %d samples (%d hypo, %d hyper)\n",
            nrow(dmr_dt), length(dmrs_filt), ncol(region), sum(dmr_dt$meandiff < 0), sum(dmr_dt$meandiff > 0)))

## ---- Panel E: Wang similarity + hclust of the shared terms ----
shared <- fread(file.path(DMGSEA_DIR, "meth_rna_shared_GO_BP.csv"))
sd_bp <- godata(annoDb = "org.Mm.eg.db", ont = "BP", computeIC = FALSE)
sim <- termSim(shared$term_id, shared$term_id, semData = sd_bp, method = "Wang")
hc <- hclust(as.dist(1 - sim), method = "ward.D")
fwrite(as.data.table(sim, keep.rownames = "term_id"), file.path(out_dir, "pathway_termsim_wang.csv"))
fwrite(data.table(merge1 = hc$merge[, 1], merge2 = hc$merge[, 2], height = hc$height),
       file.path(out_dir, "pathway_tree_hclust.csv"))
fwrite(data.table(term_id = hc$labels, order = match(seq_along(hc$labels), hc$order), cluster = cutree(hc, k = 5)),
       file.path(out_dir, "pathway_tree_leaves.csv"))
cat(sprintf("Panel E: %d terms, 5 clusters (sizes %s)\n", nrow(sim), paste(table(cutree(hc, k = 5)), collapse = "/")))

## ---- Panel G: candidate gene sets, donor-corrected ----
cg <- file.path(INT_DIR, "candidate_genes")
correct <- function(wide_file, meta_file, out_file) {
  w <- fread(file.path(cg, wide_file)); mt <- fread(file.path(cg, meta_file))
  m <- as.matrix(w[, mt$sample, with = FALSE]); rownames(m) <- w$symbol
  m <- removeBatchEffect(m, batch = factor(mt$donor), design = model.matrix(~dose_nM, data = mt))
  fwrite(as.data.table(m, keep.rownames = "symbol"), file.path(out_dir, out_file))
  fwrite(mt, file.path(out_dir, sub("corrected", "samples", out_file)))
}
correct("rna_expr_wide.csv", "rna_sample_meta.csv", "geneset_rna_corrected.csv")
correct("meth_beta_wide.csv", "meth_sample_meta.csv", "geneset_meth_corrected.csv")
cat("=== Done ===\n")
