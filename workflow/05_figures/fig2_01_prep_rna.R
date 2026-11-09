#!/usr/bin/env Rscript
# Paper figure (RNA-seq x methylation integration), RNA-side plot-ready data.

#
# Writes to paper_figures/fig2_rna_methylation/data/:
#   de_heatmap_matrix.csv / de_heatmap_samples.csv -- panel A: the 111 DE genes
#     (results_full.csv sig = padj<0.05 & |log2FC|>0.5), VST (blind), donor
#     effect removed with limma::removeBatchEffect(batch = donor, design =
#     ~dose_nM) -- dose effect kept; the design matters here because EV19/EV20
#     have no 200nM sample. Rows keep log2FC for the up/down side bar.
#   pathway_rna_counts.csv -- panels E/F: per shared GO:BP term, number of DE
#     genes in the term (GOALL, i.e. incl. descendant terms, as gseGO uses;
#     restricted to the tested genes) and GSEA leading-edge size.
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(limma)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
rna_dir <- RNA_NS_DIR
out_dir <- file.path(FIG_ROOT, "fig2_rna_methylation", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Panel A: DE genes, donor-corrected VST ----
dds <- readRDS(file.path(rna_dir, "dds.rds"))
res <- fread(file.path(rna_dir, "results_full.csv"))
cd <- as.data.frame(colData(dds))
vst_mat <- assay(vst(dds, blind = TRUE))
vst_corr <- removeBatchEffect(vst_mat, batch = cd$donor, design = model.matrix(~dose_nM, data = cd))

de <- res[sig == TRUE]
m <- vst_corr[de$feature, , drop = FALSE]
lab <- ifelse(is.na(de$symbol) | de$symbol == "", de$ensembl_id, de$symbol)
out <- data.table(symbol = make.unique(lab), log2FoldChange = de$log2FoldChange, padj = de$padj, m)
fwrite(out, file.path(out_dir, "de_heatmap_matrix.csv"))
fwrite(data.table(sample = rownames(cd), donor = cd$donor, dose_nM = cd$dose_nM),
       file.path(out_dir, "de_heatmap_samples.csv"))
cat(sprintf("Panel A: %d DE genes (%d up, %d down) x %d samples\n",
            nrow(out), sum(de$log2FoldChange > 0), sum(de$log2FoldChange < 0), ncol(m)))

## ---- Panels E/F: DE genes and leading edge per shared term ----
shared <- fread(file.path(DMGSEA_DIR, "meth_rna_shared_GO_BP.csv"))
gsea <- fread(file.path(rna_dir, "enrichment_gsea/GO_BP_gsea_full.csv"))
tested <- res[!is.na(padj), ensembl_id]
go2ens <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = shared$term_id, keytype = "GOALL",
                                                 columns = "ENSEMBL"))
go2ens <- unique(as.data.table(go2ens)[!is.na(ENSEMBL) & ENSEMBL %in% tested, .(term_id = GOALL, ENSEMBL)])
cnt <- go2ens[, .(rna_n_tested = uniqueN(ENSEMBL), rna_n_DE = sum(unique(ENSEMBL) %in% de$ensembl_id)), by = term_id]
le <- gsea[ID %in% shared$term_id, .(term_id = ID, rna_n_leading_edge = lengths(strsplit(core_enrichment, "/")))]
cnt <- merge(merge(shared[, .(term_id)], cnt, by = "term_id", all.x = TRUE), le, by = "term_id", all.x = TRUE)
fwrite(cnt, file.path(out_dir, "pathway_rna_counts.csv"))
cat(sprintf("Panels E/F: %d shared terms; DE genes per term %d-%d (median %g)\n",
            nrow(cnt), min(cnt$rna_n_DE), max(cnt$rna_n_DE), median(cnt$rna_n_DE)))
cat("=== Done ===\n")
