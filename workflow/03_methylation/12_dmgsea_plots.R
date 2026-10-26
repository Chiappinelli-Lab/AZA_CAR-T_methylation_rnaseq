#!/usr/bin/env Rscript
# Publication copy of 12_dmgsea_plots.R, trimmed to its section 3 table:
# meth_rna_shared_GO_BP.csv -- methylation GO:BP terms (dmGsea gsProbe, 11_dmgsea.R,
# FDR<0.05) that are also significant in the RNA-seq dose-response GSEA
# (workflow/02_rnaseq_de/run_gsea_gene_NS.R, GO_BP_gsea_full.csv, padj<0.05), with the
# RNA direction from the NES sign. Read by figures/fig2 (panel F term set) via
# fig2_02_prep_meth.R. Separate cohorts (only EV07 shared) -- pathway-level
# concordance, not a paired analysis.
# Dropped from the original (no figure uses them): the bar charts, the goseq
# comparison (needs 07_pathway_enrichment.R), the UpSet plot, the scatter, the
# pathway x gene heatmap and the GO tree. Run in env deseq2 (as the original).
suppressPackageStartupMessages({
  library(data.table)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
in_dir <- DMGSEA_DIR
rna_gsea <- file.path(RNA_NS_DIR, "enrichment_gsea", "GO_BP_gsea_full.csv")
FDR_CUT <- 0.05

go <- fread(file.path(in_dir, "dmgsea_GO_BP.csv"))

## ---- 3. Methylation vs RNA-seq GSEA (GO:BP) ----
rna <- fread(rna_gsea)[, .(term_id = ID, rna_NES = NES, rna_FDR = p.adjust)]
shared <- merge(go[FDR < FDR_CUT, .(term_id, term, nGene, nSigGene, meth_FDR = FDR, dmr_genes = Sig.Gene)],
                rna[rna_FDR < 0.05], by = "term_id")
shared[, rna_direction := ifelse(rna_NES > 0, "up", "down")]
setorder(shared, rna_direction, meth_FDR)
fwrite(shared, file.path(in_dir, "meth_rna_shared_GO_BP.csv"))
cat(sprintf("GO:BP terms FDR<%g in methylation: %d | shared with RNA GSEA: %d (up %d, down %d)\n",
            FDR_CUT, sum(go$FDR < FDR_CUT), nrow(shared), sum(shared$rna_direction == "up"),
            sum(shared$rna_direction == "down")))
