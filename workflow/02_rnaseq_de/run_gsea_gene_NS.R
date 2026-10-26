#!/usr/bin/env Rscript
# GSEA (GO:BP) on the NS gene-level AZA dose-response model -- a
# cutoff-free alternative to the over-representation analysis in
# run_enrichment_gene_NS.R (added 2026-09-22 so pathway results don't depend
# on an arbitrary padj/log2FC threshold).
#
# Ranking metric: DESeq2 Wald statistic for dose_100nM (unshrunk LFC / SE),
# pulled from results/gene_NS/dds.rds because results_full.csv only stores the
# apeglm-shrunk LFC. The Wald stat reflects both effect size and certainty;
# ranking by the apeglm LFC instead would collapse the many low-count genes
# shrunk to ~0 into large ties.
# Gene universe = every gene with a non-NA Wald stat (the >=10-in->=3 filtered
# set, minus DESeq2 outlier/all-zero genes); no independent-filtering step.
# GO settings match the ORA script: BP, gene-set size 15-500, BH, simplify()
# at 0.7 (Wang). Output -> results/gene_NS/enrichment_gsea/, including
# GO_BP_{up,down}_simplified.csv (split by NES sign) so
# run_pathway_diverging_bar.R gsea can plot it.
#
# Publication copy: only the dose-response comparison (the pairwise 0v50/0v100/
# 0v200 runs need run_deseq2_gene_NS_pairwise_vs_0.R, which feeds no figure), and
# the KEGG GSEA (gseKEGG, needs live KEGG access, feeds no figure) is dropped.
# The GO:BP run is unchanged.
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(ggplot2)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
cfg <- list(dds = "gene_NS/dds.rds", coef = "dose_100nM", out = "gene_NS/enrichment_gsea",
            label = "dose-response (per 100nM)")
res_dir <- RNA_DE_DIR
dds_path <- file.path(res_dir, cfg$dds)
out_dir <- file.path(res_dir, cfg$out)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Ranked gene list: Wald stat for the selected coefficient ----
dds <- readRDS(dds_path)
res <- results(dds, name = cfg$coef)
rank_dt <- data.table(ensembl_id = sub("\\..*$", "", rownames(res)), stat = res$stat)
rank_dt <- rank_dt[!is.na(stat)]
# Collapse any duplicate Ensembl IDs (after version stripping) to the max-|stat| row.
rank_dt <- rank_dt[order(-abs(stat))][!duplicated(ensembl_id)]
setorder(rank_dt, -stat)
fwrite(rank_dt, file.path(out_dir, "ranked_gene_list_wald_stat.csv"))
gene_list <- setNames(rank_dt$stat, rank_dt$ensembl_id)
message(sprintf("Ranked genes: %d (stat range %.1f to %.1f)", length(gene_list), min(gene_list), max(gene_list)))

## ---- 2. GO:BP GSEA ----
# Run with no cutoff so every tested term's NES/p is kept (GO_BP_gsea_all_terms.csv,
# used for the RNA-vs-methylation scatter in Haydar_Methylation_Array/scripts/
# 12_dmgsea_plots.R; added 2026-09-23), then apply the same p.adjust <= 0.05 filter
# gseGO's pvalueCutoff would have applied (BH is over all tested terms either way).
gse_go <- gseGO(geneList = gene_list, ont = "BP", OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH",
                eps = 0, seed = TRUE, verbose = FALSE)
gse_go <- setReadable(gse_go, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
fwrite(as.data.frame(gse_go)[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")],
       file.path(out_dir, "GO_BP_gsea_all_terms.csv"))
gse_go@result <- gse_go@result[gse_go@result$p.adjust <= 0.05, ]
go_full <- as.data.frame(gse_go)
fwrite(go_full, file.path(out_dir, "GO_BP_gsea_full.csv"))

gse_go_simple <- if (nrow(go_full) > 0) simplify(gse_go, cutoff = 0.7, by = "p.adjust", select_fun = min) else gse_go
go_simple <- as.data.frame(gse_go_simple)
fwrite(go_simple, file.path(out_dir, "GO_BP_gsea_simplified.csv"))

# Split by NES sign into the up/down files run_pathway_diverging_bar.R expects
# (core_enrichment = leading-edge genes, renamed geneID to match the ORA files).
split_dt <- as.data.table(go_simple)
if ("core_enrichment" %in% names(split_dt)) split_dt[, geneID := core_enrichment]
fwrite(split_dt[NES > 0], file.path(out_dir, "GO_BP_up_simplified.csv"))
fwrite(split_dt[NES < 0], file.path(out_dir, "GO_BP_down_simplified.csv"))

if (nrow(go_simple) > 0) {
  p <- dotplot(gse_go_simple, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign) +
    ggtitle(sprintf("GO:BP GSEA — NS AZA %s (ranked by Wald stat)", cfg$label))
  ggsave(file.path(out_dir, "GO_BP_gsea_dotplot.png"), p, width = 13,
         height = max(5, 0.35 * min(30, nrow(go_simple)) + 2), dpi = 300, bg = "white")
}

## ---- 3. Summary ----
n_sign <- function(df, s) if (nrow(df) == 0) 0L else sum(sign(df$NES) == s)
summary_lines <- c(
  sprintf("GSEA, NS AZA %s (ranked by DESeq2 Wald stat, %s)", cfg$label, cfg$coef),
  sprintf("Ranked genes: %d", length(gene_list)),
  sprintf("GO BP terms (raw): %d (up/positive NES=%d, down/negative NES=%d)",
          nrow(go_full), n_sign(go_full, 1), n_sign(go_full, -1)),
  sprintf("GO BP terms (simplified, cutoff 0.7): %d (up=%d, down=%d)",
          nrow(go_simple), n_sign(go_simple, 1), n_sign(go_simple, -1))
)
writeLines(summary_lines, file.path(out_dir, "summary.txt"))
message(paste(summary_lines, collapse = "\n"))
message("=== Done ===")
