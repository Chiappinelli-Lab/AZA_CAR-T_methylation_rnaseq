#!/usr/bin/env Rscript
# Gene-set enrichment of the primary DMR hits with dmGsea::gsProbe (Xu,
# Motsinger-Reif & Niu 2025, Bioinformatics btaf422), run in goregion style.
# Added 2026-09-23; choices confirmed with the user:
#
# Tool choice -- dmGsea over missMethyl gometh/goregion (Maksimovic et al.
#   2021). missMethyl 1.44 maps every probe's gene symbol to Entrez via
#   org.Hs.eg.db + limma::alias2SymbolTable (human) inside .getFlatAnnotation(),
#   even with a custom `anno`, and truncates probe IDs to 10 characters
#   (dropping MM285's _TC21/_BC21 suffixes) -- see also 07_pathway_enrichment.R.
#   dmGsea takes a user-supplied probe->Entrez table and builds mouse GO/KEGG
#   sets (species = "Mouse").
#
# Test choice -- gsProbe with sigProbe = probes inside the 507 filtered DMRs
#   (FDR + |meandiff| > 0.10, 06_dmrcate.R), allProbe = every tested probe that
#   maps to a gene. This is goregion's design (the other top method in
#   Maksimovic et al.): over-representation of DMR genes, with each gene's
#   probability of being "hit" weighted by its probe count (noncentral
#   hypergeometric). The threshold-based alternatives were tried first and
#   rejected (2026-09-23): 38% of probes are FDR<0.05 (99.7% hypomethylated),
#   so gsPG called ~89% of genes significant (3,523/5,222 GO:BP sets FDR<0.05),
#   gsProbe at probe FDR<0.05 would call 78.5%, and gsGene ranking was
#   dominated by long, probe-rich neuronal/developmental genes. With the 507-DMR
#   probes (1,055 probes) only 1.6% of genes count as hits.
#
# dmGsea bug patched here: in dmGsea 1.0.2 (and devel 1.3.3, R/utils.R
#   p.gs(), lines 182-183) the upper-tail p-value is written as
#       p <- pFNCHypergeo(..., lower.tail = FALSE)
#       +dFNCHypergeo(...)
#   R parses the second line as a separate statement, so p = P(X > k) instead
#   of P(X >= k): p-values too small, and exactly 0 when every gene in a set
#   is a hit. p.gs is replaced in the dmGsea namespace with the one-line fix
#   (sum of both terms) before running; the patch is verified below.
#
# Probe -> gene (GeneProbeTable): promoter + gene body, Illumina-style -- a
#   probe maps to a gene if it lies between 1.5kb upstream of the gene's
#   strand-aware TSS and the gene end (genes_gr gene models from
#   results/dm/tested_gene_background.rds, as in 06f/07). Ensembl -> Entrez via
#   org.Mm.eg.db (first match). A probe can map to several genes.
#
# Gene sets: mouse GO:BP (dmGsea::getGO), restricted to the mapped-gene
#   universe and to 15-500 genes (dmGsea applies no size filter itself; 15-500
#   matches the other GO analyses). Publication copy: the KEGG collection
#   (dmGsea::getKEGG) is dropped -- it needs live KEGGREST access and feeds no
#   figure. The GO:BP run is unchanged.
#
# Output: results/methylation/dmgsea/ -- raw dmGsea CSVs, dmgsea_GO_BP.csv (tidy),
#   gene_probe_table.csv, dmr_hit_genes.csv (per hit gene: most significant DMR's
#   meandiff/FDR), gene_median_beta_change.csv (every mapped gene: median probe
#   beta change 0->200nM), summary.txt. Plots: 12_dmgsea_plots.R.
suppressPackageStartupMessages({
  library(dmGsea)
  library(BiasedUrn)
  library(sesame)
  library(sesameData)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(data.table)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- DMGSEA_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
UPSTREAM_BP <- 1500
MIN_SET <- 15; MAX_SET <- 500
MEANDIFF_CUT <- 0.10
NCORE <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

## ---- 0. Patch dmGsea's p.gs (P(X > k) -> P(X >= k)) ----
p_gs_fixed <- function(gs, gene.info) {
  gs.ng <- length(gs)
  gs.score <- sum(gene.info$score[gene.info$entrezid %in% gs])
  n.sg <- sum(gene.info$score)
  if (gs.score == 0) {
    p <- 1
  } else {
    pw.red <- with(gene.info, mean(pwf[entrezid %in% gs]))
    pw.white <- with(gene.info, mean(pwf[!(entrezid %in% gs)]))
    odds <- min((1 - pw.red) * pw.white / ((1 - pw.white) * pw.red), 1e+100)
    ng2 <- nrow(gene.info) - gs.ng
    p <- BiasedUrn::pFNCHypergeo(gs.score, gs.ng, ng2, n.sg, odds, lower.tail = FALSE) +
      BiasedUrn::dFNCHypergeo(gs.score, gs.ng, ng2, n.sg, odds)
  }
  c(gs.ng, gs.score, p)
}
# Verify the bug exists in the installed version and that the patch fixes it:
# a toy set where every gene is a hit must get p > 0 after patching.
toy <- data.frame(entrezid = as.character(1:100), pwf = 0.5, score = c(rep(1, 20), rep(0, 80)))
p_orig <- get("p.gs", asNamespace("dmGsea"))(as.character(1:5), toy)[3]
p_new <- p_gs_fixed(as.character(1:5), toy)[3]
cat(sprintf("p.gs check (5/5 hits, 20/100 overall): installed=%.3g, patched=%.3g\n", p_orig, p_new))
stopifnot(p_new > 0, p_new >= p_orig)
assignInNamespace("p.gs", p_gs_fixed, ns = "dmGsea")
stopifnot(identical(body(get("p.gs", asNamespace("dmGsea"))), body(p_gs_fixed)))

## ---- 1. Tested probes and the 507 filtered DMRs (06_dmrcate.R) ----
res <- readRDS(DMRCATE_RDS)
tt <- res$tt
dmrs_filt <- res$dmrs[abs(mcols(res$dmrs)$meandiff) > MEANDIFF_CUT]
cat(sprintf("Tested probes: %d | filtered DMRs: %d\n", nrow(tt), length(dmrs_filt)))

## ---- 2. GeneProbeTable: promoter (1.5kb upstream) + gene body ----
source(file.path(WORKFLOW, "03_methylation", "mm39_annotation.R"))  # mm39 / GENCODE vM32 annotation (since 2026-10-07)
mm285_gr <- mm285_gr_mm39()
probe_gr <- mm285_gr[intersect(rownames(tt), names(mm285_gr))]
genes_gr <- readRDS(file.path(METH_DM_DIR, "tested_gene_background.rds"))$genes_gr
window_gr <- punion(genes_gr, promoters(genes_gr, upstream = UPSTREAM_BP, downstream = 1))
names(window_gr) <- names(genes_gr)
common_levels <- union(union(seqlevels(probe_gr), seqlevels(window_gr)), seqlevels(dmrs_filt))
seqlevels(probe_gr) <- common_levels
seqlevels(window_gr) <- common_levels
seqlevels(dmrs_filt) <- common_levels

ov <- findOverlaps(probe_gr, window_gr, ignore.strand = TRUE)
ens <- sub("\\..*$", "", names(window_gr)[subjectHits(ov)])
entrez <- mapIds(org.Mm.eg.db, keys = unique(ens), column = "ENTREZID",
                 keytype = "ENSEMBL", multiVals = "first")
gpt <- unique(data.table(Name = names(probe_gr)[queryHits(ov)], entrezid = unname(entrez[ens]))[!is.na(entrezid)])
stopifnot(all(gpt$Name %in% rownames(tt)))
fwrite(gpt, file.path(out_dir, "gene_probe_table.csv"))

## ---- 3. sigProbe = tested probes inside the filtered DMRs ----
all_probe <- unique(gpt$Name)
dmr_probe <- names(probe_gr)[unique(queryHits(findOverlaps(probe_gr, dmrs_filt, ignore.strand = TRUE)))]
sig_probe <- intersect(dmr_probe, all_probe)
ppg <- gpt[, .N, by = entrezid]$N
stats_lines <- c(
  sprintf("Probes mapped to >=1 gene (promoter 1.5kb + body): %d of %d (%.1f%%)",
          length(all_probe), nrow(tt), 100 * length(all_probe) / nrow(tt)),
  sprintf("Genes (Entrez) with >=1 probe: %d; probes per gene median %g, mean %.1f, max %d",
          uniqueN(gpt$entrezid), median(ppg), mean(ppg), max(ppg)),
  sprintf("Probes in the %d filtered DMRs: %d; of those mapped to a gene (sigProbe): %d",
          length(dmrs_filt), length(dmr_probe), length(sig_probe)),
  sprintf("Genes hit (>=1 sigProbe): %d (%.1f%% of mapped genes)",
          uniqueN(gpt[Name %in% sig_probe]$entrezid), 100 * uniqueN(gpt[Name %in% sig_probe]$entrezid) / uniqueN(gpt$entrezid)))
cat(stats_lines, sep = "\n")

## ---- 3b. Per-gene DMR summary for the hit genes (used by the gene x pathway
## heatmap in 12_dmgsea_plots.R): each hit gene's most significant DMR (lowest
## min_smoothed_fdr) among the filtered DMRs its sigProbes fall in.
ov_dmr <- findOverlaps(probe_gr[sig_probe], dmrs_filt, ignore.strand = TRUE)
probe_dmr <- data.table(Name = sig_probe[queryHits(ov_dmr)], dmr_idx = subjectHits(ov_dmr),
                        meandiff = mcols(dmrs_filt)$meandiff[subjectHits(ov_dmr)],
                        dmr_fdr = mcols(dmrs_filt)$min_smoothed_fdr[subjectHits(ov_dmr)])
gene_dmr <- merge(gpt, probe_dmr, by = "Name")
gene_dmr <- gene_dmr[order(entrezid, dmr_fdr)][, .(n_dmrs = uniqueN(dmr_idx), n_dmr_probes = uniqueN(Name),
                                                   top_dmr_meandiff = meandiff[1], top_dmr_fdr = dmr_fdr[1],
                                                   mean_meandiff = mean(unique(data.table(dmr_idx, meandiff))$meandiff)),
                                               by = entrezid]
gene_dmr[, symbol := mapIds(org.Mm.eg.db, keys = entrezid, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")]
setcolorder(gene_dmr, c("entrezid", "symbol"))
fwrite(gene_dmr[order(top_dmr_fdr)], file.path(out_dir, "dmr_hit_genes.csv"))

## ---- 3c. Gene-level methylation change for every mapped gene (2026-09-23):
## median over the gene's probes (same promoter 1.5kb + body assignment as the
## GeneProbeTable) of the per-probe fitted beta change over 0->200nM
## (myannotation `diff` from 06_dmrcate.R). Used to color every gene in the
## methylation panel of 13_meth_rna_cnet.R, including genes without a DMR.
probe_diff <- data.table(Name = names(res$myannotation@ranges), diff = res$myannotation@ranges$diff,
                         fdr = res$myannotation@ranges$ind.fdr)
gene_meth <- merge(gpt, probe_diff, by = "Name")[, .(n_probes = .N, n_probes_fdr05 = sum(fdr < 0.05),
                                                     median_beta_change = median(diff)), by = entrezid]
gene_meth[, symbol := mapIds(org.Mm.eg.db, keys = entrezid, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")]
gene_meth[, dmr_hit := entrezid %in% gene_dmr$entrezid]
setcolorder(gene_meth, c("entrezid", "symbol"))
fwrite(gene_meth[order(median_beta_change)], file.path(out_dir, "gene_median_beta_change.csv"))

## ---- 4. Gene sets (mouse), restricted to universe and 15-500 genes ----
universe <- unique(gpt$entrezid)
size_filter <- function(gs) {
  gs <- lapply(gs, function(x) intersect(unique(as.character(x)), universe))
  gs[lengths(gs) >= MIN_SET & lengths(gs) <= MAX_SET]
}
gene_sets <- list(GO_BP = size_filter(getGO(subset = "BP", species = "Mouse")))
for (nm in names(gene_sets)) cat(sprintf("%s sets after %d-%d size filter: %d\n", nm, MIN_SET, MAX_SET, length(gene_sets[[nm]])))

## ---- 5. gsProbe per collection ----
summary_lines <- c("dmGsea gsProbe (goregion-style), MM285 AZA dose-response, 507 filtered DMRs as sigProbe",
                   "p.gs patched to P(X >= k) (see header)", stats_lines)
for (nm in names(gene_sets)) {
  prefix <- file.path(out_dir, paste0("gsProbe_", nm))
  message("== gsProbe ", nm)
  gsProbe(probe.p = NULL, sigProbe = sig_probe, allProbe = all_probe,
          GeneProbeTable = as.data.frame(gpt), geneSet = gene_sets[[nm]], species = "Mouse",
          outfile = prefix, ncore = NCORE)
  d <- fread(paste0(prefix, "_userSet_userSet.csv"))
  d[, c("term_id", "term") := tstrsplit(ID, "|", fixed = TRUE, keep = 1:2)]
  d[is.na(term), term := term_id]
  d[, term := sub(" - Mus musculus \\(house mouse\\)$", "", term)]
  d <- d[, .(term_id, term, nGene, nSigGene, p, FDR, Sig.Gene)]
  setorder(d, p)
  fwrite(d, file.path(out_dir, sprintf("dmgsea_%s.csv", nm)))
  summary_lines <- c(summary_lines, sprintf("%s: %d sets tested | FDR<0.05: %d | FDR<0.10: %d | min p %.2g",
                                            nm, nrow(d), sum(d$FDR < 0.05), sum(d$FDR < 0.10), min(d$p)))
}
writeLines(summary_lines, file.path(out_dir, "summary.txt"))
cat(summary_lines, sep = "\n")
cat("=== Done ===\n")
