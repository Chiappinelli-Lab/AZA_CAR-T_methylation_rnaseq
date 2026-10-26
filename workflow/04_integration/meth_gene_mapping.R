# Shared probe -> gene mapping for every gene-level methylation value in the
# integration (sourced by 02, 04a, 11). Added 2026-09-24, replacing the
# "all genes within 5kb of a DMR" rule (07_pathway_enrichment.R) so that
# gene-level values use the same assignment as the dmGsea pathway analysis
# (Haydar_Methylation_Array/scripts/11_dmgsea.R):
#   - a probe belongs to a gene if it lies between 1.5kb upstream of the TSS
#     and the gene end (promoter + gene body); a probe can map to several genes
#     (results/dmgsea/gene_probe_table.csv, Entrez -> symbol via org.Mm.eg.db
#     from gene_median_beta_change.csv)
#   - "differentially methylated gene" (is_sig_DMR_gene) = >=1 of the gene's
#     probes inside one of the 507 filtered DMRs (FDR + |meandiff| > 0.10) --
#     dmGsea's dmr_hit / sigProbe definition
#   - the gene-level value is the MEDIAN over the gene's probes, so a long,
#     probe-rich gene isn't favoured over a short one.
suppressPackageStartupMessages(library(data.table))

# DMGSEA_DIR (results/methylation/dmgsea) is set in R/paths.R
if (!exists("DMGSEA_DIR")) source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))

# data.table(Name, symbol, is_sig_DMR_gene): one row per probe-gene pair.
load_meth_gene_probe_map <- function() {
  gpt <- fread(file.path(DMGSEA_DIR, "gene_probe_table.csv"), colClasses = c(entrezid = "character"))
  genes <- fread(file.path(DMGSEA_DIR, "gene_median_beta_change.csv"),
                 colClasses = c(entrezid = "character"))[, .(entrezid, symbol, is_sig_DMR_gene = dmr_hit)]
  stopifnot(!anyNA(genes$symbol), !anyDuplicated(genes$symbol))
  m <- merge(gpt, genes, by = "entrezid")[, .(Name, symbol, is_sig_DMR_gene)]
  stopifnot(nrow(m) == nrow(gpt))
  m
}

# Per-gene, per-sample median beta over the gene's probes (NA probes dropped).
# betas: probe x sample matrix. Returns a gene x sample matrix.
gene_median_beta <- function(betas, map) {
  map <- map[Name %in% rownames(betas)]
  long <- data.table(Name = rep(rownames(betas), ncol(betas)),
                     sample = rep(colnames(betas), each = nrow(betas)),
                     beta = as.vector(betas))[Name %in% map$Name & !is.na(beta)]
  long <- merge(long, map[, .(Name, symbol)], by = "Name", allow.cartesian = TRUE)
  agg <- long[, .(beta = median(beta)), by = .(symbol, sample)]
  w <- dcast(agg, symbol ~ sample, value.var = "beta")
  mat <- as.matrix(w[, -"symbol"]); rownames(mat) <- w$symbol
  mat[, colnames(betas), drop = FALSE]
}
