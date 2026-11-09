#!/usr/bin/env Rscript
# Figure 3 cytokine panel: methylation side. Env `haydar_methylation`.
# Unstimulated 5-azacytidine B7H3 CAR samples (3 donors x 4 doses). Gene value = median beta over the gene's
# promoter (1.5 kb upstream) + gene-body probes (integration/scripts/meth_gene_mapping.R, the dmGsea mapping
# used for every gene-level methylation value). Stats: lm(gene beta ~ donor + dose_nM) -- the same design as the
# DMR analysis, on the gene-level value -- reported as the fitted beta change over 0->200 nM.
suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- file.path(FIG_ROOT, "fig3_cytokines", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(WORKFLOW, "04_integration", "meth_gene_mapping.R"))
GENES <- c("Il2", "Ifng", "Tnf")

ss <- fread(METH_SHEET)
meta <- ss[Drug == "5-Azacytidine" & CAR == "B7H3", .(sample = Code, donor = Donor, dose_nM = as.numeric(sub("nM", "", Dose)))]
stopifnot(nrow(meta) == 12)
gmap <- load_meth_gene_probe_map()[symbol %in% GENES]
betas <- readRDS(BETAS_RDS)[unique(gmap$Name), meta$sample, drop = FALSE]
gb <- gene_median_beta(betas, gmap)

long <- rbindlist(lapply(GENES, function(g) data.table(symbol = g, meta, beta = gb[g, meta$sample])))
fwrite(long, file.path(out_dir, "meth_unstimulated_long.csv"))
stats <- rbindlist(lapply(GENES, function(g) {
  f <- summary(lm(beta ~ donor + dose_nM, data = long[symbol == g]))$coefficients["dose_nM", ]
  data.table(symbol = g, n_probes = sum(gmap$symbol == g), delta_beta_0to200 = 200 * f[["Estimate"]],
             pvalue = f[["Pr(>|t|)"]], any_probe_in_507_DMRs = gmap[symbol == g, is_sig_DMR_gene[1]])
}))
fwrite(stats, file.path(out_dir, "meth_unstimulated_stats.csv"))
# per-probe values too (for a supplementary look at promoter vs body probes)
fwrite(as.data.table(betas, keep.rownames = "probe")[, symbol := gmap$symbol[match(probe, gmap$Name)]],
       file.path(out_dir, "meth_unstimulated_probe_betas.csv"))
print(stats)

## ---- type I IFNs-> meth_unstimulated_{long,stats}_ifn.csv ----
# Ifnb1, and the Ifna family pooled as one "Ifna" gene (median over all Ifna1-16 + Ifnab probes; pseudogene
# Ifna-ps1 excluded) to match the pooled RNA value.
gmap_ifn <- load_meth_gene_probe_map()[symbol == "Ifnb1" | grepl("^Ifna([0-9]+|b)$", symbol)]
gmap_ifn[symbol != "Ifnb1", symbol := "Ifna"]
gmap_ifn <- unique(gmap_ifn, by = c("Name", "symbol"))
IFN <- c("Ifnb1", "Ifna")
betas_ifn <- readRDS(BETAS_RDS)[unique(gmap_ifn$Name), meta$sample, drop = FALSE]
gb_ifn <- gene_median_beta(betas_ifn, gmap_ifn)
long_ifn <- rbindlist(lapply(IFN, function(g) data.table(symbol = g, meta, beta = gb_ifn[g, meta$sample])))
fwrite(long_ifn, file.path(out_dir, "meth_unstimulated_long_ifn.csv"))
stats_ifn <- rbindlist(lapply(IFN, function(g) {
  f <- summary(lm(beta ~ donor + dose_nM, data = long_ifn[symbol == g]))$coefficients["dose_nM", ]
  data.table(symbol = g, n_probes = sum(gmap_ifn$symbol == g), delta_beta_0to200 = 200 * f[["Estimate"]],
             pvalue = f[["Pr(>|t|)"]], any_probe_in_507_DMRs = any(gmap_ifn[symbol == g, is_sig_DMR_gene]))
}))
fwrite(stats_ifn, file.path(out_dir, "meth_unstimulated_stats_ifn.csv"))
print(stats_ifn)
