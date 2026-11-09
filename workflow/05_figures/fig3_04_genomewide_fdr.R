#!/usr/bin/env Rscript
# Do the Fig 3 cytokine / type I IFN methylation changes hold up genome-wide with BH correction? 
#  Two views, same 12 samples (5-Aza, B7H3 CAR, 3 donors x 4 doses):
#  1. gene level -- the Fig 3 test (median beta over promoter + body probes, beta ~ donor + dose) run on EVERY
#     mapped gene, BH across genes. Ordinary OLS t (exactly the Fig 3 p-values) and limma-moderated t.
#     The pooled Ifna family is added as one extra unit (as in the figure).
#  2. probe level -- the genome-wide limma fit from 06_dmrcate.R (M-values, ~ Donor + AZA), BH across all
#     tested probes; the genes' own probes looked up.
suppressPackageStartupMessages({ library(data.table); library(limma) })
source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- file.path(FIG_ROOT, "fig3_cytokines", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(WORKFLOW, "04_integration", "meth_gene_mapping.R"))
GENES <- c("Il2", "Ifng", "Tnf", "Ifnb1", "Ifna")

ss <- fread(METH_SHEET)
meta <- ss[Drug == "5-Azacytidine" & CAR == "B7H3", .(sample = Code, donor = factor(Donor),
                                                     dose_nM = as.numeric(sub("nM", "", Dose)))]
betas <- readRDS(BETAS_RDS)[, meta$sample]
gmap <- load_meth_gene_probe_map()
ifna <- gmap[grepl("^Ifna([0-9]+|b)$", symbol)][, symbol := "Ifna"]
gmap2 <- unique(rbind(gmap[, .(Name, symbol)], ifna[, .(Name, symbol)]))
gb <- gene_median_beta(betas[rownames(betas) %in% gmap2$Name, ], gmap2)
gb <- gb[rowSums(is.na(gb)) == 0, ]

## ---- 1. gene level ----
design <- model.matrix(~ donor + dose_nM, data = meta)
fit <- lmFit(gb, design)
t_ols <- fit$coefficients[, "dose_nM"] / (fit$stdev.unscaled[, "dose_nM"] * fit$sigma)
p_ols <- 2 * pt(-abs(t_ols), fit$df.residual)
eb <- eBayes(fit)
gene <- data.table(symbol = rownames(gb), delta_beta_0to200 = 200 * fit$coefficients[, "dose_nM"],
                   p_ols = p_ols, p_moderated = eb$p.value[, "dose_nM"])
gene[symbol != "Ifna", `:=`(fdr_ols = p.adjust(p_ols, "BH"), fdr_moderated = p.adjust(p_moderated, "BH"))]
# pooled Ifna: BH over all genes + this one unit
all_bh <- function(p) p.adjust(p, "BH")[gene$symbol == "Ifna"]
gene[symbol == "Ifna", `:=`(fdr_ols = all_bh(gene$p_ols), fdr_moderated = all_bh(gene$p_moderated))]
gene[, `:=`(rank_ols = frank(p_ols), rank_moderated = frank(p_moderated))]
fwrite(gene, file.path(out_dir, "meth_genewide_dose_tests.csv"))
cat(sprintf("Gene level: %d genes (+ pooled Ifna). FDR < 0.05: %d OLS, %d moderated; %.1f%% of genes lose methylation\n",
            sum(gene$symbol != "Ifna"), sum(gene$fdr_ols < 0.05), sum(gene$fdr_moderated < 0.05),
            100 * mean(gene$delta_beta_0to200 < 0)))
cat(sprintf("Genome-wide median gene delta beta: %.4f; genes with delta beta <= -0.03: %.1f%%\n",
            median(gene$delta_beta_0to200), 100 * mean(gene$delta_beta_0to200 <= -0.03)))
print(gene[symbol %in% c(GENES, grep("^Ifna[0-9]|^Ifnab$", symbol, value = TRUE))][order(p_ols)], digits = 3)

## ---- 2. probe level (06_dmrcate.R fit) ----
tt <- as.data.table(readRDS(DMRCATE_RDS)$tt, keep.rownames = "Name")
tt[, fdr := p.adjust(P.Value, "BH")]
cat(sprintf("\nProbe level: %d probes tested, %.1f%% FDR < 0.05\n", nrow(tt), 100 * mean(tt$fdr < 0.05)))
pm <- unique(rbind(gmap[symbol %in% GENES, .(Name, symbol)], ifna[, .(Name, symbol)]))
probe <- merge(pm, tt[, .(Name, logFC_M_per_nM = logFC, P.Value, fdr)], by = "Name")
fwrite(probe, file.path(out_dir, "meth_probe_dose_tests_fig3_genes.csv"))
print(probe[, .(n_probes = .N, n_fdr05 = sum(fdr < 0.05), min_fdr = signif(min(fdr), 3),
                median_fdr = signif(median(fdr), 3)), by = symbol][match(GENES, symbol)])
