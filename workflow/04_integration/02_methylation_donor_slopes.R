#!/usr/bin/env Rscript
# Per-donor and average-across-donor AZA dose-response slopes, DNA
# methylation, B7H3 CAR arm (the arm the primary DMR analysis was run on).
# Feeds the methylation side of the RNA-seq <-> methylation per-gene
# correlation (see integration/scripts/03).
#
# Mirrors 01_rna_donor_slopes.R's donor-scope logic exactly (same
# closed-form OLS slope, fit independently per donor, then either taken for
# EV07 alone or averaged across all donors) -- fit per PROBE on beta values,
# then summarized per gene as the median probe slope over the gene's promoter
# (1.5kb upstream) + gene-body probes. Median of probe slopes rather than slope
# of the median beta: the latter tracks whichever probe sits in the middle of
# the gene's baseline beta range and correlated only rho=0.70 with dmGsea's
# per-gene median probe change (gene_median_beta_change.csv). This version's
# meth_slope_avg_donors x 200 equals that column exactly, so it is the same gene value shown in 13_meth_rna_cnet.R.
#
# Gene assignment set to the dmGsea probe->gene mapping (meth_gene_mapping.R,
# same as Haydar_Methylation_Array/scripts/11_dmgsea.R), so gene-level values
# match the pathway analysis and don't scale with gene length/probe count.
# Every gene with >=1 mapped probe gets a value (not only genes near a DMR);
# is_sig_DMR_gene = >=1 of its probes inside one of the 507 filtered DMRs.


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
source(file.path(WORKFLOW, "04_integration", "meth_gene_mapping.R"))

out_dir <- INT_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Probe -> gene map (promoter + body, dmGsea) ----
gmap <- load_meth_gene_probe_map()
cat(sprintf("Probe-gene pairs: %d; genes: %d (%d with a probe in a filtered DMR)\n",
            nrow(gmap), uniqueN(gmap$symbol), uniqueN(gmap[is_sig_DMR_gene == TRUE]$symbol)))

## ---- 2. Load beta matrix, subset to AZA + B7H3 CAR arm (12 samples) ----
sample_sheet <- read.csv(METH_SHEET, stringsAsFactors = FALSE)
meta <- sample_sheet %>%
  filter(Drug == "5-Azacytidine", CAR == "B7H3") %>%
  transmute(Code, Donor = factor(Donor), dose_nM = as.numeric(sub("nM", "", Dose)))
stopifnot(nrow(meta) == 12)
cat("\nDonor x dose design (methylation, AZA+B7H3):\n"); print(table(meta$Donor, meta$dose_nM))

betas <- readRDS(BETAS_RDS)
betas_sub <- betas[intersect(unique(gmap$Name), rownames(betas)), meta$Code]
rm(betas); gc()

## ---- 3. Per-probe, per-donor slope (same closed-form OLS as script 01) ----
fit_donor_slope <- function(y_mat, x) {
  xc <- x - mean(x)
  as.vector(y_mat %*% xc) / sum(xc^2)
}

donors <- levels(meta$Donor)
slope_mat <- matrix(NA_real_, nrow = nrow(betas_sub), ncol = length(donors),
                    dimnames = list(rownames(betas_sub), donors))
for (d in donors) {
  samp <- meta$Code[meta$Donor == d]
  slope_mat[, d] <- fit_donor_slope(betas_sub[, samp, drop = FALSE], meta$dose_nM[meta$Donor == d])
}
probe_dt <- data.table(Name = rownames(slope_mat), slope_EV07 = slope_mat[, "EV07"],
                       slope_avg = rowMeans(slope_mat))
cat(sprintf("Probes with a slope in all 3 donors: %d / %d (pOOBAH-masked betas -> NA)\n",
            sum(!is.na(probe_dt$slope_avg)), nrow(probe_dt)))

## ---- 4. Gene level: median over the gene's probes (NA probes dropped per column) ----
gene_summary <- merge(gmap, probe_dt, by = "Name")[, .(
  n_probes = sum(!is.na(slope_avg)),
  is_sig_DMR_gene = is_sig_DMR_gene[1],
  meth_slope_EV07 = median(slope_EV07, na.rm = TRUE),
  meth_slope_avg_donors = median(slope_avg, na.rm = TRUE)
), by = symbol]
gene_summary <- gene_summary[!is.na(meth_slope_EV07) & !is.na(meth_slope_avg_donors)]

cat(sprintf("\nFinal methylation gene-slope table: %d genes (%d with a probe in a filtered DMR)\n",
            nrow(gene_summary), sum(gene_summary$is_sig_DMR_gene)))

## ---- 5. Save ----
fwrite(gene_summary, file.path(out_dir, "methylation_gene_slopes.csv"))
cat("Saved", file.path(out_dir, "methylation_gene_slopes.csv"), "\n")
cat("=== Done ===\n")
