#!/usr/bin/env Rscript
# GSK-3484862 vs 5-azacytidine dose response, B7H3 CAR arm (added 2026-10-19). The GSK arm is
# one donor (EV03 = figure donor 8), so both parts compare that donor's dose response with the AZA response;
# neither generalizes to GSK (agreed with the user 2026-10-19).
#
# Part 1, descriptive (no test), at the 507 AZA DMRs (06_dmrcate.R, |meandiff| > 0.10):
#   region mean beta per sample (probes complete in all 16 B7H3 samples, as in figS panel B); per donor an OLS
#   slope of region beta on dose over its 4 samples (0/50/100/200 nM), expressed per 100 nM. Count DMRs where
#   the GSK donor's slope is higher than the mean of the three AZA donors'
#   slopes (user's definition), and, stricter, higher than every AZA donor's slope.
#
# Part 2, model: limma on M-values, 16 B7H3 samples (3 AZA donors + the GSK donor, 4 doses each), autosomal
#   probes complete in all 16 (sex is a donor property: EV03/ZY03 male, EV07/EV08 female). Design
#   ~ Donor + dose_AZA + dose_GSK (= donor + dose:drug; each drug's dose, per 100 nM, is 0 in the other drug's
#   samples), so donor absorbs sex, baseline and chip, and there are no drug or sex main effects (not
#   estimable: both are functions of donor). Contrast dose_GSK - dose_AZA, 10 residual df. The GSK slope's
#   uncertainty is within-donor array noise only (between-donor variation in response can't be estimated with
#   one GSK donor), so p-values describe donor 8 vs the AZA response, not GSK. DMRcate as in 06 (lambda 1000,
#   C 2, FDR, min.cpgs 2), diff = beta-scale slope difference over 0 -> 200 nM, then |meandiff| > 0.10.
#
# Outputs: results/methylation/gsk_vs_aza/
#   part1_dmr507_slopes.csv, part1_summary.txt
#   part2_probe_stats.rds, part2_dmrs_all.csv, part2_dmrs_filtered.csv, part2_summary.txt
# Run by run_all.sh as a SLURM job (env haydar_methylation)
suppressPackageStartupMessages({
  library(limma)
  library(DMRcate)
  library(sesame)  # BetaValueToMValue()
  library(GenomicRanges)
  library(data.table)
})

set.seed(1)

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- GSK_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(WORKFLOW, "03_methylation", "mm39_annotation.R"))
mm285_gr <- mm285_gr_mm39()

ss <- fread(METH_SHEET)[CAR == "B7H3"]
ss[, dose_nM := as.numeric(sub("nM", "", Dose))]
ss[, drug := fifelse(Drug == "GSK-3484862", "GSK", "AZA")]
setorder(ss, drug, Donor, dose_nM)
stopifnot(nrow(ss) == 16, all(ss[, .N, by = Donor]$N == 4))
betas <- readRDS(BETAS_RDS)[, ss$Code]

## ---------------------------------------------------------------- Part 1 ----
res06 <- readRDS(DMRCATE_RDS)
dmrs507 <- res06$dmrs[abs(mcols(res06$dmrs)$meandiff) > 0.10]
b16 <- betas[complete.cases(betas) & rownames(betas) %in% names(mm285_gr), ]
gr16 <- mm285_gr[rownames(b16)]
ov <- findOverlaps(dmrs507, gr16)
region <- rowsum(b16[names(gr16)[subjectHits(ov)], ], group = queryHits(ov)) / as.vector(table(queryHits(ov)))
d507 <- dmrs507[as.integer(rownames(region))]

donor_slope <- function(d) {  # OLS slope of region beta on dose, per 100 nM, for one donor's 4 samples
  cols <- ss[Donor == d, Code]
  x <- ss[Donor == d, dose_nM] / 100
  xc <- x - mean(x)
  as.vector(region[, cols] %*% xc) / sum(xc^2)
}
aza_donors <- ss[drug == "AZA", unique(Donor)]
slopes <- sapply(c(aza_donors, "EV03"), donor_slope)
aza_mean <- rowMeans(slopes[, aza_donors])
gsk <- slopes[, "EV03"]
p1 <- data.table(dmr = sprintf("%s:%d-%d", seqnames(d507), start(d507), end(d507)),
                 meandiff_06 = mcols(d507)$meandiff, n_probes = as.vector(table(queryHits(ov))),
                 slope_per100nM_EV07 = slopes[, "EV07"], slope_per100nM_EV08 = slopes[, "EV08"],
                 slope_per100nM_ZY03 = slopes[, "ZY03"], slope_per100nM_AZA_mean = aza_mean,
                 slope_per100nM_GSK_EV03 = gsk,
                 gsk_above_aza_mean = gsk > aza_mean, gsk_above_all_aza = gsk > apply(slopes[, aza_donors], 1, max))
fwrite(p1, file.path(out_dir, "part1_dmr507_slopes.csv"))
p1_lines <- c(
  "Part 1: GSK donor (EV03, figure donor 8) vs AZA donors at the 507 AZA DMRs (descriptive, no test)",
  "Slope = OLS slope of region-mean beta on dose over each donor's 4 samples, delta beta per 100 nM",
  sprintf("DMRs with region beta (probes complete in all 16 B7H3 samples): %d of %d", nrow(p1), length(dmrs507)),
  sprintf("GSK slope > mean AZA slope: %d of %d", sum(p1$gsk_above_aza_mean), nrow(p1)),
  sprintf("GSK slope > every AZA donor's slope: %d of %d", sum(p1$gsk_above_all_aza), nrow(p1)),
  sprintf("Median slope per 100 nM: AZA mean %.4f (EV07 %.4f, EV08 %.4f, ZY03 %.4f); GSK %.4f",
          median(aza_mean), median(slopes[, "EV07"]), median(slopes[, "EV08"]), median(slopes[, "ZY03"]), median(gsk)),
  sprintf("GSK slope range: %.4f to %.4f; GSK slopes < -0.02 per 100 nM: %d", min(gsk), max(gsk), sum(gsk < -0.02)))
writeLines(p1_lines, file.path(out_dir, "part1_summary.txt"))
cat(p1_lines, sep = "\n")

## ---------------------------------------------------------------- Part 2 ----
auto <- names(mm285_gr)[!(as.character(seqnames(mm285_gr)) %in% c("chrX", "chrY", "chrM"))]
b <- betas[complete.cases(betas) & rownames(betas) %in% auto, ]
m <- BetaValueToMValue(b)
meta <- data.frame(Donor = factor(ss$Donor), dose_AZA = ifelse(ss$drug == "AZA", ss$dose_nM / 100, 0),
                   dose_GSK = ifelse(ss$drug == "GSK", ss$dose_nM / 100, 0))
design <- model.matrix(~ Donor + dose_AZA + dose_GSK, data = meta)
stopifnot(qr(design)$rank == ncol(design))
cm <- makeContrasts(GSK_vs_AZA = dose_GSK - dose_AZA, levels = design)
fit <- eBayes(contrasts.fit(lmFit(m, design), cm))
tt <- topTable(fit, coef = "GSK_vs_AZA", number = Inf, sort.by = "none")
fit_b <- lmFit(b, design)
slope_diff_beta <- (coef(fit_b)[, "dose_GSK"] - coef(fit_b)[, "dose_AZA"]) * 2  # beta, over 0 -> 200 nM
saveRDS(list(tt = tt, slope_AZA_beta_per100 = coef(fit_b)[, "dose_AZA"],
             slope_GSK_beta_per100 = coef(fit_b)[, "dose_GSK"], design = design),
        file.path(out_dir, "part2_probe_stats.rds"))
n_sig <- sum(tt$adj.P.Val < 0.05)
cat(sprintf("\nPart 2: %d autosomal probes complete in 16; residual df %d, prior df %.2f; FDR < 0.05: %d (GSK slope higher %d, lower %d)\n",
            nrow(tt), fit$df.residual[1], fit$df.prior, n_sig,
            sum(tt$adj.P.Val < 0.05 & tt$logFC > 0), sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)))

gr <- mm285_gr[rownames(tt)]
annotated <- GRanges(seqnames = seqnames(gr), ranges = ranges(gr), stat = tt$t,
                     diff = slope_diff_beta[rownames(tt)], rawpval = tt$P.Value,
                     ind.fdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < 0.05)
names(annotated) <- rownames(tt)
dmr <- dmrcate(new("CpGannotated", ranges = annotated), lambda = 1000, C = 2, consec = FALSE,
               pcutoff = "fdr", min.cpgs = 2)
dmrs <- extractRanges(dmr, genome = "mm10")  # gene overlaps are mm10 and dropped (coordinates are mm39; see 06)
df <- as.data.frame(dmrs)
df$overlapping.genes <- NULL
fwrite(df, file.path(out_dir, "part2_dmrs_all.csv"))
filt <- df[abs(df$meandiff) > 0.10, ]
fwrite(filt, file.path(out_dir, "part2_dmrs_filtered.csv"))
filt_gr <- GRanges(filt$seqnames, IRanges(filt$start, filt$end), meandiff = filt$meandiff)
n_ov507 <- sum(overlapsAny(filt_gr, dmrs507))
n_507_hit <- sum(overlapsAny(dmrs507, filt_gr))
p2_lines <- c(
  "Part 2: limma ~ Donor + dose_AZA + dose_GSK on M-values, contrast dose_GSK - dose_AZA (per 100 nM)",
  "16 B7H3 samples (AZA EV07/EV08/ZY03 + GSK EV03, 4 doses each); autosomal probes complete in all 16",
  "Exploratory: compares one GSK donor's dose response with the AZA response; not generalizable to GSK",
  sprintf("Probes tested: %d; residual df %d; eBayes prior df %.2f", nrow(tt), fit$df.residual[1], fit$df.prior),
  sprintf("Probes at FDR < 0.05: %d (GSK slope higher than AZA: %d; lower: %d)", n_sig,
          sum(tt$adj.P.Val < 0.05 & tt$logFC > 0), sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)),
  sprintf("DMRcate DMRs: %d; |meandiff| > 0.10: %d (GSK higher: %d, GSK lower: %d)", nrow(df), nrow(filt),
          sum(filt$meandiff > 0), sum(filt$meandiff < 0)),
  sprintf("Filtered DMRs overlapping the 507 AZA DMRs: %d of %d; AZA DMRs hit: %d of %d", n_ov507, nrow(filt),
          n_507_hit, length(dmrs507)),
  "meandiff = beta-scale (GSK slope - AZA slope) over 0 -> 200 nM; positive = GSK slope higher than AZA slope")
writeLines(p2_lines, file.path(out_dir, "part2_summary.txt"))
cat(p2_lines, sep = "\n")
cat("=== Done ===\n")
