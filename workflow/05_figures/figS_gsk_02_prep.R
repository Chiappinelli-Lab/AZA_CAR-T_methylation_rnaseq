#!/usr/bin/env Rscript
# GSK-3484862 supplementary figure, plot-ready data (added 2026-10-19; collaborator request to regenerate the GSK
# heatmap in our style, plus GSK vs AZA). Run in env `haydar_methylation`. All betas are the shared sesame QCDPBM matrix
# (results/qc/betas_QCDPBM.rds), so GSK and AZA went through identical processing. B7H3 CAR arm only (STOP
# excluded, as everywhere else).
#
# Writes to paper_figures/figS_gsk_methylation/data/:
#   top500_gsk_{matrix,samples}.csv -- panel A: the 500 most variable CpG (cg) probes across the 4 GSK samples
#     (EV03, 0/50/100/200 nM), beta values; the plot z-scores rows. Same definition as the collaborators'
#     original "Top 500 Variable Methylation Probes" heatmap.
#   dmr507_delta_{matrix,samples}.csv -- panel B: the 507 AZA DMRs (06_dmrcate.R, |meandiff| > 0.10), region
#     mean beta (complete-case probes across the 16 B7H3 samples) minus the same donor's 0 nM value, for the
#     3 AZA donors and the GSK donor. Within-donor delta replaces the donor adjustment of Fig 2B, which can't
#     be estimated for a single GSK donor.
#   top500_delta_{matrix,samples}.csv -- panel C: unsupervised AZA-vs-GSK comparison (user, 2026-10-19): top 500
#     variable CpGs after donor centring. cg probes on autosomes only, complete across the 16 B7H3 samples; each
#     donor's mean over its 4 samples is subtracted, then probes are ranked by variance across all 16. Centring
#     removes the donor/sex/chip baseline: without it, 496 of the top 500 raw-beta probes are on chrX/Y (EV03
#     and ZY03 male, EV07 and EV08 female) and 384 of 500 autosomal ones are donor-dominated. Donor-mean
#     centring rather than subtracting the 0 nM sample, so one noisy untreated array doesn't drive selection.
#     The matrix holds delta beta vs the same donor's 0 nM (the display, matching panel B).
#   gsk_vs_aza_dmr507_slopes.csv -- panel D: copy of results/gsk_vs_aza/part1_dmr507_slopes.csv (per-DMR, per-donor
#     OLS slope of region beta on dose, per 100 nM; Haydar_Methylation_Array/scripts/13_gsk_vs_aza.R part 1).
#   gsk_vs_aza_probe_volcano.csv / gsk_vs_aza_counts.csv -- panel E: per-probe limma contrast dose_GSK - dose_AZA
#     (13_gsk_vs_aza.R part 2): beta-scale slope difference over 0 -> 200 nM and BH FDR; and the DMR counts
#     (|meandiff| > 0.10, by sign, and overlap with the 507 AZA DMRs).
#   gsk_vs_aza_dmr_volcano.csv -- panel F: all DMRcate DMRs on the dose_GSK - dose_AZA contrast
#     (results/gsk_vs_aza/part2_dmrs_all.csv) with flags passes_filter (|meandiff| > 0.10) and overlaps_507.
#   gsk_vs_aza_sig_probe_slopes.csv / gsk_vs_aza_dmr_slopes.csv -- panels G/H: model dose slopes (beta per 100 nM,
#     from lmFit on beta with the 13_gsk_vs_aza.R design) for AZA (pooled over the 3 donors) and GSK (donor 8). G =
#     probes at FDR < 0.05 on the contrast; H = every DMRcate DMR (all FDR-significant), slope = mean of the slopes
#     of the contrast's probes inside the DMR (the way DMRcate averages `diff` into meandiff; identical to meandiff
#     except 10 DMRs with one extra probe, max difference 0.028), with passes_filter.
#   (A single-donor DMRcate panel was tried and dropped: no probe reached FDR < 0.05 for any donor, AZA
#   included; Haydar_Methylation_Array/archive/single_donor_dmrcate_2026-10-19/.)
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- file.path(FIG_ROOT, "figS_gsk_methylation", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(WORKFLOW, "03_methylation", "mm39_annotation.R"))

ss <- fread(METH_SHEET)
ss[, dose_nM := as.numeric(sub("nM", "", Dose))]
ss[, drug := fifelse(Drug == "GSK-3484862", "GSK", "AZA")]
betas <- readRDS(BETAS_RDS)
mm285_gr <- mm285_gr_mm39()

b7 <- ss[CAR == "B7H3"]
setorder(b7, drug, Donor, dose_nM)
stopifnot(nrow(b7) == 16)

## ---- Panel A: top 500 variable probes, GSK donor ----
gsk <- b7[drug == "GSK"][order(dose_nM)]
bg <- betas[grepl("^cg", rownames(betas)), gsk$Code]
bg <- bg[complete.cases(bg) & rownames(bg) %in% names(mm285_gr), ]
v <- apply(bg, 1, var)
top <- names(sort(v, decreasing = TRUE))[1:500]
fwrite(data.table(probe = top, chr = as.character(seqnames(mm285_gr[top])), variance = v[top], bg[top, ]),
       file.path(out_dir, "top500_gsk_matrix.csv"))
fwrite(gsk[, .(sample = Code, drug, donor = Donor, dose_nM)], file.path(out_dir, "top500_gsk_samples.csv"))
cat(sprintf("Panel A: top 500 of %d complete cg probes; variance range %.4g-%.4g; beta range across the 4 samples per probe: median %.3f\n",
            nrow(bg), min(v[top]), max(v[top]), median(apply(bg[top, ], 1, function(x) diff(range(x))))))

## ---- region mean beta helper ----
region_beta <- function(regions, b) {
  gr <- mm285_gr[rownames(b)]
  ov <- findOverlaps(regions, gr)
  m <- rowsum(b[names(gr)[subjectHits(ov)], , drop = FALSE], group = queryHits(ov)) /
    as.vector(table(queryHits(ov)))
  list(m = m, idx = as.integer(rownames(m)))
}

## ---- Panel B: 507 AZA DMRs, delta beta vs own 0 nM, AZA + GSK ----
res <- readRDS(DMRCATE_RDS)
dmrs_filt <- res$dmrs[abs(mcols(res$dmrs)$meandiff) > 0.10]
b16 <- betas[, b7$Code]
b16 <- b16[complete.cases(b16) & rownames(b16) %in% names(mm285_gr), ]
rb <- region_beta(dmrs_filt, b16)
delta <- rb$m
for (d in unique(b7$Donor)) {
  cols <- b7[Donor == d, Code]
  delta[, cols] <- delta[, cols] - delta[, b7[Donor == d & dose_nM == 0, Code]]
}
d <- dmrs_filt[rb$idx]
fwrite(data.table(dmr = sprintf("%s:%d-%d", seqnames(d), start(d), end(d)), meandiff = mcols(d)$meandiff, delta),
       file.path(out_dir, "dmr507_delta_matrix.csv"))
fwrite(b7[, .(sample = Code, drug, donor = Donor, dose_nM)], file.path(out_dir, "dmr507_delta_samples.csv"))
for (dr in c("AZA", "GSK")) {
  cols <- b7[drug == dr & dose_nM == 200, Code]
  cat(sprintf("Panel B: %s 200 nM, mean region delta beta over %d DMRs: %s; DMRs with delta <= -0.10: %s\n",
              dr, nrow(delta), paste(sprintf("%.3f", colMeans(delta[, cols, drop = FALSE])), collapse = ", "),
              paste(colSums(delta[, cols, drop = FALSE] <= -0.10), collapse = ", ")))
}

## ---- Panel C: top 500 variable CpGs after donor centring, AZA + GSK ----
auto <- !(as.character(seqnames(mm285_gr)) %in% c("chrX", "chrY", "chrM"))
cg <- rownames(betas)[grepl("^cg", rownames(betas)) & rownames(betas) %in% names(mm285_gr)[auto]]
bc <- betas[cg, b7$Code]
bc <- bc[complete.cases(bc), ]
dl <- bc
for (d in unique(b7$Donor)) {
  cols <- b7[Donor == d, Code]
  dl[, cols] <- bc[, cols] - bc[, b7[Donor == d & dose_nM == 0, Code]]
}
cen <- bc
for (d in unique(b7$Donor)) {
  cols <- b7[Donor == d, Code]
  cen[, cols] <- bc[, cols] - rowMeans(bc[, cols])
}
v_cen <- apply(cen, 1, var)
top <- names(sort(v_cen, decreasing = TRUE))[1:500]
fwrite(data.table(probe = top, chr = as.character(seqnames(mm285_gr[top])), variance_donor_centred = v_cen[top],
                  dl[top, ]),
       file.path(out_dir, "top500_delta_matrix.csv"))
fwrite(b7[, .(sample = Code, drug, donor = Donor, dose_nM)], file.path(out_dir, "top500_delta_samples.csv"))
treated <- b7[dose_nM > 0]
mean_delta <- sapply(c("AZA", "GSK"), function(dr) mean(dl[top, treated[drug == dr & dose_nM == 200, Code]]))
q <- quantile(abs(dl[top, treated$Code]), c(0.95, 0.98))
cat(sprintf(paste0("Panel C: top 500 of %d autosomal cg probes complete in all 16 (donor-centred variance); ",
                   "mean delta at 200 nM over the 500: AZA %.3f, GSK %.3f; |delta| 95th/98th pct %.2f/%.2f\n"),
            nrow(dl), mean_delta["AZA"], mean_delta["GSK"], q[1], q[2]))

## ---- Panels D/E: GSK vs AZA dose response (13_gsk_vs_aza.R) ----
gva <- GSK_DIR
p1 <- fread(file.path(gva, "part1_dmr507_slopes.csv"))
fwrite(p1, file.path(out_dir, "gsk_vs_aza_dmr507_slopes.csv"))
ps <- readRDS(file.path(gva, "part2_probe_stats.rds"))
vol <- data.table(probe = rownames(ps$tt),
                  slope_diff_beta_0to200 = 2 * (ps$slope_GSK_beta_per100[rownames(ps$tt)] - ps$slope_AZA_beta_per100[rownames(ps$tt)]),
                  t = ps$tt$t, P.Value = ps$tt$P.Value, adj.P.Val = ps$tt$adj.P.Val)
fwrite(vol, file.path(out_dir, "gsk_vs_aza_probe_volcano.csv"))
filt <- fread(file.path(gva, "part2_dmrs_filtered.csv"))
filt_gr <- GRanges(filt$seqnames, IRanges(filt$start, filt$end))
counts <- data.table(
  n_probes = nrow(vol), n_fdr05 = sum(vol$adj.P.Val < 0.05),
  n_fdr05_gsk_higher = sum(vol$adj.P.Val < 0.05 & ps$tt$logFC > 0),
  n_fdr05_gsk_lower = sum(vol$adj.P.Val < 0.05 & ps$tt$logFC < 0),
  n_dmrs_filtered = nrow(filt), n_dmrs_gsk_higher = sum(filt$meandiff > 0), n_dmrs_gsk_lower = sum(filt$meandiff < 0),
  n_dmrs_overlap_507 = sum(overlapsAny(filt_gr, dmrs_filt)))
fwrite(counts, file.path(out_dir, "gsk_vs_aza_counts.csv"))
all_dmrs <- fread(file.path(gva, "part2_dmrs_all.csv"))
all_dmrs[, passes_filter := abs(meandiff) > 0.10]
all_dmrs[, overlaps_507 := overlapsAny(GRanges(seqnames, IRanges(start, end)), dmrs_filt)]
fwrite(all_dmrs, file.path(out_dir, "gsk_vs_aza_dmr_volcano.csv"))
sl <- data.table(probe = rownames(ps$tt), slope_AZA = ps$slope_AZA_beta_per100[rownames(ps$tt)],
                 slope_GSK = ps$slope_GSK_beta_per100[rownames(ps$tt)], adj.P.Val = ps$tt$adj.P.Val)
fwrite(sl[adj.P.Val < 0.05], file.path(out_dir, "gsk_vs_aza_sig_probe_slopes.csv"))
pgr <- mm285_gr[sl$probe]
ovd <- findOverlaps(GRanges(all_dmrs$seqnames, IRanges(all_dmrs$start, all_dmrs$end)), pgr)
reg <- data.table(i = queryHits(ovd), a = sl$slope_AZA[subjectHits(ovd)], g = sl$slope_GSK[subjectHits(ovd)])[
  , .(slope_AZA = mean(a), slope_GSK = mean(g), n_probes = .N), by = i]
dmr_sl <- cbind(all_dmrs[reg$i, .(seqnames, start, end, meandiff, min_smoothed_fdr, passes_filter, overlaps_507)],
                reg[, .(slope_AZA, slope_GSK, n_probes)])
# check: 2 x (GSK - AZA) mean slope reproduces DMRcate meandiff exactly wherever our probe count equals no.cpgs;
# 10 of 3,165 DMRs contain one more contrast probe than DMRcate counted (not explained by duplicate positions;
# max |difference| 0.028), so those 10 use all probes inside the DMR.
dmr_sl[, no.cpgs := all_dmrs$no.cpgs[reg$i]]
err <- abs(2 * (dmr_sl$slope_GSK - dmr_sl$slope_AZA) - dmr_sl$meandiff)
stopifnot(nrow(dmr_sl) == nrow(all_dmrs), max(err[dmr_sl$n_probes == dmr_sl$no.cpgs]) < 1e-6,
          sum(dmr_sl$n_probes != dmr_sl$no.cpgs) <= 10, max(err) < 0.05)
fwrite(dmr_sl, file.path(out_dir, "gsk_vs_aza_dmr_slopes.csv"))
cat(sprintf("Panels G/H: %d FDR < 0.05 probes; %d DMRs (slopes reproduce meandiff except %d DMRs off by one probe)\n", sum(sl$adj.P.Val < 0.05), nrow(dmr_sl), sum(dmr_sl$n_probes != dmr_sl$no.cpgs)))
cat(sprintf("Panel F: %d DMRs (%d pass |meandiff| > 0.10; %d of all overlap the 507)\n", nrow(all_dmrs),
            sum(all_dmrs$passes_filter), sum(all_dmrs$overlaps_507)))
cat(sprintf("Panels D/E: %d DMR slopes; %d probes (%d FDR < 0.05); %d DMRs (%d GSK higher), %d overlap the 507\n",
            nrow(p1), counts$n_probes, counts$n_fdr05, counts$n_dmrs_filtered, counts$n_dmrs_gsk_higher,
            counts$n_dmrs_overlap_507))

cat("=== Done ===\n")
