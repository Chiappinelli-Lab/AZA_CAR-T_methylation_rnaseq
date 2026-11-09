# Differential methylation region calling (DMRcate) for the primary Aza
# dose-response comparison, B7H3 CAR arm only. 
#   - Sample subset: 12 samples 
#   - Design: ~donor + AZA. Donor is a fixed-effect blocking factor (3
#     levels). AZA is coded as a continuous numeric dose (linear
#     dose-response trend, not a 4-level factor) -- coefficient of interest
#     is "AZA".
#   - Per-CpG stats: limma lmFit()+eBayes() on logit M-values, NOT sesame's
#     DML() -- DML() is plain per-CpG OLS via lm() with no eBayes variance
#     shrinkage across probes, which undercuts the small-N rationale (n=3
#     donors) for preferring DMRcate in the first place. M-values (not
#     betas) are used for the linear model since they better satisfy
#     lm/limma's homoscedasticity assumption; a parallel beta-space slope is
#     fit separately, for effect-size reporting only.
#   - Individual-CpG significance: BH FDR < 0.05 (seeds the kernel-smoothed
#     region calling; final DMR calls are gated by dmrcate()'s own
#     pcutoff="fdr" step, not this threshold).
#   - cpg.annotate() does not support MM285 (arraytype is restricted to
#     EPICv2/EPICv1/EPIC/450K), so the CpGannotated object is built manually
#     from sesameData_getManifestGRanges("MM285") + the limma stats above,
#     following the documented dmrcate() example for a manual object.
#   - dmrcate(): lambda=1000, C=2 (package default for array data; MM285's
#     median inter-probe spacing is ~5.8kb, so this only aggregates the
#     ~25% of probe pairs already <1kb apart -- mostly real CGI/promoter
#     clusters -- and leaves isolated probes as single-CpG hits rather than
#     manufacturing regions across unmeasured DNA), consec=FALSE (nucleotide
#     distance, not CpG-count), pcutoff="fdr" (default, package's own
#     recommended primary defense against Type I error), min.cpgs=2
#     (package minimum; MM285's sparse spacing already limits how often a
#     region can span more CpGs than that). betacutoff is left unused --
#     DMRcate's docs restrict it to results from DSS::DMLtest(), which
#     doesn't apply here; effect-size filtering on called regions can be
#     done post-hoc directly from meandiff/maxdiff or the beta matrix.
#   - extractRanges(genome="mm10") needs internet access for its gene-overlap
#     annotation step, which this node doesn't have -- overlapping.genes
#     will come back NA. The core DMR statistics don't depend on that call.
#     Gene annotation is deferred to its own later pipeline stage.


suppressMessages({
  library(DMRcate)
  library(limma)
  library(sesame)
  library(sesameData)
  library(GenomicRanges)
  library(dplyr)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
out_dir <- METH_DM_DIR
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_sheet <- read.csv(METH_SHEET, stringsAsFactors = FALSE)
betas <- readRDS(BETAS_RDS)

# --- subset to the primary comparison group: Aza, B7H3 CAR arm ---
meta <- sample_sheet %>%
  filter(Drug == "5-Azacytidine", CAR == "B7H3") %>%
  transmute(Code, Donor = factor(Donor), AZA = as.numeric(sub("nM", "", Dose)))
stopifnot(nrow(meta) == 12)

betas_sub <- betas[, meta$Code]
cat("Subset:", ncol(betas_sub), "samples x", nrow(betas_sub), "probes\n")
cat("\nDonor x dose design:\n"); print(table(meta$Donor, meta$AZA))

# --- design matrix & limma eBayes moderated t-stats on M-values ---
design <- model.matrix(~ Donor + AZA, data = meta)
cat("\nDesign matrix:\n"); print(design)

mvals <- BetaValueToMValue(betas_sub)
keep <- rowSums(is.na(mvals)) == 0
cat("\nDropping", sum(!keep), "of", nrow(mvals), "probes with any missing value in this subset\n")
mvals <- mvals[keep, ]
betas_keep <- betas_sub[keep, ]

fit <- lmFit(mvals, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "AZA", number = Inf, sort.by = "none")

# beta-space slope (per nM), for effect-size reporting only -- significance
# comes entirely from the M-value model above.
fit_beta <- lmFit(betas_keep, design)
beta_slope <- coef(fit_beta)[, "AZA"]

cat("\nN CpGs with FDR < 0.05 (AZA coefficient):", sum(tt$adj.P.Val < 0.05), "of", nrow(tt), "\n")

# --- manual CpGannotated object (cpg.annotate() doesn't support MM285) ---
source(file.path(WORKFLOW, "03_methylation", "mm39_annotation.R"))  # mm39 / GENCODE vM32 annotation (since 2026-10-07)
mm285_gr <- mm285_gr_mm39()
probe_ids <- rownames(tt)
has_coord <- probe_ids %in% names(mm285_gr)
cat("Dropping", sum(!has_coord), "of", length(probe_ids),
    "probes with no MM285 genomic coordinate (control/masked probes)\n")

probe_ids <- probe_ids[has_coord]
tt <- tt[probe_ids, ]
gr_coord <- mm285_gr[probe_ids]

annotated <- GRanges(
  seqnames = seqnames(gr_coord),
  ranges = ranges(gr_coord),
  stat = tt$t,
  diff = beta_slope[probe_ids] * 200,  # fitted beta-value shift over the full 0-200nM dose range
  rawpval = tt$P.Value,
  ind.fdr = tt$adj.P.Val,
  is.sig = tt$adj.P.Val < 0.05
)
names(annotated) <- probe_ids
myannotation <- new("CpGannotated", ranges = annotated)

# --- DMRcate region calling (parameters confirmed with user 2026-09-09) ---
cat("\nRunning dmrcate() ...\n")
dmrcoutput <- dmrcate(myannotation, lambda = 1000, C = 2, consec = FALSE,
                       pcutoff = "fdr", min.cpgs = 2)

dmrs <- extractRanges(dmrcoutput, genome = "mm10")
# genome="mm10" only drives the overlapping.genes annotation (coordinates are
# mm39); replace it with GENCODE vM32 gene overlaps
g39 <- vm32_genes_gr()
ov_g <- findOverlaps(dmrs, g39, ignore.strand = TRUE)
og <- tapply(g39$gene_name[subjectHits(ov_g)], queryHits(ov_g), function(x) paste(unique(x), collapse = ", "))
mcols(dmrs)$overlapping.genes <- NA_character_
mcols(dmrs)$overlapping.genes[as.integer(names(og))] <- og
cat("\nCalled", length(dmrs), "DMRs\n")

dmrs_df <- as.data.frame(dmrs)
write.csv(dmrs_df, file.path(out_dir, "dmrcate_dmrs_aza_b7h3.csv"), row.names = FALSE)

saveRDS(list(fit = fit, tt = tt, myannotation = myannotation,
             dmrcoutput = dmrcoutput, dmrs = dmrs),
        file.path(out_dir, "dmrcate_aza_b7h3.rds"))

cat("\nSaved DMR table to", file.path(out_dir, "dmrcate_dmrs_aza_b7h3.csv"), "\n")
cat("Saved full result object to", file.path(out_dir, "dmrcate_aza_b7h3.rds"), "\n")

if (length(dmrs) > 0) {
  cat("\n=== Top DMRs (by min_smoothed_fdr) ===\n")
  top <- as.data.frame(dmrs)[order(dmrs$min_smoothed_fdr),
                              c("seqnames", "start", "end", "no.cpgs",
                                "min_smoothed_fdr", "Stouffer", "maxdiff", "meandiff")]
  print(head(top, 10), row.names = FALSE)
}

# --- effect-size filter, confirmed with user 2026-09-09 ---
# All dmrcate() output is already FDR-significant (pcutoff="fdr"); on top of
# that, require |meandiff| > 0.10 (region-averaged fitted beta shift across
# the full 0-200nM dose range), analogous to the >10% methylation-difference
# filter used in the WGBS pipeline. meandiff (region average), not maxdiff
# (single peak CpG), was chosen so one large-effect CpG can't carry an
# otherwise weak region.
dmrs_filt <- dmrs[abs(mcols(dmrs)$meandiff) > 0.10]
cat("\nApplying |meandiff| > 0.10 filter:", length(dmrs_filt), "of", length(dmrs), "DMRs retained\n")
write.csv(as.data.frame(dmrs_filt),
          file.path(out_dir, "dmrcate_dmrs_aza_b7h3_filtered.csv"), row.names = FALSE)
cat("Saved filtered DMR table to",
    file.path(out_dir, "dmrcate_dmrs_aza_b7h3_filtered.csv"), "\n")
