# Publication-repo split of 00b_build_mm39_annotation.R: the two parts of 00b that
# read the 06_dmrcate.R result, moved here so the stages run in order
# 00b -> 06 -> 06b. Code is 00b sections 2 and 6, unchanged apart from paths.
#   1. Cross-check of the lifted MM285 manifest against the Zhou lab's native
#      mm39 manifest (log + liftover_vs_native_mm39*.csv). Stops if more than
#      0.1% of uniquely-mapping probes disagree.
#   2. results/methylation/dm/tested_gene_background.rds (genes_gr is read by
#      11_dmgsea.R).
# Run in env haydar_methylation, after 06_dmrcate.R.

suppressMessages({
  library(GenomicRanges)
  library(data.table)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
source(file.path(WORKFLOW, "03_methylation", "mm39_annotation.R"))
ref_dir  <- MM39_REF_DIR
m10 <- sesameData::sesameData_getManifestGRanges("MM285")
m39 <- mm285_gr_mm39()
genes_gr <- vm32_genes_gr()
dm <- readRDS(DMRCATE_RDS)

## ---- 1. cross-check against the native mm39 manifest (00b section 2) ----
nat <- fread(MM285_NATIVE_MM39,
             select = c("CpG_chrm", "CpG_beg", "CpG_end", "Probe_ID", "mapQ_A", "mapQ_B"))
nat <- nat[!is.na(CpG_beg)]
# probe mapping quality in the native alignment (Infinium I: worse of the A/B probes)
nat[, mapQ := fifelse(is.na(mapQ_B), as.numeric(mapQ_A), pmin(as.numeric(mapQ_A), as.numeric(mapQ_B)))]
tested <- rownames(dm$tt)
cmp <- data.table(Probe_ID = names(m39), lift_chr = as.character(seqnames(m39)), lift_start = start(m39))
cmp <- merge(cmp, nat[, .(Probe_ID, nat_chr = CpG_chrm, nat_start = CpG_beg + 1L, mapQ)], by = "Probe_ID", all.x = TRUE)
cmp[, status := fifelse(is.na(nat_chr), "not_in_native",   # ch/rs probes: native manifest has no CpG coordinate
                fifelse(lift_chr != nat_chr, "different_chr",
                fifelse(lift_start == nat_start, "identical",
                fifelse(abs(lift_start - nat_start) <= 2, "within_2bp", "different_position"))))]
cmp[, `:=`(unique_mapper = !is.na(mapQ) & mapQ >= 10, tested = Probe_ID %in% tested)]
agree <- cmp[, .(N = .N, n_unique_mapper = sum(unique_mapper), n_tested = sum(tested)), by = status
             ][, pct := round(100 * N / sum(N), 3)][order(-N)]
cat("\nLifted vs native mm39 manifest (lifted probes):\n"); print(agree)
fwrite(agree, file.path(ref_dir, "liftover_vs_native_mm39.csv"))
fwrite(cmp[status %in% c("different_chr", "different_position")],
       file.path(ref_dir, "liftover_vs_native_mm39_disagreements.csv"))
n_native_only <- nat[!Probe_ID %in% names(m39) & Probe_ID %in% names(m10), .N]
cat("Probes with a native mm39 position but dropped by liftover:", n_native_only, "\n")
# The first run (2026-10-07) found 1.4% of lifted probes disagreeing with the native
# manifest -- every one a multi-mapping probe (native mapQ < 10, i.e. the re-alignment
# picked another copy of a repeat), and all but 23 already removed by sesame's quality
# mask before testing. Liftover keeps such probes where sesame's mm10 manifest placed
# them, which is the placement every past result used. The stop rule therefore tests
# only uniquely-mapping probes, where a disagreement would mean the lift itself is wrong.
bad_unique <- cmp[unique_mapper & status %in% c("different_chr", "different_position"), .N] / cmp[unique_mapper == TRUE, .N]
cat(sprintf("Disagreement among uniquely-mapping probes (native mapQ >= 10): %.4f%%\n", 100 * bad_unique))
stopifnot(bad_unique < 0.001)

## ---- 2. tested_gene_background.rds (00b section 6; same fields as the old mm10 file) ----
tested_ids <- intersect(rownames(dm$tt), names(m39))
anno <- m39[tested_ids]
nn <- nearest(anno, genes_gr, ignore.strand = TRUE)
anno$gene_name <- genes_gr$gene_name[nn]
ov <- findOverlaps(anno, genes_gr, maxgap = 5000, ignore.strand = TRUE)
tested_genes <- unique(genes_gr$gene_name[subjectHits(ov)])
saveRDS(list(anno = anno, tested_genes = tested_genes, genes_gr = genes_gr),
        file.path(METH_DM_DIR, "tested_gene_background.rds"))
cat(sprintf("\ntested_gene_background.rds: %d tested probes on mm39, %d genes within 5kb\n",
            length(anno), length(tested_genes)))
