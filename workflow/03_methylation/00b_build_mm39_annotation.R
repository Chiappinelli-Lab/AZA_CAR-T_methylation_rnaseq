# Build the mm39 / GENCODE vM32 annotation the methylation side uses since
# 2026-10-07, so methylation coordinates and probe -> gene assignment share the
# RNA-seq's genome and gene annotation (RNA-seq: mm39, GENCODE vM32). 
#   - Probe coordinates: sesame's mm10 MM285 manifest lifted with the UCSC
#     mm10ToMm39 chain. A probe is kept only if it lifts to exactly one range of
#     the same width on a standard chromosome (chr1-19/X/Y/M); everything else is
#     written to mm285_liftover_dropped.csv with the reason.
#   - Cross-check only (not used downstream): the Zhou lab's native mm39 MM285
#     manifest (InfiniumAnnotationV1, probes re-aligned to mm39). Agreement is
#     reported in liftover_vs_native_mm39.csv.
#   - Gene models: "gene" records of the exact GTF the RNA-seq was quantified
#     with (gencode.vM32.primary_assembly.annotation.gtf), all gene types,
#     standard chromosomes; names = versioned Ensembl IDs (match RNA-seq
#     `feature`), mcols gene_name / gene_type -- same structure as the old
#     sesameData genes_gr, so downstream code is unchanged.
#   - TxDb for exon/intron calls (06d/06e) built from the same GTF.
#   - cCREs: the existing ENCODE mm10 combined set lifted to mm39 (no mm39 cCRE
#     registry on disk); same 1:1 / same-width rule as probes.
#   - results/dm/tested_gene_background.rds is rebuilt with the same three
#     fields as before: genes_gr (above); anno = tested probes (the 06_dmrcate
#     probe set that lifted) with gene_name of the nearest vM32 gene;
#     tested_genes = genes with a tested probe within 5kb (07's window rule).
#     Only genes_gr is read by current scripts.
# Run in env haydar_methylation, before 06_dmrcate.R.
# The diagnostic count of RNA-seq genes missing from genes_gr is dropped.

suppressMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(sesameData)
  library(data.table)
  library(txdbmaker)
  library(AnnotationDbi)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
ref_dir  <- MM39_REF_DIR
chain_f  <- CHAIN_MM10_MM39
gtf_f    <- GTF_VM32
dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)

std_chrs <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
chain <- import.chain(chain_f)

# 1:1, same-width liftover; returns list(gr = lifted GRanges, dropped = data.table)
lift_1to1 <- function(gr, id) {
  lo <- liftOver(gr, chain)
  n  <- lengths(lo)
  one <- n == 1L
  lifted <- unlist(lo[one])
  ok_w   <- width(lifted) == width(gr[one])
  ok_chr <- as.character(seqnames(lifted)) %in% std_chrs
  keep <- which(one)[ok_w & ok_chr]
  out <- unlist(lo[keep])
  names(out) <- id[keep]
  reason <- rep(NA_character_, length(gr))
  reason[n == 0L] <- "unmapped"
  reason[n > 1L]  <- "split"
  reason[which(one)[!ok_w]] <- "width_changed"
  reason[which(one)[ok_w & !ok_chr]] <- "nonstandard_chr"
  dropped <- data.table(id = id, mm10_chr = as.character(seqnames(gr)), mm10_start = start(gr),
                        reason = reason)[!is.na(reason)]
  list(gr = out, dropped = dropped)
}

## ---- 1. MM285 manifest -> mm39 ----
m10 <- sesameData_getManifestGRanges("MM285")
lp <- lift_1to1(m10, names(m10))
m39 <- lp$gr
m39 <- keepSeqlevels(m39, intersect(std_chrs, seqlevels(m39)), pruning.mode = "coarse")
cat(sprintf("MM285 manifest: %d probes with mm10 coords; lifted 1:1 to mm39: %d (%.2f%%); dropped: %d\n",
            length(m10), length(m39), 100 * length(m39) / length(m10), nrow(lp$dropped)))
print(lp$dropped[, .N, by = reason])
stopifnot(length(m39) / length(m10) >= 0.99)
saveRDS(m39, file.path(ref_dir, "mm285_manifest_mm39.rds"))
fwrite(lp$dropped, file.path(ref_dir, "mm285_liftover_dropped.csv"))

## ---- 3. vM32 gene models ----
gtf_genes <- import(gtf_f, format = "gtf", feature.type = "gene",
                    colnames = c("gene_id", "gene_name", "gene_type"))
gtf_genes <- gtf_genes[as.character(seqnames(gtf_genes)) %in% std_chrs]
gtf_genes <- keepSeqlevels(gtf_genes, intersect(std_chrs, seqlevels(gtf_genes)), pruning.mode = "coarse")
genes_gr <- granges(gtf_genes)
names(genes_gr) <- gtf_genes$gene_id
mcols(genes_gr) <- DataFrame(gene_name = gtf_genes$gene_name, gene_type = gtf_genes$gene_type)
stopifnot(!anyDuplicated(names(genes_gr)))
cat(sprintf("\nvM32 genes on chr1-19/X/Y/M: %d\n", length(genes_gr)))
saveRDS(genes_gr, file.path(ref_dir, "genes_gr_vM32.rds"))

## ---- 4. TxDb from the same GTF ----
txdb_f <- file.path(ref_dir, "txdb_vM32.sqlite")
# organism/taxonomyId left NA: setting either makes txdbmaker validate it through
# GenomeInfoDbData, which isn't installed in this env; recorded as metadata instead.
txdb <- makeTxDbFromGFF(gtf_f, format = "gtf",
                        dataSource = "GENCODE vM32 primary_assembly (RNA-seq GTF)",
                        metadata = data.frame(name = c("Genome", "Species"), value = c("mm39", "Mus musculus")))
saveDb(txdb, txdb_f)
cat("TxDb saved:", txdb_f, "\n")

## ---- 5. ENCODE cCREs -> mm39 ----
cc <- fread(CCRE_MM10_BED, header = FALSE, sep = "\t")
cc_gr <- GRanges(cc$V1, IRanges(cc$V2 + 1L, cc$V3))
lc <- lift_1to1(cc_gr, cc$V4)
cc39 <- cc[match(names(lc$gr), cc$V4)]
cc39[, `:=`(V1 = as.character(seqnames(lc$gr)), V2 = start(lc$gr) - 1L, V3 = end(lc$gr),
            V7 = start(lc$gr) - 1L, V8 = end(lc$gr))]
cat(sprintf("\ncCREs: %d mm10; lifted 1:1 to mm39: %d (%.2f%%)\n", nrow(cc), nrow(cc39), 100 * nrow(cc39) / nrow(cc)))
print(lc$dropped[, .N, by = reason])
fwrite(cc39, file.path(ref_dir, "encodeCcreCombined_mm39.bed"), sep = "\t", col.names = FALSE, quote = FALSE)
fwrite(lc$dropped, file.path(ref_dir, "encodeCcre_liftover_dropped.csv"))

