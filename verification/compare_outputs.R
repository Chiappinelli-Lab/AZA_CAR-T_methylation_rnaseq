#!/usr/bin/env Rscript
# Compare every non-figure output in verification/manifest.tsv (current project = reference, publication repo
# rerun = new). Writes verification/output_comparison.tsv. Run in env deseq2; the DMRcate rds needs env
# haydar_methylation (run with a path regex for that row).
#   csv/tsv  byte-identical -> "exact"; otherwise read both, align rows (same order, else sorted by the
#            character columns), and report per column: numeric max |diff| and max relative diff, character
#            mismatches. Status "numeric" if the only differences are numeric with relative diff <= 1e-8.
#            sample_sheet.csv Basename (absolute IDAT path) is compared by basename().
#   rds      DESeqDataSet: counts, sizeFactors, dispersions, results(dose_100nM); lists: element by element;
#            other objects: all.equal. Reports the largest all.equal discrepancy.
#   other    text count/BED files by md5; FASTQ/BAM/STAR index/sqlite are not compared (large binaries whose
#            content is checked through the count matrices they produce).
# usage: Rscript compare_outputs.R [manifest.tsv] [path-regex to restrict rows]
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  if (requireNamespace("DESeq2", quietly = TRUE)) library(DESeq2)
  if (requireNamespace("DMRcate", quietly = TRUE)) library(DMRcate)
})
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
repo <- dirname(here)
man <- fread(if (length(args) >= 1) args[1] else file.path(here, "manifest.tsv"))
man <- man[!type %in% c("pdf", "svg", "png")]
if (length(args) >= 2) man <- man[grepl(args[2], new_path)]
REL_TOL <- 1e-8

rel <- function(x, y) { d <- abs(x - y); m <- pmax(abs(x), abs(y)); ifelse(d == 0, 0, d / m) }

cmp_table <- function(a, b, path) {
  A <- fread(a, sep = if (grepl("\\.tsv$", a)) "\t" else "auto")
  B <- fread(b, sep = if (grepl("\\.tsv$", b)) "\t" else "auto")
  if (!identical(names(A), names(B))) return(list("FAIL", sprintf("columns differ: %s", paste(setdiff(union(names(A), names(B)), intersect(names(A), names(B))), collapse = ","))))
  if (nrow(A) != nrow(B)) return(list("FAIL", sprintf("rows %d vs %d", nrow(A), nrow(B))))
  if (basename(path) == "sample_sheet.csv") { A[, Basename := basename(Basename)]; B[, Basename := basename(Basename)] }
  chr <- names(A)[!vapply(A, is.numeric, TRUE)]
  if (length(chr) && !identical(A[, ..chr], B[, ..chr])) {   # same rows in another order?
    setorderv(A, chr); setorderv(B, chr)
  }
  notes <- character(); worst <- 0; chr_bad <- FALSE
  for (cn in names(A)) {
    x <- A[[cn]]; y <- B[[cn]]
    if (is.numeric(x) && is.numeric(y)) {
      if (!identical(is.na(x), is.na(y))) { notes <- c(notes, sprintf("%s: NA pattern differs", cn)); chr_bad <- TRUE; next }
      ok <- !is.na(x)
      if (any(x[ok] != y[ok])) {
        r <- max(rel(x[ok], y[ok])); worst <- max(worst, r)
        notes <- c(notes, sprintf("%s: max|d| %.3g, max rel %.3g, %d rows", cn, max(abs(x[ok] - y[ok])), r, sum(x[ok] != y[ok])))
      }
    } else if (!identical(as.character(x), as.character(y))) {
      chr_bad <- TRUE
      notes <- c(notes, sprintf("%s: %d values differ", cn, sum(as.character(x) != as.character(y), na.rm = TRUE)))
    }
  }
  if (!length(notes)) return(list("exact", "identical content (formatting differs)"))
  list(if (!chr_bad && worst <= REL_TOL) "numeric" else "DIFF", paste(notes, collapse = "; "))
}

ae <- function(x, y) { r <- all.equal(x, y, tolerance = 0); if (isTRUE(r)) "identical" else paste(head(r, 3), collapse = " | ") }
ae_tol <- function(x, y) isTRUE(all.equal(x, y, tolerance = REL_TOL))

cmp_rds <- function(a, b) {
  A <- readRDS(a); B <- readRDS(b)
  if (is(A, "DESeqDataSet")) {  # DESeq2 attached when available (env deseq2)
    parts <- list(counts = list(counts(A), counts(B)), sizeFactors = list(sizeFactors(A), sizeFactors(B)),
                  dispersions = list(dispersions(A), dispersions(B)),
                  results = list(as.data.frame(results(A, name = "dose_100nM")), as.data.frame(results(B, name = "dose_100nM"))))
  } else if (is.list(A) && !is.data.frame(A) && !is(A, "GRanges")) {
    if (!identical(names(A), names(B))) return(list("FAIL", "list elements differ"))
    parts <- setNames(lapply(names(A), function(n) list(A[[n]], B[[n]])), names(A))
  } else parts <- list(object = list(A, B))
  msgs <- vapply(names(parts), function(n) {
    x <- parts[[n]][[1]]; y <- parts[[n]][[2]]
    if (is(x, "MArrayLM")) { x <- unclass(x)[c("coefficients", "stdev.unscaled", "sigma", "df.residual")]; y <- unclass(y)[names(x)] }
    if (is(x, "CpGannotated")) { x <- x@ranges; y <- y@ranges }
    if (is(x, "DMRcate")) { x <- lapply(slotNames(x), function(n) slot(x, n)); y <- lapply(slotNames(y), function(n) slot(y, n)) }
    m <- tryCatch(ae(x, y), error = function(e) paste("compare error:", conditionMessage(e)))
    if (m == "identical") "identical" else sprintf("%s [%s]", m, if (tryCatch(ae_tol(x, y), error = function(e) FALSE)) "within 1e-8" else "beyond 1e-8")
  }, "")
  st <- if (all(msgs == "identical")) "exact" else if (all(msgs == "identical" | grepl("within 1e-8\\]$", msgs))) "numeric" else "DIFF"
  list(st, paste(sprintf("%s: %s", names(msgs), msgs), collapse = "; "))
}

out <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  r <- man[i]; a <- r$current_path; b <- file.path(repo, r$new_path)
  res <- if (!file.exists(b) && !nzchar(Sys.readlink(b))) list("MISSING", "not produced by the rerun") else
    tryCatch({
      if (r$type %in% c("csv", "tsv")) {
        if (tools::md5sum(a) == tools::md5sum(b)) list("exact", "byte-identical") else cmp_table(a, b, r$new_path)
      } else if (r$type == "rds") cmp_rds(a, b)
      else if (grepl("\\.(fq\\.gz|fastq\\.gz|bam|sqlite)$|/Genome$", b)) list("skipped", "large binary; checked through downstream counts")
      else if (nzchar(Sys.readlink(b))) list(if (basename(Sys.readlink(a)) == basename(Sys.readlink(b))) "exact" else "DIFF", paste("symlink ->", Sys.readlink(b)))
      else if (tools::md5sum(a) == tools::md5sum(b)) list("exact", "byte-identical") else list("DIFF", "md5 differs")
    }, error = function(e) list("ERROR", conditionMessage(e)))
  message(sprintf("[%s] %s", res[[1]], r$new_path))
  data.table(new_path = r$new_path, type = r$type, status = res[[1]], detail = res[[2]])
}))
fwrite(out, file.path(here, if (length(args) >= 2) "output_comparison_subset.tsv" else "output_comparison.tsv"), sep = "\t")
print(out[, .N, by = status])
