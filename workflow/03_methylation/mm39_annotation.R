# Shared mm39 / GENCODE vM32 annotation accessors for the methylation side and
# every downstream script (built by 00b_build_mm39_annotation.R).
# Replaces sesameData_getManifestGRanges("MM285") (mm10), sesameData mm10 gene
# models, TxDb.Mmusculus.UCSC.mm10.knownGene and the mm10 ENCODE cCRE BED.

# MM39_REF_DIR is set in R/paths.R (results/methylation/reference/mm39)
if (!exists("MM39_REF_DIR")) source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))

# MM285 probe coordinates on mm39 (liftover of sesame's mm10 manifest, 1:1 probes only)
mm285_gr_mm39 <- function() readRDS(file.path(MM39_REF_DIR, "mm285_manifest_mm39.rds"))

# GENCODE vM32 gene models (the RNA-seq GTF): names = versioned Ensembl ID, mcols gene_name, gene_type
vm32_genes_gr <- function() readRDS(file.path(MM39_REF_DIR, "genes_gr_vM32.rds"))

# TxDb built from the same vM32 GTF
vm32_txdb <- function() AnnotationDbi::loadDb(file.path(MM39_REF_DIR, "txdb_vM32.sqlite"))

# ENCODE combined cCREs, lifted mm10 -> mm39 (same columns as reference/encodeCcreCombined.bed)
ccre_bed_path <- function() file.path(MM39_REF_DIR, "encodeCcreCombined_mm39.bed")
