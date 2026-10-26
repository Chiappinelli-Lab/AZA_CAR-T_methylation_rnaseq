# Software, versions and reference data

Versions are the ones installed in the environments that produced the results (exported 2026-10-20).
Full pins: `envs/<env>.yaml` (`conda env export --no-builds`), `envs/<env>.explicit.txt` (`conda list --explicit`,
exact builds), `envs/snakemake_rule_envs/` (the per-rule Snakemake envs), `envs/<env>.R_installed_packages.tsv`
(every R package in the env), and `envs/R_packages_non_conda.txt` (R packages installed outside conda, so not in
the conda exports).

## Conda environments

| Env | Used by | Main contents |
|---|---|---|
| `snakemake.RNAseq.standard.with.SQuiRE` | stage 1 controller | Snakemake 8.20.6, snakemake-executor-plugin-slurm 0.12.1, Python 3.11 |
| per-rule envs (`envs/snakemake_rule_envs/`, created by Snakemake from `workflow/01_rnaseq_pipeline/envs/`) | stage 1 rules | see RNA-seq table |
| `deseq2` | stage 2; 12_dmgsea_plots, integration 01/03/10; fig2_01, fig3_01 | R 4.5.3, Bioconductor 3.22 |
| `haydar_methylation` | stage 3; integration 02/11/12; fig2_02, fig3_02, fig3_04, figS_gsk_02 | R 4.5.3, Bioconductor 3.22 |
| `data_analysis` | figure plotting (`*_plot.py`), `make_fastq_links.py`, config reads in `run_all.sh` | Python 3.11.9 |

## RNA-seq processing (stage 1)

| Step | Tool | Version | Where it runs |
|---|---|---|---|
| Raw-read QC | FastQC | 0.12.1 | rule env `fastqc` |
| Adapter/quality trimming (+ FastQC on trimmed reads) | Trim Galore | 0.6.10 (cutadapt 5.2) | rule env `trimGalore` |
| Alignment for QC metrics | STAR / samtools sort | 2.7.0e / 1.21 | rule env `star` |
| RNA-seq metrics | Picard CollectRnaSeqMetrics | 3.2.0 | rule env `picard` |
| QC report | MultiQC | 1.25.2 | rule env `multiqc` |
| **Gene counts** (Build STAR index, Clean, Map, Count) | **SQuIRE** | **0.9.9.92** (Apptainer image `squire_0.9.9.92.sif`; inside: STAR 2.5.3a, samtools 1.1, bedtools 2.25.0, StringTie 1.3.3, Python 2.7.14) | Apptainer 1.3.0 |
| `join` bound into the SQuIRE image (missing from it) | GNU coreutils | 9.5 | host binary |
| Count matrix merge (`combine_squire_counts.py`) | Python / pandas | 3.11 / 2.2.2 | rule env `combine_counts` |

The gene counts used downstream (`squire_gene_counts.tsv`) come from SQuIRE's own STAR 2.5.3a alignment of the
trimmed reads. The STAR 2.7.0e alignment feeds only Picard/MultiQC QC.

Parameters (`workflow/01_rnaseq_pipeline/Snakefile`, `config.yaml`): paired-end 151 bp reads, stranded
(dUTP-type, reverse); Trim Galore `--paired --gzip --fastqc`; QC STAR `--outFilterMultimapNmax 100
--winAnchorMultimapNmax 200`; SQuIRE build `--sjdbOverhang 150`, `-r 151`, Count strandedness `-s 1`, build label mm39.

## RNA-seq differential expression and GSEA (stage 2, env `deseq2`)

| Package | Version | Use |
|---|---|---|
| DESeq2 | 1.50.2 | `~ donor + dose_100nM` (dose continuous, nM/100); NS and S1 fits |
| apeglm | 1.32.0 | `lfcShrink(type = "apeglm")` |
| clusterProfiler | 4.18.4 | `gseGO` (BP, 15–500 genes, BH, `eps = 0`, `seed = TRUE`), `simplify` (0.7, Wang) |
| fgsea | 1.36.2 | GSEA engine under clusterProfiler |
| enrichplot / DOSE / GOSemSim | 1.30.4 / 4.4.0 / 2.36.0 | `setReadable`, simplify |
| org.Mm.eg.db | 3.22.0 | Ensembl → symbol/Entrez, GO annotation |
| GO.db | 3.22.0 | GO term data |
| limma | 3.66.0 | `removeBatchEffect` for the donor-corrected Fig 2 heatmap (fig2_01); installed outside conda |
| data.table / ggplot2 / pheatmap | 1.17.8 / 4.0.3 / 1.0.13 | |

Gene filter: ≥10 counts in ≥3 samples.

## Methylation array (stage 3, env `haydar_methylation`)

| Package | Version | Use |
|---|---|---|
| sesame | 1.28.1 | `readIDATpair`, `prepSesame(prep = "QCDPBM")` with custom quality mask (`ref_issue`, `design_issue`), pOOBAH p < 0.05, noob, `matchDesign`; `getBetas` |
| sesameData | 1.28.0 | MM285 mm10 manifest (`sesameData_getManifestGRanges("MM285")`), masks; ExperimentHub snapshot 2026-10-07 (local cache) |
| ExperimentHub / AnnotationHub | 3.0.0 / 4.0.0 | |
| limma | 3.66.0 | `lmFit`/`eBayes` on M-values, `~ Donor + AZA` |
| DMRcate | 3.6.0 | `dmrcate(lambda = 1000, C = 2, consec = FALSE, pcutoff = "fdr", min.cpgs = 2)`, then \|meandiff\| > 0.10 |
| dmGsea | 1.0.2 | `gsProbe` (goregion-style), GO:BP via `getGO`, 15–500 genes; `p.gs` patched to P(X ≥ k) (see `11_dmgsea.R`) |
| BiasedUrn | 2.0.12 | Fisher noncentral hypergeometric (dmGsea) |
| rtracklayer | 1.70.1 | `liftOver` with the UCSC chain, GTF import |
| txdbmaker | 1.6.2 | TxDb from the vM32 GTF (installed outside conda) |
| GenomicRanges / GenomeInfoDb | 1.62.1 / 1.46.2 | |
| org.Mm.eg.db / GO.db | 3.22.0 / 3.22.0 | Ensembl → Entrez/symbol; GO (org.Mm.eg.db installed outside conda in this env) |
| GOSemSim | 2.36.0 | Wang semantic similarity of the shared terms (fig2_02) |
| ComplexHeatmap / circlize | 2.26.1 / 0.4.18 | candidate-gene heatmaps (integration 12) |
| readxl / dplyr / data.table | 1.5.0 / 1.2.1 / 1.18.6.1 | |

## Figures (env `data_analysis`)

Python 3.11.9, numpy 2.4.6, pandas 3.0.3, matplotlib 3.11.0, matplotlib-venn 1.1.2, scipy 1.17.1, pillow 12.2.0.
The paper figures use Arial 2.82 (Microsoft core web fonts, conda-forge `mscorefonts`), which is not shipped here
because of its license. Without it the scripts use Liberation Sans (shipped in `workflow/05_figures/fonts/`;
metric-compatible, so the layout is the same and only the glyph shapes differ). To get Arial, put `arial*.ttf` in
that folder.

## Reference data and annotation

| Item | Version / source | Used for |
|---|---|---|
| Genome | GRCm39 primary assembly (GENCODE), `GRCm39.primary_assembly.genome.fa` | SQuIRE STAR index |
| Gene annotation | GENCODE vM32 primary_assembly GTF | STAR, SQuIRE, methylation gene models + TxDb |
| QC STAR index | STAR 2.7.0e index on GRCm39 + vM32, `--sjdbOverhang 150` | QC alignment |
| refFlat / rRNA interval list | derived from vM32 | Picard |
| RepeatMasker | UCSC goldenPath mm39 `mm39.fa.out` (3-line header removed) | SQuIRE Clean |
| MM285 probe coordinates | sesameData MM285 manifest (mm10), lifted to mm39 with the UCSC `mm10ToMm39.over.chain` (1:1, same width, standard chromosomes) | all methylation coordinates |
| Native mm39 MM285 manifest | Zhou lab InfiniumAnnotationV1 `MM285.mm39.manifest.tsv.gz` | liftover cross-check only (06b) |
| ENCODE cCREs | UCSC `encodeCcreCombined` (mm10), lifted to mm39 | built by 00b, not used by the figures |
| Gene sets | GO:BP from org.Mm.eg.db / GO.db 3.22.0 | RNA GSEA, dmGsea |
