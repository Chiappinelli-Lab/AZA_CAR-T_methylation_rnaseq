# Azacitidine dose response in B7H3 CAR-T cells: RNA-seq × DNA methylation

Code that regenerates Figure 2, Figure 3 and the GSK-3484862 supplementary figure from the raw sequencing
reads (FASTQ) and methylation array intensities (IDAT). It covers:
- bulk RNA-seq of mouse B7H3 CAR-T cells treated with 0/50/100/200 nM azacitidine (AZA);
- Infinium Mouse Methylation BeadChip (MM285) profiling of the same treatment;
- their gene- and pathway-level integration.

> Citation, data accessions and license are placeholders until publication: see `CITATION.md` and `LICENSE`.

## Study design, and why the two assays are integrated at the gene level

| | RNA-seq | Methylation array (MM285) |
|---|---|---|
| Donors | EV07, EV17, EV18, EV19, EV20 | EV07, EV08, ZY03 (AZA); EV03 (GSK-3484862) |
| Treatment | AZA 0/50/100/200 nM | AZA or GSK-3484862, 0/50/100/200 nM |
| CAR | B7H3 | B7H3 (primary) and STOP control CAR |
| Conditions | non-stimulated (NS, 18 samples) and stimulated (S1, 12 samples; GL261-B7H3 co-culture) | no stimulation |
| Primary model | DESeq2 `~ donor + dose_100nM` (NS) | limma/DMRcate `~ Donor + AZA`, AZA + B7H3 CAR, 12 samples |

Only donor EV07 is in both cohorts, so the two assays are joined per gene (dose-response effect vs dose-response
effect) or per pathway, never per sample.

## Repository layout

```
config/paths.example.yaml  template for config/paths.yaml: every external input (raw data, references,
                           containers, SLURM settings)
metadata/                  RNA-seq sample manifest (.tsv) and methylation array manifest (.xlsx); not tracked
R/paths.R, py/paths.py     all in/out paths, derived from REPO_ROOT + config/paths.yaml
run_all.sh                 driver: submits stages 1-5 as chained SLURM jobs
workflow/
  01_rnaseq_pipeline/      Snakefile, config.yaml, rule conda envs, SLURM profile template,
                           make_fastq_links.py, scripts/combine_squire_counts.py
  02_rnaseq_de/            DESeq2 NS + S1 fits, GO:BP GSEA
  03_methylation/          sample sheet, mm39 annotation, sesame preprocessing, DMRcate, dmGsea, GSK vs AZA
  04_integration/          gene-level slopes, multi-omic table, candidate gene sets
  05_figures/              per-figure data prep (R) and plotting (Python), fonts/
  run_stage.sh             runs one of stages 2-5 (script order + conda env per script)
envs/                      conda exports, R package lists (see SOFTWARE.md)
results/  figures/  logs/  outputs (not tracked)
verification/             rerun-vs-original comparison: compare scripts, REPORT.md
```

## Setup

0. **Config and metadata.** `cp config/paths.example.yaml config/paths.yaml` and fill it in (steps 1, 2 and 4).
   The sample metadata is not in this repository. Put the two manifests, with the file names the scripts expect,
   in `metadata/`: `haydar_rnaseq_sample_manifest.tsv` (RNA-seq libraries, donors, doses, FASTQ names) and
   `Methyl Array Sample Manifest EV.xlsx` (array sample codes, donors, CAR, drug, dose, IDAT basenames). The
   sample information is in the GEO records (see `CITATION.md`).
1. **Data.** Put the 60 FASTQs (30 libraries × R1/R2, names as delivered) in one folder. Put the 24 IDAT pairs in
   one folder with a subfolder per Sentrix chip. Set both in `config/paths.yaml` under `data:`. (Accessions: see
   `CITATION.md`.)
2. **References** (`config/paths.yaml` under `references:` and `rnaseq_pipeline:`):
   - GRCm39 primary-assembly FASTA and GENCODE vM32 primary-assembly GTF;
   - a STAR index built from both with `--sjdbOverhang 150`;
   - a vM32 refFlat and rRNA interval list (for Picard);
   - the UCSC mm39 RepeatMasker `.fa.out` with its 3-line header removed;
   - the UCSC `mm10ToMm39.over.chain`;
   - the UCSC `encodeCcreCombined` mm10 BED;
   - the Zhou lab InfiniumAnnotation `MM285.mm39.manifest.tsv.gz`.

   Versions are in `SOFTWARE.md`.
3. **Software.**
   - SQuIRE 0.9.9.92 as an Apptainer image (`SQuIRE_container`), plus a GNU coreutils `join` binary
     (`slurm.squire_join_binary`): the image lacks `join`, and the binary is bound in.
   - Conda envs, from `envs/<name>.yaml` (or `envs/<name>.explicit.txt` for exact builds):
     ```
     conda env create -f envs/deseq2.yaml
     conda env create -f envs/haydar_methylation.yaml
     conda env create -f envs/data_analysis.yaml
     conda env create -f envs/snakemake.RNAseq.standard.with.SQuiRE.yaml
     ```
   - Install the R packages in `envs/R_packages_non_conda.txt` with BiocManager (Bioconductor 3.22). The
     per-rule pipeline envs are built by Snakemake from `workflow/01_rnaseq_pipeline/envs/`.
   - `run_all.sh` sources `~/.bashrc` for `conda activate`, and loads Apptainer with `module load` (module name
     in `slurm.apptainer_module`).
4. **Cluster.** Set the SLURM account/partition (and optional email) under `slurm:`. The scripts assume SLURM. On
   another scheduler, run stage 1's `snakemake` command and `workflow/run_stage.sh 2…5` by hand.

## Running

```
./run_all.sh          # everything: stage 1, then 2 -> 5 (3 runs in parallel with 1-2)
./run_all.sh 1        # FASTQ -> counts only (~1-2 days, ~450 GB in results/rnaseq)
./run_all.sh 2-5      # downstream only, on results/rnaseq/.../squire_gene_counts.tsv
RNASEQ_COUNTS=/path/to/squire_gene_counts.tsv ./run_all.sh 2-5   # downstream on another count matrix
START_AT=06_dmrcate.R REPO_ROOT=$PWD workflow/run_stage.sh 3     # resume one stage at a script
SLURM_EXCLUDE=cpu096 ./run_all.sh 1   # keep stage-1 rule jobs off a bad node; Snakemake reruns only what's missing

# reproduce the paper's numbers exactly, from the deposited intermediates (see Reproducibility)
RNASEQ_COUNTS=/path/to/squire_gene_counts.tsv METH_PREPROCESSED=/path/to/meth_qc_dir ./run_all.sh 2-5
```

### Reproducibility

Every step is deterministic, so a rerun of this repo gives the same output twice. Two steps of the original
analysis were not, though, so a fresh run from raw data matches the paper's conclusions but not every digit.
The paper's numbers are reproduced exactly by starting from the two deposited intermediates (command above).
All of this was checked on 2026-10-20/21; see `verification/REPORT.md`.

- **SQuIRE gene counts (stage 1).** In the original run, SQuIRE Count used 8 threads, which only StringTie uses.
  Multithreaded StringTie splits reads among overlapping transcripts in a thread-order-dependent way, so
  two runs differed in 42 of 1.7 M gene-count cells (23 genes). The TE counts were identical. The repo now runs
  Count with 1 thread. Deposited: `squire_gene_counts.tsv` (GEO processed data) → `RNASEQ_COUNTS`.
  - Effect of a fresh FASTQ run: the DE calls are unchanged (111 genes) and so are the 26 shared GO terms.
    RNA GSEA goes from 679 to 689 terms, and Fig 2D reads 609 / OR 7.2 / p 1.4e-10 instead of 600 / 7.3 / 9.6e-11.
- **Methylation preprocessing (stage 3, `01_load_idats_qc.R`).** sesame's `matchDesign` step is stochastic. The
  original run set no seed; the repo now calls `set.seed(1)` before each sample's `prepSesame()`. Deposited:
  `sdfs_QCDPBM.rds` and `betas_QCDPBM.rds` → `METH_PREPROCESSED` (or `methylation_preprocessed_override` in
  `config/paths.yaml`), which 01 then copies in instead of reading the IDATs.
  - Effect of a fresh IDAT run: 517 filtered DMRs instead of 507, and 452 DMR genes instead of 443. The 26
    shared GO terms are unchanged.
  - The betas are made with `openSesame(sdfs, func = getBetas)`. This applies sesame's default `QCDPB` prep a second
    time on top of the `QCDPBM` SigDFs. The step is deterministic, and the repo keeps it because the paper's
    betas were made this way.
- **NS DESeq2/apeglm** shows floating-point noise across machines (log2FC ≤ 1e-4). No call or figure changes.

Stage logs go to `logs/stage<N>_<jobid>.log`. Stage 1's Snakemake logs are in `results/rnaseq/logs/`.

| Stage | Scripts, in order (env) | Main outputs |
|---|---|---|
| 1 RNA-seq pipeline | Snakefile: FastQC → Trim Galore → STAR + Picard + MultiQC (QC) → SQuIRE Build/Clean/Map/Count → `combine_squire_counts.py` | `results/rnaseq/output/results/counts/squire_gene_counts.tsv` |
| 2 RNA-seq DE (deseq2) | `run_deseq2_gene_NS.R`, `run_deseq2_gene_S1.R`, `run_gsea_gene_NS.R` | `results/rnaseq_de/gene_{NS,S1}/` (`results_full.csv`, `dds.rds`), `gene_NS/enrichment_gsea/` |
| 3 Methylation (haydar_methylation) | `00_build_sample_sheet.R`, `00b_build_mm39_annotation.R`, `01_load_idats_qc.R`, `06_dmrcate.R`, `06b_gene_background.R`, `11_dmgsea.R`, `13_gsk_vs_aza.R` | `results/methylation/{metadata,reference/mm39,qc,dm,dmgsea,gsk_vs_aza}/` |
| 4 Integration | `12_dmgsea_plots.R` (deseq2), `01_rna_donor_slopes.R` (deseq2), `02_methylation_donor_slopes.R` (haydar_methylation), `03_correlate_and_plot.R` (deseq2), `10_candidate_gene_rna_data.R` (deseq2), `11_candidate_gene_meth_data.R`, `12_candidate_gene_heatmaps.R` (haydar_methylation) | `results/methylation/dmgsea/meth_rna_shared_GO_BP.csv`, `results/integration/` |
| 5 Figures | `fig2_01` (deseq2), `fig2_02` (haydar_methylation), `fig3_01` (deseq2), `fig3_02`, `fig3_04`, `figS_gsk_02` (haydar_methylation); then `fig2_03_plot.py` (default, and `--font-bump 1 --tag v3`), `fig3_03_plot.py` (default, `--absolute`, `--bars`, `--absolute --bars`), `figS_gsk_03_plot.py` (data_analysis) | `figures/` |

Each script's header records the analysis choices behind it (dated, as agreed while the analysis was done).

## Figures → code → outputs

| Figure | Panels | Data prep → plot | Output |
|---|---|---|---|
| Fig 2, RNA-seq × methylation | A DE-gene heatmap (111 genes); B DMR heatmap (507 DMRs); C DE genes vs DMR genes; D GO:BP RNA GSEA vs methylation dmGsea; E the 26 shared terms (Wang tree + dot plot); F effector-function candidate genes | `fig2_01_prep_rna.R`, `fig2_02_prep_meth.R` → `fig2_03_plot.py` | `figures/fig2_rna_methylation/release_v3/` (as shared, `--font-bump 1 --tag v3`); working copy + `panels/` one level up |
| Fig 3, cytokines | IL-2, IFN-γ, TNF-α: RNA (NS and S1) and promoter/body methylation vs dose | `fig3_01_prep_rna.R`, `fig3_02_prep_meth.R`, `fig3_04_genomewide_fdr.R` (genome-wide FDRs quoted in the legend) → `fig3_03_plot.py` | `figures/fig3_cytokines/Fig3_cytokines_meth_rna_{FC,absolute,FC_bars,absolute_bars}.{pdf,svg,png}` |
| Supplementary, GSK-3484862 | top-variable-probe and 507-DMR heatmaps, GSK vs AZA slopes, limma/DMRcate slope contrast | `13_gsk_vs_aza.R` → `figS_gsk_02_prep.R` → `figS_gsk_03_plot.py` | `figures/figS_gsk_methylation/` |

Legend numbers (gene/DMR/term counts, FDRs, odds ratios) are printed by the prep and plot scripts and land in
`logs/stage5_*.log`. `verification/REPORT.md` lists them.


