# Verification: publication repo rerun vs the current project

Run on 2026-10-20 on the GWU HPC cluster. Reference = the outputs in the current project (`Haydar_RNAseq/`,
`Haydar_Methylation_Array/`, `integration/`) that made the shared figures (Fig 2 release_v3, Fig 3, Fig S GSK).
Files compared = the 261-row `manifest.tsv` (figure-path outputs, from the inventory agent).
Per-file results: `output_comparison.tsv` (`compare_outputs.R`) and `figure_comparison.tsv` (`compare_figures.py`).

## Summary

| What | Result |
|---|---|
| **All 59 figure files**, on the current count matrix (Fig 2 working + release_v3 composites and panels; Fig 3 FC/absolute/bars variants; Fig S composite + panels; PDF/SVG/PNG) | **59/59 PASS**: 32 PDF/PNG pixel-identical at 300 dpi; 15 SVGs text-identical after id/date canonicalization; 12 SVGs whose text differs only in float coordinates have pixel-identical PDF twins. Checked 2026-10-20 22:55 with the corrected `compare_figures.py` (see the note below) |
| Figures on the rerun's own count matrix (Finding 3) | 31 PASS / 28 FAIL. Fig 2A, 2D, 2E, the Fig 2 composites and all Fig 3 variants change slightly; Fig 2B/C/F and Fig S are identical |
| Data files, stages 2–5 (77) | 66 byte-identical or identical content; 3 within 1e-8 relative; 7 differ by floating-point noise traced to the NS DESeq2/apeglm fit (below); 1 TxDb sqlite not compared |
| Legend numbers | all reproduced (table below) |
| Stage 1, FASTQ → counts | SQuIRE subfamily and locus (TE) matrices **byte-identical**. Gene matrix: 1,710,258 of 1,710,300 cells identical; **42 cells (23 genes) differ**, caused by multithreaded StringTie inside SQuIRE Count (Finding 3). Downstream: same DE calls; GSEA 679 → 689 terms; the Fig 2D legend changes |
| Methylation from raw IDATs | **not exactly reproducible**: stochastic preprocessing in the original `01_load_idats_qc.R` (below). Everything from the preprocessed SigDFs onward reproduces exactly |

How the run was done:
- Stages 2–5 used the current count matrix (`RNASEQ_COUNTS`), so the downstream checks didn't wait for stage 1.
- Stage 3 ran twice:
  1. from the IDATs (jobs 138531600–602; snapshot in `idat_rerun_unseeded/`);
  2. from `06_dmrcate.R` onward on the original `sdfs_QCDPBM.rds` / `betas_QCDPBM.rds` (jobs 138532951–53).

  The comparisons above are for run 2.

> **Correction (2026-10-20 22:45).** Until 22:45 the figure checker reported every PDF/PNG as pixel-identical,
> whatever the content. `pixel_diff` called `Image.getbbox()` on an RGBA difference image, and Pillow then
> looks at the alpha channel only. It now compares all channels with numpy. The fix was checked on a known-different
> pair (the old vs new Fig 2D Venn, 574 vs 583 terms): it FAILs that pair and PASSes identical renders. The 59/59
> result above comes from a rerun of stages 2/4/5 on the current matrix (jobs 138533790–92), checked with the fixed
> script. The earlier 59/59 claim was not valid evidence; the corrected check confirms the same conclusion.

## Finding 1: methylation preprocessing can't be regenerated exactly from IDATs

Run 1 (from IDATs) gave a different result set:

| | paper | IDAT rerun |
|---|---|---|
| tested probes (06) | 277,292 | 278,805 |
| DMRs (all / \|meandiff\| > 0.10) | 25,372 / 507 | 25,875 / 517 |
| DMR-hit genes | 443 | 452 |
| dmGsea GO:BP FDR < 0.05 | 48 | 50 |
| terms shared with RNA GSEA | 26 | 26 |
| GSK contrast: probes FDR < 0.05 / DMRs / filtered | 11,739 / 3,165 / 616 | 11,939 / 3,192 / 571 |

Causes. Both are in the original script; the repo did not introduce them.

1. **sesame `matchDesign` (prep step "M") is random, and the script sets no seed.**
   - Two `prepSesame(prep = "QCDPBM")` runs on the same IDAT differ by up to 0.012 in beta.
   - Steps Q–B are deterministic.
   - With `set.seed(1)` before each run, the two runs are identical.
   - sesame 1.28.1 / sesameData 1.28.0 are unchanged since the env was created (2026-09-03), so a package update
     is not the cause.
2. **The saved betas are preprocessed twice.** `openSesame(sdfs, func = getBetas)` applies its default
   `prep = "QCDPB"` to SigDFs that are already QCDPBM-processed. So the betas behind every methylation result are
   QCDPBM followed by QCDPB.
   - This second pass is deterministic: `openSesame(original sdf)` reproduces the original betas exactly.
   - It moves betas by up to 0.81 compared with a single QCDPBM pass (`getBetas(sdf)`).
   - The script header, and so far the methods notes, describe a single QCDPBM pass.

Run-to-run size of cause 1: 507 → 517 filtered DMRs, 443 → 452 DMR genes, 48 → 50 GO terms; the 26 shared terms
are stable.

**Decision needed (not made here):**
- **(a) Keep the published results.** Describe the actual preprocessing in the methods, and deposit
  `sdfs_QCDPBM.rds` (or the betas) as the starting point for reproduction.
- **(b) Fix the script.** Add `set.seed`, and/or call `openSesame(..., prep = "")` or `getBetas` directly. Then
  rerun the methylation side; all methylation numbers and Fig 2B–F / Fig 3 / Fig S values shift slightly.

## Finding 2: NS DESeq2 floating-point noise (no effect on any call or figure)

The `gene_NS` refit matches the original in everything that determines the results:
- **Identical:** counts, size factors, gene-wise and fitted dispersions (≤ 1e-12 relative), the 712 padj < 0.05
  genes, the 111 `sig` genes, and the 679 GSEA terms at padj ≤ 0.05.
- **Small differences:**
  - final MAP dispersions: ≤ 5e-5 relative;
  - apeglm log2FC: ≤ 9.4e-5 absolute;
  - padj: ≤ 2e-6 absolute;
  - fgsea p-values: 4–5 terms move by up to 19% relative (≤ 1.6e-5 absolute), with no term crossing 0.05.
- **Likely cause:** multithreaded BLAS / node differences.
- **Propagation:** into `rna_gene_slopes.csv`, `gene_level_multiomic_table.csv`, Fig 2 `de_heatmap_matrix.csv`
  (≤ 7e-8) and Fig 3 `rna_unstimulated_stats.csv` (≤ 3e-6). No figure pixel changes.
- `gene_S1` is byte-identical.

## Finding 3: SQuIRE gene counts are not exactly reproducible (multithreaded StringTie)

Stage 1 (job 138531605, FASTQ → counts) compared with `Haydar_RNAseq/output/results/counts/`:
- `squire_subfamily_counts.tsv` and `squire_locus_counts.tsv`: byte-identical.
- Manifest rows for stage 1 (`output_comparison_stage1_counts.tsv`): SQuIRE clean/fetch outputs exact; the 30 per-sample `*_refGenecounts.txt` differ; Trim Galore FASTQs, SQuIRE BAMs and the STAR index were not compared (large binaries).
- `squire_gene_counts.tsv`: same 57,010 genes × 30 samples. 42 cells differ, in 23 genes and 20 samples (23 cells in
  16 genes are NS); the full list is in `stage1_gene_count_diffs.tsv`.
  - 29 cells differ by 1–5 counts.
  - The larger differences are in genes that overlap other genes or have many paralogs:
    - Gm54610, a miRNA inside Rn7s2: S1_20_50 57,421 → 2;
    - Actg1: NS_18_200 76,184 → 103,071;
    - Rps27a: NS_07_0 +1,882;
    - Rbm3: +715 / +1,265;
    - Uba52: −1,016;
    - Atp5l: −158;
    - Gm10275: −304.
- **Cause:** SQuIRE Count's TE outputs (`*_TEcounts.txt`, `*_subFcounts.txt`) are identical per sample. Only the
  StringTie gene-abundance outputs differ (`*_abund.txt`, `*.gtf`, `*_refGenecounts.txt`). SQuIRE runs StringTie
  1.3.3 `-e` with `-p` = the pipeline's 8 threads (`squire/Count.py`, `Stringtie()`). Multithreaded StringTie
  splits reads among overlapping transcripts in a thread-order-dependent way.
- **Upstream steps:** Trim Galore output and SQuIRE Map BAMs are not checked here. The identical TE counts show that
  the alignments agree.
- **STAR (QC only):** 26 rule attempts failed on node cpu096 with STAR `std::bad_alloc`. Four of them failed on their
  retry too: NS_07_100, NS_18_0, S1_07_0, S1_20_0. They are being redone with `SLURM_EXCLUDE=cpu096 ./run_all.sh 1`
  (job 138533767). These jobs don't feed the counts.
- **Effect downstream.** Stages 2, 4 and 5 were rerun on the new matrix (jobs 138533772 → 75 → 76). Stage 3 was
  not rerun (original-SigDF outputs), which isolates the count change. Snapshot: `stage1_counts_run/`. Tables:
  `{output,figure}_comparison_stage1_counts.tsv`.
  - Size factors move by ~1e-5 relative, so every DESeq2 statistic shifts a little: median |ΔLFC| 1e-5, 99th
    percentile 0.002, max 0.065.
  - NS: 712 → 713 genes at padj < 0.05 (0 lost, 1 gained). The **111 sig genes are identical**.
  - S1: 617 / 53, identical.
  - GO:BP GSEA: 679 → 689 terms at padj ≤ 0.05 (8 lost, 18 gained). The **26 shared terms are identical** (same IDs).
  - Legend numbers:
    - Fig 2C is unchanged.
    - **Fig 2D changes:** 609 RNA-seq terms (was 600), OR = 7.2 (was 7.3), p = 1.4e-10 (was 9.6e-11). The Venn's
      RNA-only count changes from 574 to 583.
    - Fig 2E: same 26 terms and 5 clusters (6/5/4/8/3); the RNA FDR values move.
    - Fig 2A: same 111 genes; the values move slightly.
    - Fig 3: the RNA values move by ≤ 1e-4 relative.
    - Fig S: identical.
  - So the conclusions don't change, but a figure or legend made from a fresh FASTQ run would not match the paper
    figures digit for digit.
- **Options:**
  - **(a)** Deposit the count matrix with the paper (GEO processed data) and document stage 1 as reproducible up to
    StringTie thread noise.
  - **(b)** Change SQuIRE Count to `-p 1` for the published pipeline. Its StringTie step is fast compared with the
    TE EM. The paper figures would still come from the existing matrix, so in practice this is (a) plus a
    deterministic pipeline for others.

## Legend numbers (reproduced, run 2 logs: `logs/stage{3,4,5}_13853295*.log`)

| Number | Value | Where |
|---|---|---|
| DE genes (NS, padj < 0.05 & \|log2FC\| > 0.5) | 111 (92 up, 19 down); 712 at padj < 0.05 | stage 2, Fig 2A |
| DMRs | 25,372 called; 507 with \|meandiff\| > 0.10, all hypomethylated | 06, Fig 2B |
| DMR-hit genes (dmGsea mapping) | 443; 28,546 Entrez genes with ≥ 1 probe; 1,063 sigProbes | 11, Fig 3 legend |
| Liftover vs native mm39 manifest | 0.0011% disagreement among uniquely-mapping probes | 06b |
| Fig 2C | 13,782 genes on both platforms; 102 DE, 385 DMR, 4 both; one-sided Fisher OR = 1.42, p = 0.32 | fig2_03 |
| Fig 2D | 4,166 GO:BP terms in both; 600 RNA, 48 methylation, 26 both; OR = 7.3, p = 9.6e-11 | fig2_03 |
| Fig 2E | 26 shared terms (all RNA-up), 5 clusters (6/5/4/8/3) | fig2_02 |
| Fig 3 genome-wide | 28,546 genes; FDR < 0.05: 16,393 OLS / 17,523 moderated; 94.9% lose methylation; median Δβ −0.0248 | fig3_04 |
| Fig 3 genes (Δβ 0→200 nM, p) | Il2 −0.030 (0.0066), Ifng −0.041 (0.0044), Tnf −0.003 (0.66); none in a DMR | fig3_02 |
| Fig S | 507-DMR mean Δβ at 200 nM: AZA −0.121/−0.146/−0.109, GSK 0.001; GSK contrast 11,739 probes FDR < 0.05, 3,165 DMRs (616 filtered, all GSK-higher; 316 overlap the 507) | figS_gsk_02, 13 |

Every one of these that is quoted in `Fig2_legend_notes.md` or the Fig 3 / Fig S legend notes/drafts in `integration/paper_figures/` matches (checked: Fig 2C/D counts, OR, p; 28,546; 11,739; 3,165; 616; 316; −0.121). The Fig 3 genome-wide FDR counts are not quoted there.

## Decisions (user, 2026-10-21)
- **Findings 1 and 3: option (a)** plus deterministic code. The paper keeps its published numbers. The deposited
  intermediates, the SQuIRE gene count matrix and `sdfs_QCDPBM.rds` / `betas_QCDPBM.rds`, are the documented
  start for an exact reproduction (`RNASEQ_COUNTS`, `METH_PREPROCESSED`; README, Reproducibility).
  - `01_load_idats_qc.R` calls `set.seed(1)` before each `prepSesame()`. Tested: two runs on 3 IDATs give identical
    betas (max diff 0). The double `QCDPB` pass in `openSesame` is kept and documented.
  - `METH_PREPROCESSED` copy branch tested: md5s are unchanged.
  - SQuIRE Count runs with 1 thread (`cpus_per_task = 1`). In `Count.py`, `-p` is used only by StringTie.
    Not rerun: this is single-threaded StringTie.
- **Fonts:** Arial is not redistributed (`.gitignore`). The scripts fall back to the shipped Liberation Sans.
- **Names:** the collaborator's name was removed from the code comments.

## Still open
- Methods section.
- The S3/SMB paths in the RNA manifest, LICENSE, GEO accessions.
