#!/usr/bin/env Rscript
# Per-gene RNA-seq (dose_nM slope of log2 expression) vs methylation
# (median over the gene's promoter + gene-body probes of each probe's dose_nM
# beta slope; see 02) correlation, NS/AZA-treated only.
#
# Two donor scopes (EV07-only, average-across-donors; see 01/02) x three
# gene-set splits:
#   - All genes: every gene with both an RNA slope and >=1 mapped probe
#     (promoter 1.5kb + body, dmGsea mapping; was ">=1 DMR within 5kb" before
#     2026-09-24).
#   - DE genes: restricted to DE_sig==TRUE (existing pooled 5-donor DESeq2
#     padj<0.05 call from results_full.csv).
#   - Sig DMR genes: restricted to is_sig_DMR_gene==TRUE (>=1 promoter/body
#     probe inside one of the 507 FDR+|meandiff|>0.10 DMRs from 06_dmrcate.R --
#     the dmGsea hit-gene definition).
# The plotted RNA/methylation VALUE is always the freshly-computed
# per-donor-slope (EV07 or avg-across-donors) -- DE_sig/is_sig_DMR_gene are
# used only to pick which genes are included, not what value is shown.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
in_dir <- out_dir <- INT_DIR
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Load + join ----
rna <- fread(file.path(in_dir, "rna_gene_slopes.csv"))
meth <- fread(file.path(in_dir, "methylation_gene_slopes.csv"))

dt <- merge(rna, meth, by = "symbol")
cat(sprintf("Joined per-gene table: %d genes (RNA universe %d, methylation universe %d)\n",
            nrow(dt), nrow(rna), nrow(meth)))
fwrite(dt, file.path(out_dir, "gene_level_multiomic_table.csv"))

## ---- 2. Define gene-set filters + donor-scope column pairs ----
gene_sets <- list(
  "All genes"     = quote(rep(TRUE, nrow(dt))),
  "DE genes"      = quote(dt$DE_sig),
  "Sig DMR genes" = quote(dt$is_sig_DMR_gene)
)
donor_scopes <- list(
  "EV07"        = c(rna = "rna_slope_EV07",        meth = "meth_slope_EV07"),
  "avg_donors"  = c(rna = "rna_slope_avg_donors",  meth = "meth_slope_avg_donors")
)

## ---- 3. Compute correlations, collect plotting data ----
summary_rows <- list()
plot_rows <- list()

for (gs_name in names(gene_sets)) {
  keep <- eval(gene_sets[[gs_name]])
  for (ds_name in names(donor_scopes)) {
    cols <- donor_scopes[[ds_name]]
    sub <- dt[keep, .(x = get(cols["rna"]), y = get(cols["meth"]))]
    sub <- sub[is.finite(x) & is.finite(y)]

    pear <- cor.test(sub$x, sub$y, method = "pearson")
    spear <- suppressWarnings(cor.test(sub$x, sub$y, method = "spearman"))

    summary_rows[[paste(gs_name, ds_name)]] <- data.table(
      gene_set = gs_name, donor_scope = ds_name, n = nrow(sub),
      pearson_r = unname(pear$estimate), pearson_p = pear$p.value,
      spearman_rho = unname(spear$estimate), spearman_p = spear$p.value
    )
    plot_rows[[paste(gs_name, ds_name)]] <- sub[, `:=`(gene_set = gs_name, donor_scope = ds_name)]

    cat(sprintf("%-14s x %-11s: n=%5d  Pearson r=%6.3f (p=%.2e)  Spearman rho=%6.3f (p=%.2e)\n",
                gs_name, ds_name, nrow(sub), pear$estimate, pear$p.value,
                spear$estimate, spear$p.value))
  }
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(out_dir, "correlation_summary.csv"))
cat("\nSaved", file.path(out_dir, "correlation_summary.csv"), "\n")

## ---- 4. Faceted scatter plot: gene_set (rows) x donor_scope (columns) ----
plot_dt <- rbindlist(plot_rows)
plot_dt$gene_set <- factor(plot_dt$gene_set, levels = names(gene_sets))
plot_dt$donor_scope <- factor(plot_dt$donor_scope, levels = names(donor_scopes))

label_dt <- merge(summary_dt, plot_dt[, .(x_max = max(x), y_max = max(y)), by = .(gene_set, donor_scope)],
                   by = c("gene_set", "donor_scope"))
label_dt$label <- sprintf("r=%.2f\nrho=%.2f\nn=%d", label_dt$pearson_r, label_dt$spearman_rho, label_dt$n)

p <- ggplot(plot_dt, aes(x, y)) +
  geom_point(alpha = 0.35, size = 0.7, color = "#2A78D6") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 0.6) +
  geom_text(data = label_dt, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.1, size = 3, lineheight = 0.9) +
  facet_grid(gene_set ~ donor_scope, scales = "free") +
  labs(x = "RNA-seq dose slope (log2 expr / nM AZA)",
       y = "Methylation dose slope (median probe beta / nM AZA)",
       title = "RNA-seq vs. methylation AZA dose-response slopes, per gene",
       subtitle = "NS (non-stimulated), AZA-treated; methylation = median over promoter + gene-body probes\nrows = gene set, columns = donor scope") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(plot_dir, "gene_level_correlation_panels.png"), p, width = 9, height = 10, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "gene_level_correlation_panels.pdf"), p, width = 9, height = 10, device = cairo_pdf)
cat("Saved", file.path(plot_dir, "gene_level_correlation_panels.png/.pdf"), "\n")
cat("=== Done ===\n")
