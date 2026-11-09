#!/usr/bin/env python3
"""Supplementary figure: GSK-3484862 vs 5-azacytidine methylation (B7H3 CAR-T, MM285).

Reuses the Fig 2 style and helpers (fig2_03_plot.py:
palette, Ward row clustering, annotation strips, vector colorbars, Arial), so these heatmaps match Fig 2B.

Panels (one plot_* function each, drawing into one `ax` container; the composite is drawn in one pass):
  A  plot_a_top500          500 most variable CpG probes across the GSK donor's 4 doses, row z-score of beta
                            (rows clustered only; columns in dose order). Regenerates the collaborators' heatmap.
  B  plot_b_dmr507_delta    the 507 AZA DMRs (Fig 2B rows): region delta beta vs the same donor's 0 nM, AZA donors
                            and the GSK donor side by side on one fixed scale (not z-scored, so a flat response
                            stays white)
  C  plot_c_top500_delta    unsupervised comparison: top 500 variable autosomal CpGs after donor centring (each
                            donor's mean removed), AZA and GSK pooled; shown as delta beta vs the donor's 0 nM on its
                            own +-0.3 scale. Probes are not chosen by any AZA statistic, so GSK-specific changes would
                            surface here if there were any.
  D  plot_d_slope_scatter   GSK vs AZA at the 507 AZA DMRs (13_gsk_vs_aza.R part 1, descriptive): per-DMR slope of
                            region beta on dose (delta beta per 100 nM), mean of the 3 AZA donors (x) vs the GSK donor
                            (y), y = x line; equal axes
  E  plot_e_contrast_volcano probe-level limma contrast dose_GSK - dose_AZA (~ Donor + dose_AZA + dose_GSK, M-values;
                            13_gsk_vs_aza.R part 2): beta-scale slope difference over 0 -> 200 nM vs -log10 FDR, plus
                            DMRcate counts. The point layer is rasterized (262k probes); axes and text stay vector.
  F  plot_f_dmr_volcano      DMR-level version of E: every DMRcate DMR on the same contrast, meandiff (beta-scale
                            slope difference over 0 -> 200 nM) vs -log10 min smoothed FDR (the 08_dmr_pathway_plots.R
                            DMR volcano convention); |meandiff| > 0.10 highlighted, dashed lines at +-0.10. All DMRs are
                            FDR-significant by construction, so there is no FDR line.
  G  plot_g_sig_probe_slopes  the 11,739 FDR < 0.05 probes of E: model AZA dose slope (pooled over 3 donors, x) vs
                            GSK dose slope (donor 8, y), beta per 100 nM; y = x line, equal axes; point layer
                            rasterized
  H  plot_h_dmr_slopes      the 3,165 DMRs of F (all FDR-significant): mean probe slope per DMR, AZA (x) vs GSK (y);
                            the 616 with |meandiff| > 0.10 highlighted


Outputs: paper_figures/figS_gsk_methylation/FigS_gsk_methylation.{pdf,svg} (+ _preview.png) and
panels/FigS_gsk<letter>_<name>.{pdf,svg}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_03_plot as f2  # noqa: E402  (style + helpers; its figure code only runs under __main__)

OUT = f2.FIG_DIR / "figS_gsk_methylation"
DATA = OUT / "data"
PANEL_DIR = OUT / "panels"
STEM = "FigS_gsk_methylation"

# The GSK donor (EV03) gets the next donor number after Fig 2's 1-7 (confirm with the user) and its own hue.
DONOR_COLORS = {**f2.DONOR_COLORS, "EV03": "#56606b"}
DONOR_LABELS = {**f2.DONOR_LABELS, "EV03": "8"}
DRUG_COLORS = {"AZA": "#1c5cab", "GSK": "#eda100"}
DRUG_NAMES = {"AZA": "5-Azacytidine", "GSK": "GSK-3484862"}
DELTA_LIM = 0.2      # panel B: region (DMR) mean delta
DELTA_LIM_C = 0.3    # panel C: single-probe delta, larger than region means (|delta| 95th pct 0.26, 98th 0.29)
fs = f2.fs


# --------------------------------------------------------------- helpers ----
def read(prefix):
    return pd.read_csv(DATA / f"{prefix}_matrix.csv"), pd.read_csv(DATA / f"{prefix}_samples.csv")


def order_cols(meta, donor_order=("EV07", "EV08", "ZY03", "EV03")):
    meta = meta.copy()
    meta["drug_rank"] = meta["drug"].map({"AZA": 0, "GSK": 1})
    meta["donor_rank"] = meta["donor"].map({d: i for i, d in enumerate(donor_order)})
    return meta.sort_values(["drug_rank", "dose_nM", "donor_rank"]).reset_index(drop=True)


def legends(ax, items):
    """Stacked legends in one axes, top-left anchored. items = [(title, [(label, color, outlined)])]."""
    ax.set_axis_off()
    kw = dict(frameon=False, title_fontsize=fs(5.5), fontsize=fs(5.5), handlelength=0.9, handleheight=0.9,
              handletextpad=0.4, labelspacing=0.3, borderaxespad=0, borderpad=0, alignment="left")
    h_pts = f2.container_height_pts(ax)
    drop = 0.0
    for title, entries in items:
        leg = ax.legend(handles=[Patch(facecolor=c, edgecolor=f2.INK2 if o else "none", linewidth=0.3, label=l)
                                 for l, c, o in entries],
                        loc="upper left", bbox_to_anchor=(0, 1 - drop / h_pts), title=title, **kw)
        ax.add_artist(leg)
        drop += 1.33 * fs(5.5) * (len(entries) + 1) + 7


def dose_entries(doses):
    return [(str(d), f2.DOSE_COLORS[d], d == 0) for d in doses]


def heatmap(ax, values, meta, strips, title, cbar_label, legend_items, vmin=-2, vmax=2, ticks=(-2, 0, 2),
            subtitle=None):
    """Fig 2B layout: row dendrogram (Ward) | matrix | legends, annotation strips on top, colorbar below."""
    f2.container(ax)
    x_dend, w_dend = 0.0, 0.10
    x_hm, w_hm = 0.11, 0.60
    y_hm, h_hm = 0.12, 0.70
    strip_h, strip_gap = 0.028, 0.008
    if values.shape[0] >= 2:
        Z = f2.cluster_rows(values)
        values = values[leaves_list(Z)]
        f2.draw_dendrogram(f2.inset(ax, x_dend, y_hm, w_dend, h_hm), Z, values.shape[0])
    ax_hm = f2.inset(ax, x_hm, y_hm, w_hm, h_hm)
    f2.draw_matrix(ax_hm, values, vmin=vmin, vmax=vmax)
    f2.dose_block_gaps(ax_hm, meta, values.shape[0])
    y = y_hm + h_hm + strip_gap
    for key, colors, label in strips:
        f2.annotation_strip(f2.inset(ax, x_hm, y, w_hm, strip_h), meta, key, colors, label, label_side="left")
        y += strip_h + strip_gap
    f2.panel_title(ax, title, x=x_hm + w_hm / 2)
    if subtitle:
        ax.annotate(subtitle, xy=(x_hm + w_hm / 2, 1), xycoords="axes fraction", xytext=(0, -1.45 * fs(7)),
                    textcoords="offset points", ha="center", va="top", fontsize=fs(5.5), color=f2.INK2)
    legends(f2.inset(ax, x_hm + w_hm + 0.04, y_hm, 0.24, y - y_hm), legend_items)
    f2.z_colorbar(f2.inset(ax, x_hm + w_hm * 0.2, 0.065, w_hm * 0.6, 0.022), label=cbar_label, vmin=vmin,
                  vmax=vmax, ticks=ticks, orientation="horizontal")
    return ax


# --------------------------------------------------------------- panels ----
def plot_a_top500(ax):
    mat, meta = read("top500_gsk")
    meta = order_cols(meta)
    z = f2.row_z(mat[meta["sample"]].to_numpy())
    return heatmap(ax, z, meta, [("dose_nM", f2.DOSE_COLORS, "Dose")],
                   "GSK-3484862: top 500 variable CpGs", "Row z-score of β",
                   [("GSK (nM)", dose_entries([0, 50, 100, 200]))],
                   subtitle=f"Donor {DONOR_LABELS['EV03']}, B7H3 CAR")


def plot_b_dmr507_delta(ax):
    mat, meta = read("dmr507_delta")
    meta = order_cols(meta[meta["dose_nM"] > 0])  # 0 nM is the reference (delta = 0 by definition)
    vals = mat[meta["sample"]].to_numpy()
    donors = ["EV07", "EV08", "ZY03", "EV03"]
    return heatmap(ax, vals, meta,
                   [("donor", DONOR_COLORS, "Donor"), ("drug", DRUG_COLORS, "Drug"), ("dose_nM", f2.DOSE_COLORS, "Dose")],
                   f"AZA DMRs (n = {len(mat)}): AZA vs GSK", "Δβ vs same donor's 0 nM",
                   [("Donor", [(DONOR_LABELS[d], DONOR_COLORS[d], False) for d in donors]),
                    ("Drug", [(DRUG_NAMES[k], c, False) for k, c in DRUG_COLORS.items()]),
                    ("Dose (nM)", dose_entries([50, 100, 200]))],
                   vmin=-DELTA_LIM, vmax=DELTA_LIM, ticks=(-DELTA_LIM, 0, DELTA_LIM))


def plot_c_top500_delta(ax):
    mat, meta = read("top500_delta")
    meta = order_cols(meta[meta["dose_nM"] > 0])
    vals = mat[meta["sample"]].to_numpy()
    donors = ["EV07", "EV08", "ZY03", "EV03"]
    return heatmap(ax, vals, meta,
                   [("donor", DONOR_COLORS, "Donor"), ("drug", DRUG_COLORS, "Drug"), ("dose_nM", f2.DOSE_COLORS, "Dose")],
                   "Top 500 variable CpGs (donor-centred), AZA + GSK", "Δβ vs same donor's 0 nM",
                   [("Donor", [(DONOR_LABELS[d], DONOR_COLORS[d], False) for d in donors]),
                    ("Drug", [(DRUG_NAMES[k], c, False) for k, c in DRUG_COLORS.items()]),
                    ("Dose (nM)", dose_entries([50, 100, 200]))],
                   vmin=-DELTA_LIM_C, vmax=DELTA_LIM_C, ticks=(-DELTA_LIM_C, 0, DELTA_LIM_C))


def _axes(ax, rect):
    """A plain axes inset in a panel container (fractions of the container), Fig 2 tick style."""
    a = f2.inset(ax, *rect)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
    a.tick_params(labelsize=fs(5.5), width=0.5, length=2)
    return a


def plot_d_slope_scatter(ax):
    f2.container(ax)
    d = pd.read_csv(DATA / "gsk_vs_aza_dmr507_slopes.csv")
    x, y = d["slope_per100nM_AZA_mean"], d["slope_per100nM_GSK_EV03"]
    a = _axes(ax, (0.17, 0.14, 0.66, 0.70))
    lo = np.floor(min(x.min(), y.min()) / 0.05) * 0.05
    hi = max(0.05, np.ceil(max(x.max(), y.max()) / 0.05) * 0.05)
    a.plot([lo, hi], [lo, hi], ls=(0, (3, 2)), lw=0.6, color=f2.INK2, zorder=1)
    a.axhline(0, lw=0.4, color="#d6d5d0", zorder=0)
    a.axvline(0, lw=0.4, color="#d6d5d0", zorder=0)
    a.scatter(x, y, s=4, color=DRUG_COLORS["GSK"], edgecolors="#2e2d2a", linewidths=0.15, alpha=0.85, zorder=3)
    a.set_xlim(lo, hi); a.set_ylim(lo, hi); a.set_aspect("equal")
    a.set_xlabel("AZA, mean of donors 1, 6, 7\n(Δβ per 100 nM)", fontsize=fs(6))
    a.set_ylabel(f"GSK, donor {DONOR_LABELS['EV03']} (Δβ per 100 nM)", fontsize=fs(6))
    n, above, above_all = len(d), int(d["gsk_above_aza_mean"].sum()), int(d["gsk_above_all_aza"].sum())
    # summary sits below the y = x line, where there are no points
    a.text(0.97, 0.20, f"{above} of {n} DMRs above y = x\n"
           f"Median: AZA {x.median():.3f},\nGSK {y.median():+.3f}",
           transform=a.transAxes, ha="right", va="bottom", fontsize=fs(5.5), color=f2.INK2, linespacing=1.3)
    a.text(hi - 0.01 * (hi - lo), hi - 0.12 * (hi - lo), "y = x", ha="right", va="top", fontsize=fs(5),
           color=f2.INK2)
    f2.panel_title(ax, "Dose slopes at the 507 AZA DMRs", x=0.5)
    return ax


def plot_e_contrast_volcano(ax):
    f2.container(ax)
    v = pd.read_csv(DATA / "gsk_vs_aza_probe_volcano.csv")
    c = pd.read_csv(DATA / "gsk_vs_aza_counts.csv").iloc[0]
    a = _axes(ax, (0.15, 0.14, 0.78, 0.70))
    nlp = -np.log10(v["adj.P.Val"].clip(lower=1e-300))
    sig = v["adj.P.Val"] < 0.05
    a.scatter(v.loc[~sig, "slope_diff_beta_0to200"], nlp[~sig], s=1, color="#c9c8c3", linewidths=0, rasterized=True,
              zorder=2)
    a.scatter(v.loc[sig, "slope_diff_beta_0to200"], nlp[sig], s=1.5, color=DRUG_COLORS["GSK"], linewidths=0,
              rasterized=True, zorder=3)
    a.axhline(-np.log10(0.05), ls=(0, (3, 2)), lw=0.6, color=f2.INK2, zorder=1)
    a.axvline(0, lw=0.4, color="#d6d5d0", zorder=0)
    lim = np.ceil(v["slope_diff_beta_0to200"].abs().quantile(0.9999) * 10) / 10
    a.set_xlim(-lim, lim)
    a.set_ylim(0, nlp.max() * 1.3)  # headroom above the points for the count labels
    a.set_xlabel("Δβ slope difference, GSK − AZA (0→200 nM)", fontsize=fs(6))
    a.set_ylabel("−log10 FDR", fontsize=fs(6))
    a.text(0.02, 0.97, f"FDR < 0.05: {int(c.n_fdr05_gsk_lower):,} probes\nGSK slope lower",
           transform=a.transAxes, ha="left", va="top", fontsize=fs(5.5), color=f2.INK2, linespacing=1.3)
    a.text(0.98, 0.97, f"FDR < 0.05: {int(c.n_fdr05_gsk_higher):,} probes\nGSK slope higher",
           transform=a.transAxes, ha="right", va="top", fontsize=fs(5.5), color=f2.INK2, linespacing=1.3)
    f2.panel_title(ax, "GSK vs AZA dose slope (~ donor + dose:drug)", x=0.5)
    return ax


def plot_f_dmr_volcano(ax):
    f2.container(ax)
    d = pd.read_csv(DATA / "gsk_vs_aza_dmr_volcano.csv")
    a = _axes(ax, (0.15, 0.14, 0.78, 0.70))
    nlp = -np.log10(d["min_smoothed_fdr"].clip(lower=1e-300))
    f = d["passes_filter"].astype(bool)
    a.scatter(d.loc[~f, "meandiff"], nlp[~f], s=3, color="#c9c8c3", linewidths=0, zorder=2)
    a.scatter(d.loc[f, "meandiff"], nlp[f], s=4, color=DRUG_COLORS["GSK"], edgecolors="#2e2d2a", linewidths=0.15,
              zorder=3)
    for xv in (-0.10, 0.10):
        a.axvline(xv, ls=(0, (3, 2)), lw=0.6, color=f2.INK2, zorder=1)
    a.axvline(0, lw=0.4, color="#d6d5d0", zorder=0)
    lim = max(0.3, np.ceil(d["meandiff"].abs().max() * 10) / 10)
    a.set_xlim(-lim, lim)
    a.set_ylim(0, nlp.max() * 1.35)  # headroom for the count labels
    a.set_xlabel("Δβ slope difference, GSK − AZA (0→200 nM)", fontsize=fs(6))
    a.set_ylabel("−log10 min smoothed FDR", fontsize=fs(6))
    n_lo, n_hi = int((f & (d["meandiff"] < 0)).sum()), int((f & (d["meandiff"] > 0)).sum())
    a.text(0.02, 0.97, f"|Δβ| > 0.10: {n_lo} DMRs\nGSK slope lower", transform=a.transAxes, ha="left", va="top",
           fontsize=fs(5.5), color=f2.INK2, linespacing=1.3,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=4)
    a.text(0.98, 0.97, f"|Δβ| > 0.10: {n_hi} DMRs\nGSK slope higher", transform=a.transAxes,
           ha="right", va="top", fontsize=fs(5.5), color=f2.INK2, linespacing=1.3,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=4)
    # wrapped so it stays left of the -0.10 line
    a.text(0.02, 0.80, f"{len(d):,} DMRs called\n{n_hi + n_lo} with |Δβ| > 0.10\n"
           f"{int((f & d['overlaps_507']).sum())} overlap the\n507 AZA DMRs",
           transform=a.transAxes, ha="left", va="top", fontsize=fs(5.5), color=f2.INK, linespacing=1.3,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=4)
    f2.panel_title(ax, "GSK vs AZA dose slope, DMRs", x=0.5)
    return ax


def _slope_scatter(ax, x, y, highlight, title, unit_label, note, rasterize=False, size=(2, 4)):
    """AZA slope (x) vs GSK slope (y), equal axes, y = x; `highlight` points in the GSK colour, the rest grey."""
    f2.container(ax)
    a = _axes(ax, (0.17, 0.14, 0.66, 0.70))
    lo = np.floor(min(x.min(), y.min()) / 0.05) * 0.05
    hi = max(0.05, np.ceil(max(x.max(), y.max()) / 0.05) * 0.05)
    a.plot([lo, hi], [lo, hi], ls=(0, (3, 2)), lw=0.6, color=f2.INK2, zorder=1)
    a.axhline(0, lw=0.4, color="#d6d5d0", zorder=0)
    a.axvline(0, lw=0.4, color="#d6d5d0", zorder=0)
    a.scatter(x[~highlight], y[~highlight], s=size[0], color="#c9c8c3", linewidths=0, rasterized=rasterize, zorder=2)
    a.scatter(x[highlight], y[highlight], s=size[1], color=DRUG_COLORS["GSK"], edgecolors="#2e2d2a",
              linewidths=0.1 if rasterize else 0.15, alpha=0.85, rasterized=rasterize, zorder=3)
    a.set_xlim(lo, hi); a.set_ylim(lo, hi); a.set_aspect("equal")
    a.set_xlabel(f"AZA slope, donors 1, 6, 7 ({unit_label})", fontsize=fs(6))
    a.set_ylabel(f"GSK slope, donor {DONOR_LABELS['EV03']} ({unit_label})", fontsize=fs(6))
    above, below = int((y > x).sum()), int((y < x).sum())
    a.text(0.03, 0.97, f"Above y = x (GSK higher): {above:,}\nBelow: {below:,}", transform=a.transAxes, ha="left",
           va="top", fontsize=fs(5.5), color=f2.INK2, linespacing=1.3,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=4)
    a.text(0.97, 0.03, note + f"\nMedian slope: AZA {np.median(x):.3f}, GSK {np.median(y):+.3f}",
           transform=a.transAxes, ha="right", va="bottom", fontsize=fs(5.5), color=f2.INK2, linespacing=1.3,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=4)
    a.text(hi - 0.01 * (hi - lo), hi - 0.12 * (hi - lo), "y = x", ha="right", va="top", fontsize=fs(5),
           color=f2.INK2)
    f2.panel_title(ax, title, x=0.5)
    return ax


def plot_g_sig_probe_slopes(ax):
    d = pd.read_csv(DATA / "gsk_vs_aza_sig_probe_slopes.csv")
    return _slope_scatter(ax, d["slope_AZA"].to_numpy(), d["slope_GSK"].to_numpy(), np.ones(len(d), bool),
                          f"Dose slopes, FDR < 0.05 probes (n = {len(d):,})", "Δβ per 100 nM",
                          "Model slopes (~ donor + dose:drug)", rasterize=True, size=(1.5, 1.5))


def plot_h_dmr_slopes(ax):
    d = pd.read_csv(DATA / "gsk_vs_aza_dmr_slopes.csv")
    f = d["passes_filter"].astype(bool).to_numpy()
    return _slope_scatter(ax, d["slope_AZA"].to_numpy(), d["slope_GSK"].to_numpy(), f,
                          f"Dose slopes, DMRs (n = {len(d):,})", "Δβ per 100 nM",
                          f"Highlighted: |Δβ| > 0.10 (n = {int(f.sum())})")


PANELS = {"A": ("FigS_gskA_top500_variable", plot_a_top500),
          "B": ("FigS_gskB_aza_dmrs_delta", plot_b_dmr507_delta),
          "C": ("FigS_gskC_top500_delta_aza_gsk", plot_c_top500_delta),
          "D": ("FigS_gskD_dmr_slope_scatter", plot_d_slope_scatter),
          "E": ("FigS_gskE_slope_contrast_volcano", plot_e_contrast_volcano),
          "F": ("FigS_gskF_slope_contrast_dmr_volcano", plot_f_dmr_volcano),
          "G": ("FigS_gskG_sig_probe_slopes", plot_g_sig_probe_slopes),
          "H": ("FigS_gskH_dmr_slopes", plot_h_dmr_slopes)}
SLOT = (3.6, 3.4)


def composite():
    n_rows = (len(PANELS) + 1) // 2
    fig = plt.figure(figsize=(2 * SLOT[0], n_rows * SLOT[1]))
    for i, (k, (_, fn)) in enumerate(PANELS.items()):
        r, c = divmod(i, 2)
        ax = fig.add_axes([c / 2, (n_rows - 1 - r) / n_rows, 0.5, 0.94 / n_rows])
        fn(ax)
        f2.panel_letter(ax, k)
    return fig


if __name__ == "__main__":
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    fig = composite()
    for ext in ("pdf", "svg"):
        fig.savefig(OUT / f"{STEM}.{ext}")
    fig.savefig(OUT / f"{STEM}_preview.png", dpi=200)
    plt.close(fig)
    for k, (name, fn) in PANELS.items():
        fig = plt.figure(figsize=SLOT)
        fn(fig.add_axes([0, 0, 1, 0.94]))
        for ext in ("pdf", "svg"):
            fig.savefig(PANEL_DIR / f"{name}.{ext}")
        plt.close(fig)
    print("Saved", OUT / f"{STEM}.{{pdf,svg}}", "and panels in", PANEL_DIR)
