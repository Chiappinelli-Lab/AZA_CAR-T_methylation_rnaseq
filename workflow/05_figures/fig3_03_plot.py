#!/usr/bin/env python3
"""Figure 3 cytokine panel for the collaborators (ask of 2026-09-24): IL-2, IFN-g, TNF-a across 0/50/100/200 nM
AZA -- rows: DNA methylation in UNSTIMULATED cells (array, 3 donors), RNA expression in UNSTIMULATED cells
(RNA-seq, 5 donors; added at the user's request) and in STIMULATED cells (RNA-seq, 3 donors).

Run in env `data_analysis` after fig3_01_prep_rna.R (deseq2) and fig3_02_prep_meth.R (haydar_methylation).

Design choices: (1) values are each donor's CHANGE from its own 0 nM sample (delta beta; log2 fold change) --
donor baselines differ far more than the dose effect (e.g. Il2 beta 0.56 in donor 7 vs 0.70 in donor 1), which
otherwise hides the trend and puts the mean between donors; (2) with n = 3 per dose a box plot is just three
points drawn as a box, so every donor is shown as its own points joined across doses (dose is within-donor).
Mean bars and a "low expression" note on Il2 were removed at the user's request (the low IL-2 RNA is covered in
the legend notes).
`--absolute` draws the alternative version requested 2026-09-24: absolute beta and DESeq2-normalized counts
instead of change from 0 nM (-> Fig3_cytokines_meth_rna_absolute.{pdf,svg}).
`--bars` (2026-09-25, the joined lines were hard to follow across doses): one bar per dose = mean over
donors +/- SEM, bar shaded by the main figure's dose greys, each donor as its own colored point on the bar, no
lines (-> Fig3_cytokines_meth_rna_{FC,absolute}_bars.{pdf,svg,png}; combines with --absolute).
`--ifn` (2026-09-25, collaborators asked whether IL-2 is as low as IFN-a/b): adds IFN-b (Ifnb1) and IFN-a (Ifna family
pooled) columns from the *_ifn.csv prep files; they miss the RNA-seq expression filter, so their RNA panels show
mean counts instead of a model estimate (-> ..._ifn.png only, for now).
Style (fonts, donor colors and 1-7 numbering, +1 pt font size of the v1 figure) comes from fig2_03_plot.py.
Output: paper_figures/fig3_cytokines/Fig3_cytokines_meth_rna_FC.{pdf,svg,png} (change from 0 nM; "_FC" added
2026-09-25, was unsuffixed)
"""
import importlib.util
import sys
from pathlib import Path

if "--font-bump" not in sys.argv:          # match the v1 main figure (+1 pt) unless told otherwise
    sys.argv += ["--font-bump", "1"]
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fig2", HERE / "fig2_03_plot.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = F.FIG_DIR / "fig3_cytokines"
DATA = OUT / "data"
GENES = [("Il2", "IL-2"), ("Ifng", "IFN-γ"), ("Tnf", "TNF-α")]
IFN = "--ifn" in sys.argv
if IFN:
    GENES += [("Ifnb1", "IFN-β"), ("Ifna", "IFN-α")]
DOSES = F.DOSES
fs, INK, INK2, MUTED = F.fs, F.INK, F.INK2, F.MUTED
ABSOLUTE = "--absolute" in sys.argv
BARS = "--bars" in sys.argv
SUFFIX = ("_absolute" if ABSOLUTE else "_FC") + ("_bars" if BARS else "") + ("_ifn" if IFN else "")   # _FC = change from 0 nM


def change_from_zero(df, value):
    """Each donor's value minus its own 0 nM value."""
    base = df[df.dose_nM == 0].set_index(["symbol", "donor"])[value]
    return df[value] - [base[(g, d)] for g, d in zip(df.symbol, df.donor)]


def signed(x, nd):
    x = round(x, nd) + 0.0          # avoid "-0.00"
    return f"{x:+.{nd}f}".replace("-", "−")


def fmt_p(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"


def donor_series(ax, df, value, donors):
    """One colored point per sample, joined across doses per donor."""
    offsets = dict(zip(donors, np.linspace(-0.12, 0.12, len(donors))))
    for d in donors:
        sub = df[df.donor == d].sort_values("dose_nM")
        x = sub.dose_nM.map(DOSES.index) + offsets[d]
        ax.plot(x, sub[value], color=F.DONOR_COLORS[d], lw=0.8, alpha=0.6, zorder=2)
        ax.scatter(x, sub[value], s=16, color=F.DONOR_COLORS[d], edgecolors="white", linewidths=0.5, zorder=3)
    if not ABSOLUTE:
        ax.axhline(0, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_xticks(range(len(DOSES)))
    ax.set_xticklabels([str(d) for d in DOSES])
    ax.set_xlim(-0.5, len(DOSES) - 0.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="both", labelsize=fs(5.5), length=2, pad=1.5)
    ax.margins(y=0.15)


def donor_bars(ax, df, value, donors):
    """Mean +/- SEM bar per dose (over the donors that have that dose), one colored point per donor on top."""
    offsets = dict(zip(donors, np.linspace(-0.22, 0.22, len(donors))))
    for k, dose in enumerate(DOSES):
        v = df.loc[df.dose_nM == dose, value]
        sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0
        ax.bar(k, v.mean(), width=0.7, color=F.DOSE_COLORS[dose], alpha=0.45, edgecolor=INK2, linewidth=0.5,
               zorder=1)
        ax.errorbar(k, v.mean(), yerr=sem, color=INK, elinewidth=0.6, capsize=2, capthick=0.6, zorder=2)
    for d in donors:
        sub = df[df.donor == d]
        ax.scatter(sub.dose_nM.map(DOSES.index) + offsets[d], sub[value], s=14, color=F.DONOR_COLORS[d],
                   edgecolors="white", linewidths=0.5, zorder=3)
    if not ABSOLUTE:   # absolute bars start at the x-axis already
        ax.axhline(0, color=INK2, lw=0.5, zorder=1.5)
    ax.set_xticks(range(len(DOSES)))
    ax.set_xticklabels([str(d) for d in DOSES])
    ax.set_xlim(-0.6, len(DOSES) - 0.4)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="both", labelsize=fs(5.5), length=2, pad=1.5)
    ax.margins(y=0.15)


def stat_label(ax, text):
    ax.text(0.5, 1.0, text, transform=ax.transAxes, fontsize=fs(5), color=INK2, ha="center", va="bottom")


def read(name):
    """A prep CSV, plus its type I IFN counterpart (<name>_ifn.csv) under --ifn."""
    df = pd.read_csv(DATA / f"{name}.csv")
    return pd.concat([df, pd.read_csv(DATA / f"{name}_ifn.csv")], ignore_index=True) if IFN else df


def load_rows():
    """(row label, long data, value column, stats table, stat formatter, donors, y label) per row."""
    meth = read("meth_unstimulated_long")
    meth["value"] = meth.beta if ABSOLUTE else change_from_zero(meth, "beta")
    rows = [("DNA methylation\nunstimulated", meth, read("meth_unstimulated_stats"),
             lambda st: f"Δβ {signed(st.delta_beta_0to200, 3)}, {fmt_p(st.pvalue)}",
             ["EV07", "EV08", "ZY03"], "Methylation (β)" if ABSOLUTE else "Δβ vs 0 nM")]
    for cond, donors in (("unstimulated", ["EV07", "EV17", "EV18", "EV19", "EV20"]),
                         ("stimulated", ["EV07", "EV17", "EV20"])):
        rna = read(f"rna_{cond}_long")
        rna["log2_norm"] = np.log2(rna.norm_count + 1)
        rna["value"] = rna.norm_count if ABSOLUTE else change_from_zero(rna, "log2_norm")
        rows.append((f"RNA expression\n{cond}", rna, pd.read_csv(DATA / f"rna_{cond}_stats.csv"),
                     lambda st: f"LFC {signed(st.log2FC_per100nM, 2)}, {fmt_p(st.pvalue)}",
                     donors, "Normalized counts" if ABSOLUTE else "log2 fold change\nvs 0 nM"))
    return rows


def plot_cytokine_panel(fig):
    rows = load_rows()
    # absolute version: shared y-axes (tick labels only on the left) allow much tighter spacing
    hspace, wspace = (0.42, 0.12) if ABSOLUTE else (0.62, 0.45)
    W = fig.get_figwidth()   # margins in inches (0.85 left, 0.75 right) so extra gene columns only widen the grid
    gs = fig.add_gridspec(len(rows), len(GENES), left=0.85 / W, right=1 - 0.75 / W, top=0.91, bottom=0.08, hspace=hspace, wspace=wspace)
    axes = {}
    for i, (row_label, df, stats, fmt, donors, ylabel) in enumerate(rows):
        stats = stats.set_index("symbol")
        for j, (sym, protein) in enumerate(GENES):
            ax = axes[i, j] = fig.add_subplot(gs[i, j])
            (donor_bars if BARS else donor_series)(ax, df[df.symbol == sym], "value", donors)
            if sym in stats.index:
                stat_label(ax, fmt(stats.loc[sym]))
            else:   # type I IFNs: below the RNA-seq filter, no dose model -- report how little there is
                stat_label(ax, f"not tested, mean {df.loc[df.symbol == sym, 'norm_count'].mean():.1f} counts")
            if i == 0:
                ax.set_title(f"{protein} ($\\it{{{sym}}}$)", fontsize=fs(7), fontweight="bold", pad=11)
            if i == len(rows) - 1:
                ax.set_xlabel("AZA (nM)", fontsize=fs(6), labelpad=1.5)
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=fs(6))
                box = ax.get_position()
                fig.text(0.015, (box.y0 + box.y1) / 2, row_label, rotation=90, fontsize=fs(6.5), fontweight="bold",
                         color=INK, ha="left", va="center", multialignment="center", linespacing=1.2)
    # one shared donor legend (all donors in any row, in number order), top-aligned with the first row
    all_donors = sorted({d for row in rows for d in row[4]}, key=lambda d: int(F.DONOR_LABELS[d]))
    handles = [Line2D([], [], marker="o", color=F.DONOR_COLORS[d], markeredgecolor="white", markersize=4.5, lw=0 if BARS else 0.8,
                      label=F.DONOR_LABELS[d]) for d in all_donors]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(1 - 0.65 / W, axes[0, 0].get_position().y1),
               frameon=False, title="Donor", title_fontsize=fs(5.5), fontsize=fs(5.5), handlelength=1.4,
               labelspacing=0.25, borderaxespad=0, alignment="left")
    if ABSOLUTE:
        # Shared y-axes from 0 so levels compare across genes (user request 2026-09-24): one for the methylation
        # row, one across both RNA rows (counts use joint size factors, so stimulated vs unstimulated compare too).
        for group in ([(0, j) for j in range(len(GENES))],
                      [(i, j) for i in (1, 2) for j in range(len(GENES))]):
            top = max(axes[k].get_ylim()[1] for k in group)
            for k in group:
                axes[k].set_ylim(-0.03 * top, top)   # small margin so points at 0 aren't hidden by the axis
                if k[1] > 0:
                    axes[k].tick_params(labelleft=False)


if __name__ == "__main__":
    fig = plt.figure(figsize=(5.0 + 1.25 * (len(GENES) - 3), 5.0))
    plot_cytokine_panel(fig)
    for ext in () if IFN else ("pdf", "svg"):
        fig.savefig(OUT / f"Fig3_cytokines_meth_rna{SUFFIX}.{ext}")
    fig.savefig(OUT / f"Fig3_cytokines_meth_rna{SUFFIX}.png", dpi=600, facecolor="white")   # raster copy
    if "--preview" in sys.argv:
        fig.savefig(sys.argv[sys.argv.index("--preview") + 1], dpi=200)
    print("Saved", OUT / f"Fig3_cytokines_meth_rna{SUFFIX}.{{{'png' if IFN else 'pdf,svg,png'}}}")
