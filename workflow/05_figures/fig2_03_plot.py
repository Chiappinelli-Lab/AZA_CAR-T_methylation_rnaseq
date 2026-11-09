#!/usr/bin/env python3
"""Paper figure: AZA dose-response in B7H3 CAR-T cells, RNA-seq x DNA methylation.

Panels (one plot_* function each, all drawing into a single `ax` container, so
the same code builds the composite and the standalone panels -- the composite is
NOT stitched from saved files):
  A  plot_a_de_heatmap         111 DE genes, donor-corrected VST, row z-score
  B  plot_b_dmr_heatmap        507 DMRs, donor-adjusted region beta, row z-score
  C  plot_c_gene_venn          DE genes vs genes with a promoter/body probe in a DMR
  D  plot_d_pathway_venn       GO:BP terms: RNA-seq GSEA vs methylation dmGsea
  E  plot_e_shared_pathways    the 26 shared terms: Wang semantic-similarity tree + term labels
                               (_pathway_tree) and RNA vs methylation p / gene-count dotplot
                               (_pathway_dotplot), one panel sharing rows
  F  plot_f_geneset_heatmap    one candidate gene set (default effector function;
                               `--set <name from candidate_gene_sets.R>` to swap)
Each standalone panel is drawn at the same size as its slot in the composite, so
they are identical apart from the panel letter.

Outputs (vector; fonts embedded as TrueType so text stays editable):
  figures/fig2_rna_methylation/Fig2_rna_methylation.{pdf,svg}
  figures/fig2_rna_methylation/panels/Fig2<letter>_<name>.{pdf,svg}

"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.patches import Rectangle, Patch
from matplotlib_venn import venn2
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "py"))
import paths  # noqa: E402  (py/paths.py: repo-relative inputs/outputs)
FIG_DIR = paths.FIG_ROOT                              # shared with fig3_03 / figS_gsk_03: <repo>/figures
FIG2_DIR = FIG_DIR / "fig2_rna_methylation"
DATA = FIG2_DIR / "data"
PANEL_DIR = FIG2_DIR / "panels"
STEM = "Fig2_rna_methylation"
RNA_DIR = paths.RNA_NS_DIR
DMG_DIR = paths.DMGSEA_DIR
INT_RES = paths.INT_DIR

GENESET = sys.argv[sys.argv.index("--set") + 1] if "--set" in sys.argv else "01_Effector_Function"
FDR_CUT = 0.05
# --font-bump N adds N pt to every font size (for trying sizes); such runs write
# Fig2_rna_methylation_font+N.{pdf,svg} only and leave the standard outputs alone.
FONT_BUMP = float(sys.argv[sys.argv.index("--font-bump") + 1]) if "--font-bump" in sys.argv else 0.0


# --tag NAME writes a frozen, shareable copy: fig2_rna_methylation/release_NAME/Fig2_rna_methylation_NAME.{pdf,svg}
# and .../release_NAME/panels/Fig2<panel>_NAME.{pdf,svg} (composite + all standalone panels, at the current --font-bump).
TAG = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else None


def fs(size):
    return size + FONT_BUMP

# ---------------------------------------------------------------- style ----
for f in paths.FONT_DIR.glob("*.ttf"):
    font_manager.fontManager.addfont(str(f))
# The paper figures use Arial 2.82 (Microsoft core web fonts, conda-forge `mscorefonts`). Arial isn't shipped in
# this repo (license); drop arial*.ttf into fonts/ to get it. Otherwise Liberation Sans (metric-compatible, shipped)
# is used, so the layout is the same and only the glyph shapes differ.
FONT = "Arial" if any(f.name == "Arial" for f in font_manager.fontManager.ttflist) else "Liberation Sans"
mpl.rcParams.update({
    "font.family": [FONT],
    "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic", "mathtext.bf": f"{FONT}:bold",
    "font.size": fs(6), "axes.titlesize": fs(7), "axes.labelsize": fs(6),
    "xtick.labelsize": fs(5.5), "ytick.labelsize": fs(5.5), "legend.fontsize": fs(5.5),
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.titleweight": "bold", "savefig.dpi": 300,
})
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
NA_COLOR = "#e4e3df"
# Diverging blue <-> red, grey midpoint (same poles as the earlier project heatmaps).
Z_CMAP = LinearSegmentedColormap.from_list("z", ["#2166ac", "#6ea3cf", "#f0efec", "#e0806e", "#b2182b"])
SEQ_CMAP = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
# Donors: EV07 is the one donor in both cohorts and keeps one color; every other
# donor gets its own hue (validated with dataviz validate_palette.js: RNA five
# adjacent-pairs pass; methylation three all-pairs pass).
DONOR_COLORS = {"EV07": "#e87ba4", "EV17": "#2a78d6", "EV18": "#008300", "EV19": "#eda100",
                "EV20": "#4a3aa7", "EV08": "#a8541a", "ZY03": "#16866d"}
# Figure labels: one number per person across both cohorts (confirmed with the user 2026-09-24), so the shared
# donor EV07 is "1" in both the RNA-seq and methylation panels.
DONOR_LABELS = {"EV07": "1", "EV17": "2", "EV18": "3", "EV19": "4", "EV20": "5", "EV08": "6", "ZY03": "7"}
DOSES = [0, 50, 100, 200]
DOSE_COLORS = dict(zip(DOSES, ["#f0efec", "#b9b8b2", "#7f7e79", "#2e2d2a"]))
OMIC_COLORS = {"rna": "#2a78d6", "meth": "#eb6834"}


# --------------------------------------------------------------- helpers ----
def container(ax):
    """Turn a panel's ax into an invisible container for inset axes."""
    ax.set_axis_off()
    return ax


def inset(ax, x0, y0, w, h):
    a = ax.inset_axes([x0, y0, w, h])
    a._inset_parent, a._inset_h = ax, h
    return a


def row_z(m):
    m = np.asarray(m, float)
    mu = np.nanmean(m, axis=1, keepdims=True)
    sd = np.nanstd(m, axis=1, ddof=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (m - mu) / sd
    z[~np.isfinite(z)] = np.nan
    return z


def order_samples(meta, donor_order):
    meta = meta.copy()
    meta["donor_rank"] = meta["donor"].map({d: i for i, d in enumerate(donor_order)})
    return meta.sort_values(["dose_nM", "donor_rank"]).reset_index(drop=True)


def draw_matrix(ax, z, cmap=Z_CMAP, vmin=-2, vmax=2):
    """Vector cells (pcolormesh quads), row 0 at the top."""
    n_r, n_c = z.shape
    zm = np.ma.masked_invalid(z)
    cm = cmap.copy()
    cm.set_bad(NA_COLOR)
    ax.pcolormesh(np.arange(n_c + 1) - 0.5, np.arange(n_r + 1) - 0.5, zm, cmap=cm, vmin=vmin, vmax=vmax,
                  edgecolors="face", linewidth=0.05, antialiased=False)
    ax.set_xlim(-0.5, n_c - 0.5)
    ax.set_ylim(n_r - 0.5, -0.5)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def dose_block_gaps(ax, meta, n_rows, lw=0.8):
    for i in np.where(np.diff(meta["dose_nM"].to_numpy()) != 0)[0]:
        ax.plot([i + 0.5, i + 0.5], [-0.5, n_rows - 0.5], color="white", lw=lw, solid_capstyle="butt")


def dose_block_labels(ax, meta, y=-0.04, unit=True):
    for d in DOSES:
        idx = np.where(meta["dose_nM"].to_numpy() == d)[0]
        if len(idx):
            ax.text(idx.mean(), y, f"{d}", transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=fs(5.5), color=INK2)
    if unit:
        ax.text(1.0, y - 0.07, "AZA (nM)", transform=ax.transAxes, ha="right", va="top", fontsize=fs(5.5), color=INK2)


def annotation_strip(ax, meta, key, colors, label=None, label_side="right"):
    for i, v in enumerate(meta[key]):
        ax.add_patch(Rectangle((i - 0.5, 0), 1, 1, facecolor=colors[v], edgecolor="none"))
    for i in np.where(np.diff(meta["dose_nM"].to_numpy()) != 0)[0]:
        ax.plot([i + 0.5, i + 0.5], [0, 1], color="white", lw=0.8, solid_capstyle="butt")
    ax.set_xlim(-0.5, len(meta) - 0.5); ax.set_ylim(0, 1)
    ax.set_axis_off()
    if label:
        x, ha = (1.02, "left") if label_side == "right" else (-0.02, "right")
        ax.text(x, 0.5, label, transform=ax.transAxes, ha=ha, va="center", fontsize=fs(5.5), color=INK2)


def draw_dendrogram(ax, Z, n_leaves, orientation="left", color=INK2, lw=0.4, link_color_func=None):
    """Draw a scipy linkage with leaf i at row i (row 0 on top, matching draw_matrix).
    link_color_func(link_id) -> color colors individual links (default: all `color`)."""
    d = dendrogram(Z, no_plot=True, orientation=orientation,
                   link_color_func=link_color_func or (lambda k: color))
    hmax = max(max(x) for x in d["dcoord"])
    for ic, dc, col in zip(d["icoord"], d["dcoord"], d["color_list"]):
        rows = (np.array(ic) - 5) / 10
        heights = np.array(dc) / hmax
        if orientation == "left":        # root on the left, leaves on the right
            ax.plot(1 - heights, rows, color=col, lw=lw, solid_capstyle="butt")
        else:
            ax.plot(heights, rows, color=col, lw=lw, solid_capstyle="butt")
    ax.set_ylim(n_leaves - 0.5, -0.5)
    ax.set_xlim(-0.02, 1.0)
    ax.set_axis_off()
    return d["leaves"]


def z_colorbar(ax, label="Row z-score", cmap=Z_CMAP, vmin=-2, vmax=2, ticks=(-2, 0, 2), orientation="vertical"):
    cb = plt.colorbar(mpl.cm.ScalarMappable(Normalize(vmin, vmax), cmap), cax=ax, orientation=orientation,
                      ticks=list(ticks))
    cb.solids.set_rasterized(False)   # matplotlib rasterizes colorbar gradients by default; keep vector
    cb.solids.set_edgecolor("face")   # ...and without this the vector segments show hairline seams
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=fs(5), width=0.4, length=1.5, pad=1)
    if orientation == "vertical":
        cb.ax.set_title(label, fontsize=fs(5.5), pad=3, loc="left", fontweight="normal")
    else:
        cb.set_label(label, fontsize=fs(5.5), labelpad=1)
    return cb


def donor_legend(ax, donors, ncol=1):
    ax.set_axis_off()
    ax.legend(handles=[Patch(facecolor=DONOR_COLORS[d], edgecolor="none", label=DONOR_LABELS[d]) for d in donors],
              loc="upper left", frameon=False, title="Donor", title_fontsize=fs(5.5), fontsize=fs(5.5), ncol=ncol,
              handlelength=0.9, handleheight=0.9, handletextpad=0.4, labelspacing=0.3, columnspacing=0.8,
              borderaxespad=0, borderpad=0, alignment="left")


def sample_legends(ax, donors):
    """Donor legend with the AZA dose legend stacked underneath (both in one axes, top-left anchored)."""
    ax.set_axis_off()
    kw = dict(frameon=False, title_fontsize=fs(5.5), fontsize=fs(5.5), handlelength=0.9, handleheight=0.9,
              handletextpad=0.4, labelspacing=0.3, borderaxespad=0, borderpad=0, alignment="left")
    leg1 = ax.legend(handles=[Patch(facecolor=DONOR_COLORS[d], edgecolor="none", label=DONOR_LABELS[d]) for d in donors],
                     loc="upper left", title="Donor", **kw)
    ax.add_artist(leg1)
    h_pts = container_height_pts(ax)
    drop = 1.33 * fs(5.5) * (len(donors) + 1) + 7   # legend row pitch ~1.33 x font size, plus a gap
    ax.legend(handles=[Patch(facecolor=DOSE_COLORS[d], edgecolor=INK2 if d == 0 else "none", linewidth=0.3,
                             label=str(d)) for d in DOSES],
              loc="upper left", bbox_to_anchor=(0, 1 - drop / h_pts), title="AZA (nM)", **kw)


def container_height_pts(inset_ax):
    """Height in points of an inset axes (its parent container has a fixed figure position)."""
    pos = inset_ax._inset_parent.get_position()
    return inset_ax._inset_h * pos.height * inset_ax.figure.get_figheight() * 72


LETTER_PT = 11   # titles are indented this many points so the composite letter fits in front


def panel_title(ax, title, subtitle=None, center=True, x=None):
    """Centered on the panel by default; `x` (axes fraction) centers it elsewhere, or left-anchors it there
    when center=False."""
    if x is None:
        x, dx = (0.5, 0) if center else (0, LETTER_PT)
    else:
        dx = 0
    ax.annotate(title, xy=(x, 1), xycoords="axes fraction", xytext=(dx, 0), textcoords="offset points",
                fontsize=fs(7), fontweight="bold", va="top", ha="center" if center else "left", color=INK)
    if subtitle:
        ax.annotate(subtitle, xy=(0, 1), xycoords="axes fraction", xytext=(LETTER_PT, -1.35 * fs(7)),
                    textcoords="offset points", fontsize=fs(5.5), va="top", ha="left", color=INK2, linespacing=1.25)


def panel_letter(ax, letter):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(0, 1), textcoords="offset points",
                fontsize=fs(9), fontweight="bold", va="top", ha="left", color=INK)


def cluster_rows(z):
    filled = np.where(np.isnan(z), 0, z)
    return linkage(filled, method="ward", metric="euclidean")


# ---------------------------------------------------------- data loaders ----
def load_heatmap(prefix, id_col, drop_cols):
    mat = pd.read_csv(DATA / f"{prefix}_matrix.csv")
    meta = pd.read_csv(DATA / f"{prefix}_samples.csv")
    return mat, meta


def load_shared_pathways():
    shared = pd.read_csv(DMG_DIR / "meth_rna_shared_GO_BP.csv")
    rna = pd.read_csv(RNA_DIR / "enrichment_gsea/GO_BP_gsea_full.csv").set_index("ID")
    meth = pd.read_csv(DMG_DIR / "dmgsea_GO_BP.csv").set_index("term_id")
    counts = pd.read_csv(DATA / "pathway_rna_counts.csv").set_index("term_id")
    shared["rna_p"] = rna.loc[shared.term_id, "pvalue"].to_numpy()
    shared["meth_p"] = meth.loc[shared.term_id, "p"].to_numpy()
    shared["rna_n_DE"] = counts.loc[shared.term_id, "rna_n_DE"].to_numpy()
    leaves = pd.read_csv(DATA / "pathway_tree_leaves.csv").set_index("term_id")
    shared["cluster"] = leaves.loc[shared.term_id, "cluster"].to_numpy()
    return shared.set_index("term_id")


def pathway_tree():
    """R hclust (ward.D on 1 - Wang similarity) -> scipy linkage; returns (Z, term_ids in hclust label order)."""
    sim = pd.read_csv(DATA / "pathway_termsim_wang.csv")
    labels = sim["term_id"].tolist()
    hc = pd.read_csv(DATA / "pathway_tree_hclust.csv")
    n = len(labels)
    conv = lambda x: -x - 1 if x < 0 else n + x - 1
    size = {i: 1 for i in range(n)}
    Z = []
    for k, (a, b, h) in enumerate(hc[["merge1", "merge2", "height"]].itertuples(index=False)):
        a, b = conv(int(a)), conv(int(b))
        size[n + k] = size[a] + size[b]
        Z.append([a, b, h, size[n + k]])
    return np.array(Z, float), labels


# ------------------------------------------------------------ panel A/B ----
def _omics_heatmap(ax, z, meta, donor_order, title, z_label):
    container(ax)
    Z = cluster_rows(z)
    order = leaves_list(Z)
    z = z[order]
    # layout (fractions of the container)
    x_dend, w_dend = 0.0, 0.10
    x_hm, w_hm = 0.11, 0.62
    y_hm, h_hm = 0.12, 0.76
    strip_h, strip_gap = 0.028, 0.008
    ax_hm = inset(ax, x_hm, y_hm, w_hm, h_hm)
    draw_matrix(ax_hm, z)
    dose_block_gaps(ax_hm, meta, z.shape[0])
    ax_dend = inset(ax, x_dend, y_hm, w_dend, h_hm)
    draw_dendrogram(ax_dend, Z, z.shape[0])
    y = y_hm + h_hm + strip_gap
    annotation_strip(inset(ax, x_hm, y, w_hm, strip_h), meta, "donor", DONOR_COLORS, "Donor", label_side="left")
    annotation_strip(inset(ax, x_hm, y + strip_h + strip_gap, w_hm, strip_h), meta, "dose_nM", DOSE_COLORS, "Dose",
                     label_side="left")
    panel_title(ax, title, x=x_hm + w_hm / 2)   # centered on the heatmap, not heatmap + legends
    # right: donor + dose legends; bottom: horizontal z colorbar
    top = y + 2 * strip_h + strip_gap
    sample_legends(inset(ax, x_hm + w_hm + 0.04, y_hm, 0.24, top - y_hm), donor_order)
    z_colorbar(inset(ax, x_hm + w_hm * 0.2, 0.065, w_hm * 0.6, 0.022), label=z_label, orientation="horizontal")
    ax_hm.set_rasterization_zorder(None)
    return ax


def plot_a_de_heatmap(ax):
    mat = pd.read_csv(DATA / "de_heatmap_matrix.csv")
    meta = pd.read_csv(DATA / "de_heatmap_samples.csv")
    donors = ["EV07", "EV17", "EV18", "EV19", "EV20"]
    meta = order_samples(meta, donors)
    z = row_z(mat[meta["sample"]].to_numpy())
    n_up, n_dn = int((mat.log2FoldChange > 0).sum()), int((mat.log2FoldChange < 0).sum())
    return _omics_heatmap(ax, z, meta, donors, "Gene expression (RNA-seq)",
                          "Row z-score of donor-corrected VST expression")


def plot_b_dmr_heatmap(ax):
    mat = pd.read_csv(DATA / "dmr_heatmap_matrix.csv")
    meta = pd.read_csv(DATA / "dmr_heatmap_samples.csv")
    donors = ["EV07", "EV08", "ZY03"]
    meta = order_samples(meta, donors)
    z = row_z(mat[meta["sample"]].to_numpy())
    n_hypo = int((mat.meandiff < 0).sum())
    return _omics_heatmap(ax, z, meta, donors, "DNA methylation (Infinium MM285)",
                          "Row z-score of donor-adjusted DMR β")


# ------------------------------------------------------------ panel C/D ----
def _venn(ax, sets_counts, labels, title):
    """Area-proportional 2-set Venn. Counts of small regions move outside on a leader line."""
    container(ax)
    a_only, b_only, both = sets_counts
    total = a_only + b_only + both
    v_ax = inset(ax, 0.06, 0.16, 0.74, 0.64)
    v = venn2(subsets=(a_only, b_only, both), set_labels=labels, ax=v_ax,
              set_colors=(OMIC_COLORS["rna"], OMIC_COLORS["meth"]), alpha=0.4)
    for pid in ("10", "01", "11"):
        p = v.get_patch_by_id(pid)
        if p is not None:
            p.set_edgecolor("white"); p.set_linewidth(0.6)
    for t in v.set_labels:
        t.set_fontsize(fs(7)); t.set_color(INK)
    (cxa, cya), (cxb, cyb) = [(float(c.x), float(c.y)) for c in v.centers]
    ra, rb = v.radii
    bottom = min(cya - ra, cyb - rb)
    if rb < 0.5 * ra:   # small second circle: its set label goes above-right, clear of the leader lines
        v.set_labels[1].set_position((cxb + 0.6 * rb, cyb + rb + 0.06))
        v.set_labels[1].set_ha("left"); v.set_labels[1].set_va("bottom")
    for pid, n in (("10", a_only), ("01", b_only), ("11", both)):
        t = v.get_label_by_id(pid)
        if t is None:
            continue
        t.set_fontsize(fs(7.5)); t.set_color(INK)
        if n / total >= 0.06:
            continue
        x, y = t.get_position()
        if pid == "11":
            tx, ty, ha, va = x, bottom - 0.12, "center", "top"
        else:
            tx, ty, ha, va = cxb + rb + 0.1, y, "left", "center"
        t.set_visible(False)
        v_ax.annotate(str(n), xy=(x, y), xytext=(tx, ty), fontsize=fs(7.5), color=INK, ha=ha, va=va,
                      arrowprops=dict(arrowstyle="-", lw=0.4, color=INK2, shrinkA=1, shrinkB=0))
    v_ax.set_axis_off()
    panel_title(ax, title)
    return ax


def plot_c_gene_venn(ax):
    df = pd.read_csv(INT_RES / "gene_level_multiomic_table.csv")
    de, dmr = df.DE_sig.astype(bool), df.is_sig_DMR_gene.astype(bool)
    both, de_only, dmr_only = int((de & dmr).sum()), int((de & ~dmr).sum()), int((~de & dmr).sum())
    neither = len(df) - both - de_only - dmr_only
    odds, p = fisher_exact([[both, de_only], [dmr_only, neither]], alternative="greater")
    print(f"[legend C] {len(df):,} genes measured on both platforms; {de_only + both} DE, {dmr_only + both} with a "
          f"probe in a DMR, {both} both; one-sided Fisher OR = {odds:.2f}, p = {p:.2g}")
    return _venn(ax, (de_only, dmr_only, both), ("DE genes", "Genes with a\nprobe in a DMR"),
                 "Gene overlap")


def plot_d_pathway_venn(ax):
    rna_all = pd.read_csv(RNA_DIR / "enrichment_gsea/GO_BP_gsea_all_terms.csv")
    meth = pd.read_csv(DMG_DIR / "dmgsea_GO_BP.csv")
    universe = set(rna_all.ID) & set(meth.term_id)
    rna_sig = set(rna_all.loc[rna_all["p.adjust"] < FDR_CUT, "ID"]) & universe
    meth_sig = set(meth.loc[meth.FDR < FDR_CUT, "term_id"]) & universe
    both = len(rna_sig & meth_sig)
    a, b = len(rna_sig) - both, len(meth_sig) - both
    odds, p = fisher_exact([[both, a], [b, len(universe) - both - a - b]], alternative="greater")
    print(f"[legend D] {len(universe):,} GO:BP terms tested in both; FDR < {FDR_CUT}: {a + both} RNA-seq, "
          f"{b + both} methylation, {both} both; one-sided Fisher OR = {odds:.1f}, p = {p:.2g}")
    return _venn(ax, (a, b, both), ("RNA-seq GSEA", "Methylation\ndmGsea"),
                 "GO:BP pathway overlap")


# -------------------------------------------------------------- panel E ----
# One hue per semantic cluster, assigned top-to-bottom in row order (validated categorical palette slots 1-5,
# whose adjacent pairs pass the CVD checks -- clusters are only ever adjacent to their row neighbours).
CLUSTER_PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]


def cluster_colors(clusters_in_row_order):
    seen = list(dict.fromkeys(clusters_in_row_order))
    return {c: CLUSTER_PALETTE[i] for i, c in enumerate(seen)}


def tint(color, alpha):
    """Opaque tint of `color` over white (opaque so overlapping marks don't change it)."""
    r, g, b = to_rgb(color)
    return (1 - alpha + alpha * r, 1 - alpha + alpha * g, 1 - alpha + alpha * b)


def _cluster_bands(ax, clusters_in_row_order, x0=0, x1=1, alpha=0.12):
    """Color block per semantic cluster (same rows and colors in E and F)."""
    cols = cluster_colors(clusters_in_row_order)
    start = 0
    for i in range(1, len(clusters_in_row_order) + 1):
        if i == len(clusters_in_row_order) or clusters_in_row_order[i] != clusters_in_row_order[start]:
            ax.add_patch(Rectangle((x0, start - 0.5), x1 - x0, i - start, facecolor=tint(cols[clusters_in_row_order[start]], alpha),
                                   edgecolor="none", zorder=0, clip_on=False))
            start = i


def shared_row_order():
    Z, labels = pathway_tree()
    d = dendrogram(Z, no_plot=True)
    return Z, [labels[i] for i in d["leaves"]]


# Cluster names: cleaned-up versions of the enrichplot treeplot word-frequency labels (12_dmgsea_plots.R),
# attached to each cluster via one anchor term it contains.
CLUSTER_LABELS = {
    "regulation of T cell activation": "Lymphocyte\nactivation and\nproliferation",
    "regulation of immune effector process": "Immune effector\nactivation",
    "mononuclear cell proliferation": "Leukocyte\nproliferation",
    "leukocyte cell-cell adhesion": "Leukocyte cell-cell\nadhesion",
    "regulation of hemopoiesis": "Hemopoiesis and\nleukocyte\ndifferentiation",
}


def _pathway_tree(ax, cluster_frac=0.22, tree_frac=0.19, fontsize=fs(6.5)):
    """[cluster name | bar][dendrogram][term labels]; fractions are of the panel width."""
    container(ax)
    sh = load_shared_pathways()
    Z, rows = shared_row_order()
    n = len(rows)
    y0, h = 0.02, 0.905
    x_tree = cluster_frac
    x_lab = cluster_frac + tree_frac
    clusters = sh.loc[rows, "cluster"].tolist()
    cols = cluster_colors(clusters)
    # a link takes its cluster's color when every leaf under it is in that cluster
    _, labels = pathway_tree()
    leaf_cluster = [sh.loc[t, "cluster"] for t in labels]
    members = {i: {leaf_cluster[i]} for i in range(len(labels))}
    for k, (a, b) in enumerate(Z[:, :2].astype(int)):
        members[len(labels) + k] = members[a] | members[b]
    link_color = lambda k: cols[next(iter(members[k]))] if len(members[k]) == 1 else INK2
    ax_t = inset(ax, x_tree, y0, tree_frac - 0.01, h)
    draw_dendrogram(ax_t, Z, n, lw=0.6, link_color_func=link_color)
    ax_l = inset(ax, 0.0, y0, 1.0, h)
    ax_l.set_xlim(0, 1); ax_l.set_ylim(n - 0.5, -0.5); ax_l.set_axis_off()
    x_bar = x_tree - 0.012
    _cluster_bands(ax_l, clusters, x0=x_bar + 0.006, x1=x_lab, alpha=0.28)
    _cluster_bands(ax_l, clusters, x0=x_lab, x1=1.01, alpha=0.12)   # overlaps F's blocks: no seam
    for i, t in enumerate(rows):
        ax_l.text(x_lab + 0.005, i, sh.loc[t, "term"], fontsize=fontsize, color=INK, va="center", ha="left")
    # cluster names, each followed by a bar spanning its rows
    term_to_cluster = dict(zip(sh.loc[rows, "term"], clusters))
    label_of = {term_to_cluster[k]: v for k, v in CLUSTER_LABELS.items()}
    assert sorted(label_of) == sorted(set(clusters)), "every cluster needs exactly one CLUSTER_LABELS anchor"
    for c, lab in label_of.items():
        idx = [i for i, cc in enumerate(clusters) if cc == c]
        lo, hi = min(idx), max(idx)
        ax_l.plot([x_bar, x_bar], [lo - 0.35, hi + 0.35], color=cols[c], lw=1.0, solid_capstyle="butt")
        ax_l.text(x_bar - 0.012, (lo + hi) / 2, lab, fontsize=fontsize, color=INK, va="center", ha="right",
                  linespacing=1.2)
    ax_l.set_zorder(ax_t.get_zorder() - 1)
    return x_lab   # left edge of the light term-label blocks (fraction of ax)


def _pathway_dotplot(ax, y0=0.02, h=0.905):
    container(ax)
    sh = load_shared_pathways()
    _, rows = shared_row_order()
    n = len(rows)
    sh = sh.loc[rows]
    x_dot, w_dot = 0.0, 0.34
    ax_d = inset(ax, x_dot, y0, w_dot, h)
    ax_d.set_xlim(-0.6, 1.6); ax_d.set_ylim(n - 0.5, -0.5)
    _cluster_bands(ax_d, sh["cluster"].tolist(), x0=-0.6, x1=1.6)
    logp = np.concatenate([-np.log10(sh.rna_p), -np.log10(sh.meth_p)])
    norm = Normalize(0, np.ceil(logp.max()))
    counts = np.concatenate([sh.rna_n_DE, sh.nSigGene])
    s_max = counts.max()
    size = lambda c: 4 + 56 * (np.asarray(c) / s_max)
    for j, (pcol, ccol) in enumerate((("rna_p", "rna_n_DE"), ("meth_p", "nSigGene"))):
        ax_d.scatter(np.full(n, j), np.arange(n), s=size(sh[ccol]), c=-np.log10(sh[pcol]), cmap=SEQ_CMAP, norm=norm,
                     edgecolors="white", linewidths=0.3, zorder=3)
    ax_d.set_xticks([0, 1]); ax_d.set_xticklabels(["RNA-seq", "Methyl."], fontsize=fs(5.5))
    ax_d.xaxis.tick_top()
    ax_d.tick_params(axis="x", length=0, pad=1.5)
    ax_d.set_yticks([])
    for s in ax_d.spines.values():
        s.set_visible(False)
    # legends: colorbar + size key
    x_leg = x_dot + w_dot + 0.08
    cax = inset(ax, x_leg, y0 + h - 0.26, 0.05, 0.20)
    z_colorbar(cax, label="−log10 p", cmap=SEQ_CMAP, vmin=norm.vmin, vmax=norm.vmax,
               ticks=np.arange(0, norm.vmax + 1, 2))
    kax = inset(ax, x_leg, y0 + h - 0.66, 0.5, 0.34)
    kax.set_xlim(0, 1); kax.set_ylim(0, 1); kax.set_axis_off()
    kax.text(0, 1.0, "Genes in term\nDE (RNA) or\nDMR (methyl.)", fontsize=fs(5.5), va="top", ha="left", color=INK)
    for k, c in enumerate([5, 10, 20]):
        yk = 0.52 - k * 0.2
        kax.scatter([0.08], [yk], s=size(c), color=MUTED, edgecolors="white", linewidths=0.3, clip_on=False)
        kax.text(0.25, yk, str(c), fontsize=fs(5.5), va="center", color=INK2)
    return x_dot + w_dot   # right edge of the dot columns / cluster blocks (fraction of ax)


def plot_e_shared_pathways(ax, tree_w=4.65, dot_w=2.55):
    """Tree + term labels on the left, dotplot on the right, sharing rows and cluster blocks.
    tree_w / dot_w: relative widths of the two parts (inches at the composite size)."""
    container(ax)
    f = tree_w / (tree_w + dot_w)
    x_left = _pathway_tree(inset(ax, 0, 0, f, 1)) * f
    x_right = f + _pathway_dotplot(inset(ax, f, 0, 1 - f, 1)) * (1 - f)
    panel_title(ax, "Shared GO:BP pathways", x=(x_left + x_right) / 2)   # centered on term labels + dotplot
    return ax


# -------------------------------------------------------------- panel F ----
def plot_f_geneset_heatmap(ax, set_name=GENESET):
    container(ax)
    genes = load_candidate_set(set_name)
    rna = pd.read_csv(DATA / "geneset_rna_corrected.csv").set_index("symbol")
    meth = pd.read_csv(DATA / "geneset_meth_corrected.csv").set_index("symbol")
    rmeta = order_samples(pd.read_csv(DATA / "geneset_rna_samples.csv"), ["EV07", "EV17", "EV18", "EV19", "EV20"])
    mmeta = order_samples(pd.read_csv(DATA / "geneset_meth_samples.csv"), ["EV07", "EV08", "ZY03"])
    corr = pd.read_csv(INT_RES / "candidate_genes/rna_meth_dose_trajectory_correlation.csv").set_index("symbol")
    de = pd.read_csv(INT_RES / "rna_gene_slopes.csv").set_index("symbol")["DE_sig"]
    dmr = pd.read_csv(INT_RES / "methylation_gene_slopes.csv").set_index("symbol")["is_sig_DMR_gene"]

    rz = row_z(rna.reindex(genes)[rmeta["sample"]].to_numpy())
    mz = row_z(meth.reindex(genes)[mmeta["sample"]].to_numpy())
    r = corr["r"].reindex(genes).to_numpy()[:, None]

    n_r, n_c, n_m = len(rmeta), len(mmeta), 1
    x0, y0, h = 0.085, 0.04, 0.70
    unit = 0.575 / (n_r + n_c + 1.6)                       # width per sample column
    gap = 0.03
    w_r, w_m, w_c = n_r * unit, n_c * unit, 1.6 * unit
    xr, xm = x0, x0 + w_r + gap
    xc = xm + w_m + gap
    strip_h, sg = 0.045, 0.012
    h_pts = ax.get_position().height * ax.figure.get_figheight() * 72
    y_hdr = y0 + h + 3 * sg + 2 * strip_h + 0.01 + 1.25 * fs(6) / h_pts      # top of the block headers
    for xx, ww, mat, meta, title in ((xr, w_r, rz, rmeta, "Gene expression (RNA-seq)"),
                                     (xm, w_m, mz, mmeta, "DNA methylation")):
        a = inset(ax, xx, y0, ww, h)
        draw_matrix(a, mat)
        dose_block_gaps(a, meta, len(genes))
        annotation_strip(inset(ax, xx, y0 + h + sg, ww, strip_h), meta, "donor", DONOR_COLORS,
                         "Donor" if meta is rmeta else None, label_side="left")
        annotation_strip(inset(ax, xx, y0 + h + 2 * sg + strip_h, ww, strip_h), meta, "dose_nM", DOSE_COLORS,
                         "Dose" if meta is rmeta else None, label_side="left")
        ax.text(xx + ww / 2, y_hdr, title, transform=ax.transAxes, fontsize=fs(6), color=INK, va="top", ha="center")
    ac = inset(ax, xc, y0, w_c, h)
    draw_matrix(ac, r, vmin=-1, vmax=1)
    for i, v in enumerate(r[:, 0]):
        if np.isfinite(v):
            ac.text(0, i, f"{v:.2f}", ha="center", va="center", fontsize=fs(5),
                    color="white" if abs(v) > 0.6 else INK)
    ax.text(xc + w_c / 2, y_hdr, "Dose trend\ncorrelation", transform=ax.transAxes,   # top aligned with the others
            fontsize=fs(6), color=INK, va="top", ha="center")
    # gene labels with DE (*) / DMR (†) marks
    for i, g in enumerate(genes):
        mark = ("*" if de.get(g, False) else "") + ("†" if dmr.get(g, False) else "")
        bold = bool(mark)
        ax.text(x0 - 0.006, y0 + h * (1 - (i + 0.5) / len(genes)), g + mark, transform=ax.transAxes,
                fontsize=fs(7), ha="right", va="center", style="italic", fontweight="bold" if bold else "normal")
    name = GENESET_NAMES.get(set_name, set_name)
    panel_title(ax, name[0].upper() + name[1:], x=(xr + xc + w_c) / 2)   # centered over the heatmaps, not legends
    # legends
    xl = xc + w_c + 0.035
    top = y0 + h + 2 * (sg + strip_h)
    cb_top = top - (fs(5.5) + 4) / h_pts   # colorbar title (text + 3pt pad) then starts level with "Donor"
    z_colorbar(inset(ax, xl, cb_top - 0.27, 0.012, 0.27), label="Row z-score")
    z_colorbar(inset(ax, xl, y0 + 0.02, 0.012, 0.27), label="Pearson r", vmin=-1, vmax=1, ticks=(-1, 0, 1))
    sample_legends(inset(ax, xl + 0.085, y0 + 0.13, 0.13, top - y0 - 0.13),
                   ["EV07", "EV17", "EV18", "EV19", "EV20", "EV08", "ZY03"])
    if any(de.get(g, False) or dmr.get(g, False) for g in genes):   # key only when some gene is marked
        ax.text(xl + 0.085, y0, "* DE gene\n† probe in a DMR\nbold = either",
                transform=ax.transAxes, fontsize=fs(5), color=INK2, va="bottom", ha="left")
    return ax


GENESET_NAMES = {
    "01_Effector_Function": "effector function", "02_TCell_Signaling_Activation": "T-cell signaling and activation",
    "03_Exhaustion_Chronic_Stimulation": "exhaustion / chronic stimulation",
    "04_Memory_Stemness_Differentiation": "memory, stemness and differentiation", "05_Cytokines": "cytokines",
    "06_Proliferation_Survival": "proliferation and survival", "07_Immune_Regulation": "immune regulation",
    "08_Epigenetics": "epigenetic regulators", "09_Metabolic_Fitness": "metabolic fitness",
    "10_Trafficking": "trafficking",
}


def load_candidate_set(set_name):
    txt = paths.CANDIDATE_GENE_SETS_R.read_text()
    import re
    m = re.search(rf'"{set_name}"\s*=\s*c\(([^)]*)\)', txt)
    return re.findall(r'"([^"]+)"', m.group(1))


# ------------------------------------------------------------- assembly ----
# Slot sizes in inches (width, height); standalone panels use the same size.
W = 7.2
SLOTS = {
    "A": (2.75, 3.1), "B": (2.45, 3.1), "C": (2.0, 1.5), "D": (2.0, 1.6),
    "E": (7.2, 3.5), "F": (7.2, 2.2),
}

PANELS = {
    "A": ("Fig2A_de_heatmap", plot_a_de_heatmap), "B": ("Fig2B_dmr_heatmap", plot_b_dmr_heatmap),
    "C": ("Fig2C_gene_venn", plot_c_gene_venn), "D": ("Fig2D_pathway_venn", plot_d_pathway_venn),
    "E": ("Fig2E_shared_pathways", plot_e_shared_pathways), "F": ("Fig2F_geneset_heatmap", plot_f_geneset_heatmap),
}


def composite():
    top, mid, bot = SLOTS["A"][1], SLOTS["E"][1], SLOTS["F"][1]
    pad = 0.12
    margin = 0.06
    H = top + mid + bot + 2 * pad + margin
    fig = plt.figure(figsize=(W, H))

    def slot(x, y_top, w, h):
        return fig.add_axes([x / W, 1 - (margin + y_top + h) / H, w / W, h / H])

    axes = {
        "A": slot(0, 0, *SLOTS["A"]),
        "B": slot(SLOTS["A"][0], 0, *SLOTS["B"]),
        "C": slot(W - SLOTS["C"][0], 0, *SLOTS["C"]),
        "D": slot(W - SLOTS["D"][0], SLOTS["C"][1] + 0.05, SLOTS["D"][0], SLOTS["D"][1] - 0.05),
        "E": slot(0, top + pad, *SLOTS["E"]),
        "F": slot(0, top + mid + 2 * pad, *SLOTS["F"]),
    }
    for k, (_, fn) in PANELS.items():
        fn(axes[k])
        panel_letter(axes[k], k)
    out_dir, suffix = output_location()
    for ext in ("pdf", "svg"):
        fig.savefig(out_dir / f"{STEM}{suffix}.{ext}")
    return fig


def output_location():
    """(folder, filename suffix) for this run: a --tag release folder, a --font-bump trial, or the working files."""
    if TAG:
        d = FIG2_DIR / f"release_{TAG}"
        (d / "panels").mkdir(parents=True, exist_ok=True)
        return d, f"_{TAG}"
    return FIG2_DIR, (f"_font+{FONT_BUMP:g}" if FONT_BUMP else "")


def standalone():
    out_dir, suffix = output_location()
    panel_dir = out_dir / "panels" if TAG else PANEL_DIR
    panel_dir.mkdir(parents=True, exist_ok=True)
    for k, (name, fn) in PANELS.items():
        fig = plt.figure(figsize=SLOTS[k])
        ax = fig.add_axes([0, 0, 1, 0.94])
        fn(ax)
        for ext in ("pdf", "svg"):
            fig.savefig(panel_dir / f"{name}{suffix if TAG else ''}.{ext}")
        plt.close(fig)


if __name__ == "__main__":
    fig = composite()
    if "--preview" in sys.argv:
        fig.savefig(sys.argv[sys.argv.index("--preview") + 1], dpi=200)
    plt.close(fig)
    out_dir, suffix = output_location()
    if FONT_BUMP and not TAG:
        print("Saved", out_dir / f"{STEM}{suffix}.{{pdf,svg}} (standalone panels not written)")
    else:
        standalone()
        print("Saved", out_dir / f"{STEM}{suffix}.{{pdf,svg}}", "and panels in",
              out_dir / "panels" if TAG else PANEL_DIR)
