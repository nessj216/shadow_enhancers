#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from pathlib import Path


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, VPacker, HPacker, TextArea, DrawingArea
from pathlib import Path
from scipy.stats import norm
#from __future__ import annotations
from matplotlib.offsetbox import AnchoredOffsetbox, VPacker, HPacker, TextArea, DrawingArea
import matplotlib.patches as patches

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "input"
OUTPUT_FILES = ROOT / "outputs" / "files"
OUTPUT_FIGURES = ROOT / "outputs" / "figures"
OUTPUT_FILES.mkdir(parents=True, exist_ok=True)
OUTPUT_FIGURES.mkdir(parents=True, exist_ok=True)

# =========================================================
# FILE INPUTS — MOUSE
# =========================================================
genome_file = str(INPUT_DIR / "FINAL_TE_genomefile_merged_dedup.bed")

shadow_file = str(INPUT_DIR / "Shadows_mouse_TEcooption_lastcol_Final.bed")

single_file = str(INPUT_DIR / "singlesmouse_TEcooption.bed")

outdir = OUTPUT_FIGURES


# =========================================================
# PARAMETERS
# =========================================================
FLANK_BP = 20000      # ±20 kb neighborhood
USE_FDR = True
EPS = 1e-9

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# No Helitron/RC for mouse here
cat_order = ["LTR", "LINE", "SINE", "DNA-TIR"]
te_colors = ["#c6dbef", "#6baed6", "#1f78b4", "#f16913"]


# =========================================================
# TE CLASS MAPPING
# =========================================================
def map_te_bin(te_type):
    te_type = str(te_type)

    if te_type.startswith("LTR"):
        return "LTR"

    if te_type.startswith("LINE"):
        return "LINE"

    if te_type.startswith("SINE"):
        return "SINE"

    # Conservative DNA-TIR mouse families
    if te_type.startswith((
        "DNA/hAT",
        "DNA/TcMar",
        "DNA/PiggyBac",
        "DNA/MULE",
        "DNA/MuDR"
    )):
        return "DNA-TIR"

    return None


# =========================================================
# LOADING
# =========================================================
def load_background(path):
    """
    Mouse genome TE BED:
    chrom, start, end, TE_name, TE_class
    """
    bg = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "te_name", "te_type"]
    )

    bg["start"] = pd.to_numeric(bg["start"], errors="coerce")
    bg["end"] = pd.to_numeric(bg["end"], errors="coerce")

    bg = bg.dropna(subset=["chrom", "start", "end", "te_type"]).copy()
    bg["bp"] = bg["end"] - bg["start"]
    bg = bg[bg["bp"] > 0].copy()

    bg["te_bin"] = bg["te_type"].map(map_te_bin)
    bg = bg[bg["te_bin"].notna()].copy()

    return bg


def load_mouse_overlap(path):
    """
    Works for both mouse shadow and single overlap files.

    Assumes:
      enhancer = cols 0,1,2,3
      TE class = second-to-last column
      overlap bp = last column

    No 50 bp filter is applied here.
    """
    raw = pd.read_csv(path, sep="\t", header=None)

    df = pd.DataFrame({
        "enh_chrom": raw.iloc[:, 0],
        "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "enh_id": raw.iloc[:, 3].astype(str),
        "te_type": raw.iloc[:, -2],
        "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
    })

    df = df.dropna(subset=[
        "enh_chrom",
        "enh_start",
        "enh_end",
        "enh_id",
        "te_type",
        "overlap_bp"
    ]).copy()

    df = df[df["te_type"] != "."].copy()

    df["te_bin"] = df["te_type"].map(map_te_bin)
    df = df[df["te_bin"].notna()].copy()

    return df


def unique_enhancers(df):
    enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
    enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
    enh["win_end"] = enh["enh_end"] + FLANK_BP
    return enh


# =========================================================
# BP-BASED COMPOSITIONS
# =========================================================
def observed_overlap_bp_by_category(df):
    """
    Observed TE contribution inside enhancers, measured by overlap bp.
    """
    return (
        df.groupby("te_bin")["overlap_bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def observed_overlap_composition(df):
    x = observed_overlap_bp_by_category(df)
    return x / x.sum() * 100 if x.sum() > 0 else x


def genome_bp_by_category(bg):
    """
    Genome-wide TE background, measured by TE bp.
    """
    return (
        bg.groupby("te_bin")["bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def genome_composition(bg):
    x = genome_bp_by_category(bg)
    return x / x.sum() * 100 if x.sum() > 0 else x


def neighborhood_bp_per_enhancer(enhancers, bg):
    """
    For each enhancer window, calculate how many TE bp from each class
    overlap the ±FLANK_BP neighborhood.
    """
    rows = []

    for chrom, enh_chr in enhancers.groupby("enh_chrom"):
        te_chr = bg[bg["chrom"] == chrom].copy()

        te_st = te_chr["start"].to_numpy()
        te_en = te_chr["end"].to_numpy()
        te_bin = te_chr["te_bin"].to_numpy()

        for _, row in enh_chr.iterrows():
            ws = row["win_start"]
            we = row["win_end"]

            mask = (te_en > ws) & (te_st < we)

            sub_st = te_st[mask]
            sub_en = te_en[mask]
            sub_bin = te_bin[mask]

            if len(sub_st) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            ov = np.minimum(sub_en, we) - np.maximum(sub_st, ws)
            valid = ov > 0

            ov = ov[valid]
            sub_bin = sub_bin[valid]

            if len(ov) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            tmp = pd.DataFrame({
                "te_bin": sub_bin,
                "bp": ov
            })

            sums = (
                tmp.groupby("te_bin")["bp"]
                .sum()
                .reindex(cat_order, fill_value=0)
                .astype(float)
            )

            sums.name = row["enh_id"]
            rows.append(sums)

    out = pd.DataFrame(rows)
    out.index = enhancers["enh_id"].tolist()
    return out.reindex(columns=cat_order, fill_value=0)


def mean_local_bp_composition(enhancers, bg):
    """
    Mean per-enhancer local neighborhood composition, measured by bp.

    This first converts each enhancer's local neighborhood to percentages,
    then averages those percentages across enhancers.
    """
    per_enh_bp = neighborhood_bp_per_enhancer(enhancers, bg)

    per_enh_pct = per_enh_bp.div(
        per_enh_bp.sum(axis=1).replace(0, np.nan),
        axis=0
    ) * 100

    return per_enh_pct.mean(axis=0).reindex(cat_order, fill_value=0)


def pooled_local_bp_by_category(enhancers, bg):
    """
    Pooled local neighborhood bp across all enhancer windows.
    Useful for statistics.
    """
    per_enh_bp = neighborhood_bp_per_enhancer(enhancers, bg)
    return per_enh_bp.sum(axis=0).reindex(cat_order, fill_value=0).astype(float)


# =========================================================
# STATS
# =========================================================
def two_prop_ztest(count1, n1, count2, n2):
    if n1 <= 0 or n2 <= 0:
        return np.nan

    p1 = count1 / n1
    p2 = count2 / n2

    p_pool = (count1 + count2) / (n1 + n2)
    se = np.sqrt(p_pool * (1 - p_pool) * ((1 / n1) + (1 / n2)))

    if se == 0 or np.isnan(se):
        return np.nan

    z = (p1 - p2) / se
    p = 2 * norm.sf(abs(z))

    return p


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(len(pvals), np.nan)

    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return out

    p = pvals[valid]
    order = np.argsort(p)
    ranked = p[order]

    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    tmp = np.empty(len(q))
    tmp[order] = q
    out[valid] = tmp

    return out


def p_to_star(p):
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def category_bp_tests(obs_bp, local_bp, group_name):
    obs_total = obs_bp.sum()
    local_total = local_bp.sum()

    rows = []

    for cat in cat_order:
        p = two_prop_ztest(
            obs_bp[cat],
            obs_total,
            local_bp[cat],
            local_total
        )

        obs_prop = obs_bp[cat] / obs_total if obs_total > 0 else np.nan
        local_prop = local_bp[cat] / local_total if local_total > 0 else np.nan

        rows.append({
            "group": group_name,
            "te_class": cat,
            "obs_bp": obs_bp[cat],
            "local_bp": local_bp[cat],
            "obs_prop": obs_prop,
            "local_prop": local_prop,
            "log2_obs_over_local": np.log2((obs_prop + EPS) / (local_prop + EPS)),
            "p_value": p
        })

    out = pd.DataFrame(rows)
    out["p_adj"] = bh_fdr(out["p_value"].values) if USE_FDR else out["p_value"]
    out["stars"] = out["p_adj"].map(p_to_star)

    return out


# =========================================================
# PLOTTING
# =========================================================
def plot_stacked_bar(ax, values, xpos, width=0.62):
    bottom = 0

    for val, color in zip(values, te_colors):
        ax.bar(
            xpos,
            val,
            width=width,
            bottom=bottom,
            color=color,
            edgecolor="black"
        )

        if val > 0:
            if val < 5:
                ax.text(
                    xpos + 0.34,
                    bottom + val / 2,
                    f"{val:.1f}%",
                    ha="left",
                    va="center",
                    fontsize=10
                )
            else:
                ax.text(
                    xpos,
                    bottom + val / 2,
                    f"{val:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=10
                )

        bottom += val


def make_grouped_te_legend(ax):
    box_w = 34
    box_text_sep = 6

    def legend_row(color, label):
        da = DrawingArea(box_w, 14, 0, 0)

        rect = patches.Rectangle(
            (0, 2),
            30,
            10,
            facecolor=color,
            edgecolor="black",
            linewidth=1.2
        )

        da.add_artist(rect)

        txt = TextArea(
            label,
            textprops=dict(size=12, family="Arial")
        )

        return HPacker(
            children=[da, txt],
            align="center",
            pad=0,
            sep=box_text_sep
        )

    def header_row(label):
        spacer = DrawingArea(box_w + box_text_sep, 1, 0, 0)

        txt = TextArea(
            label,
            textprops=dict(size=14, weight="bold", family="Arial")
        )

        return HPacker(
            children=[spacer, txt],
            align="center",
            pad=0,
            sep=0
        )

    legend_box = VPacker(
        children=[
            header_row("RNA"),
            legend_row(te_colors[0], "LTR"),
            legend_row(te_colors[1], "LINE"),
            legend_row(te_colors[2], "SINE"),
            header_row("DNA"),
            legend_row(te_colors[3], "TIR"),
        ],
        align="left",
        pad=0,
        sep=4
    )

    anchored_box = AnchoredOffsetbox(
        loc="upper left",
        child=legend_box,
        pad=0.3,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        bbox_transform=ax.transAxes,
        borderpad=0
    )

    ax.add_artist(anchored_box)


def make_plot(summary):
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.grid(False)

    all_values = [
        summary["shadow_observed_pct"].values,
        summary["shadow_local_pct"].values,
        summary["single_observed_pct"].values,
        summary["single_local_pct"].values,
        summary["genome_pct"].values
    ]

    all_labels = [
        "shadows",
        "shadow\nneighborhood",
        "singles",
        "single\nneighborhood",
        "genome"
    ]

    x_positions = [0, 1, 2.4, 3.4, 4.8]

    for xpos, vals in zip(x_positions, all_values):
        plot_stacked_bar(ax, vals, xpos, width=0.62)

    ax.set_ylabel("TE bp %")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(all_labels, fontsize=14)
    ax.set_ylim(0, 100)

    make_grouped_te_legend(ax)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)

    ax.tick_params(
        axis="y",
        which="both",
        left=True,
        right=False,
        length=4,
        width=1,
        color="black",
        labelcolor="black"
    )

    ax.tick_params(
        axis="x",
        which="both",
        bottom=True,
        top=False,
        length=4,
        width=1,
        color="black",
        labelcolor="black"
    )

    plt.tight_layout()

    fig.savefig(
        outdir / "mouse_stacked_TE_observed_local_genome_BP_based.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.1
    )

    fig.savefig(
        outdir / "mouse_stacked_TE_observed_local_genome_BP_based.pdf",
        bbox_inches="tight",
        pad_inches=0.1
    )

    plt.show()


# =========================================================
# MAIN
# =========================================================
def main():
    shadow_df = load_mouse_overlap(shadow_file)
    single_df = load_mouse_overlap(single_file)
    genome_df = load_background(genome_file)

    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    print("Unique shadow enhancers:", len(shadow_enh))
    print("Unique single enhancers:", len(single_enh))

    # Percent compositions for plotting
    shadow_obs_pct = observed_overlap_composition(shadow_df)
    single_obs_pct = observed_overlap_composition(single_df)

    shadow_local_pct = mean_local_bp_composition(shadow_enh, genome_df)
    single_local_pct = mean_local_bp_composition(single_enh, genome_df)

    genome_pct = genome_composition(genome_df)

    # Raw bp counts for stats
    shadow_obs_bp = observed_overlap_bp_by_category(shadow_df)
    single_obs_bp = observed_overlap_bp_by_category(single_df)

    shadow_local_bp = pooled_local_bp_by_category(shadow_enh, genome_df)
    single_local_bp = pooled_local_bp_by_category(single_enh, genome_df)

    genome_bp = genome_bp_by_category(genome_df)

    shadow_stats = category_bp_tests(
        shadow_obs_bp,
        shadow_local_bp,
        group_name="shadow"
    )

    single_stats = category_bp_tests(
        single_obs_bp,
        single_local_bp,
        group_name="single"
    )

    stats = pd.concat([shadow_stats, single_stats], ignore_index=True)

    summary = pd.DataFrame({
        "shadow_observed_pct": shadow_obs_pct,
        "shadow_local_pct": shadow_local_pct,
        "single_observed_pct": single_obs_pct,
        "single_local_pct": single_local_pct,
        "genome_pct": genome_pct,
        "shadow_observed_bp": shadow_obs_bp,
        "shadow_local_bp": shadow_local_bp,
        "single_observed_bp": single_obs_bp,
        "single_local_bp": single_local_bp,
        "genome_bp": genome_bp
    }).reindex(cat_order)

    print("\nBP-based composition summary:")
    print(summary.round(4))

    print("\nObserved vs local BP tests:")
    print(stats.round(6))

    summary.to_csv(
        outdir / "mouse_stacked_TE_BP_based_summary.tsv",
        sep="\t"
    )

    stats.to_csv(
        outdir / "mouse_TEclass_observed_vs_local_BP_ztests.tsv",
        sep="\t",
        index=False
    )

    make_plot(summary)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, chisquare, power_divergence, chi2 as chi2_dist
from pathlib import Path

from matplotlib.offsetbox import AnchoredOffsetbox, VPacker, HPacker, TextArea, DrawingArea
import matplotlib.patches as patches
import matplotlib.patches as mpatches


# =========================================================
# FILE INPUTS — MOUSE
# =========================================================
genome_file = str(INPUT_DIR / "FINAL_TE_genomefile_merged_dedup.bed")

shadow_file = str(INPUT_DIR / "Shadows_mouse_TEcooption_lastcol_Final.bed")

single_file = str(INPUT_DIR / "singlesmouse_TEcooption.bed")

outdir = OUTPUT_FIGURES


# =========================================================
# PARAMETERS
# =========================================================
FLANK_BP = 20000
USE_FDR = True
EPS = 1e-9

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

cat_order = ["LTR", "LINE", "SINE", "DNA-TIR"]
te_colors = ["#c6dbef", "#6baed6", "#1f78b4", "#f16913"]

te_color_map = {
    "LTR": "#c6dbef",
    "LINE": "#6baed6",
    "SINE": "#1f78b4",
    "DNA-TIR": "#f16913",
}


# =========================================================
# TE CLASS MAPPING
# =========================================================
def map_te_bin(te_type):
    te_type = str(te_type)

    if te_type.startswith("LTR"):
        return "LTR"

    if te_type.startswith("LINE"):
        return "LINE"

    if te_type.startswith("SINE"):
        return "SINE"

    # Mouse DNA transposons treated as DNA-TIR.
    # This intentionally excludes RC/Helitron-like DNA classes.
    if te_type.startswith((
        "DNA/hAT",
        "DNA/TcMar",
        "DNA/PiggyBac",
        "DNA/MULE",
        "DNA/MuDR"
    )):
        return "DNA-TIR"

    return None


# =========================================================
# LOADING
# =========================================================
def load_background(path):
    """
    Mouse genome TE BED:
        chrom, start, end, TE_name, TE_class

    BP-based analysis:
        No insertion deduplication is done here.
    """
    bg = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "te_name", "te_type"]
    )

    bg["start"] = pd.to_numeric(bg["start"], errors="coerce")
    bg["end"] = pd.to_numeric(bg["end"], errors="coerce")

    bg = bg.dropna(subset=["chrom", "start", "end", "te_type"]).copy()
    bg["bp"] = bg["end"] - bg["start"]
    bg = bg[bg["bp"] > 0].copy()

    bg["te_bin"] = bg["te_type"].map(map_te_bin)
    bg = bg[bg["te_bin"].notna()].copy()

    return bg


def load_mouse_overlap(path):
    """
    Works for both mouse shadow and single overlap files.

    Assumes:
        enhancer = cols 0,1,2,3
        TE class = second-to-last column
        overlap bp = last column

    This analysis is BP-based.
    No TE insertion deduplication is done.
    No extra 50 bp filter is applied here because your files are already filtered.
    """
    raw = pd.read_csv(path, sep="\t", header=None)

    df = pd.DataFrame({
        "enh_chrom": raw.iloc[:, 0],
        "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "enh_id": raw.iloc[:, 3].astype(str),
        "te_type": raw.iloc[:, -2],
        "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
    })

    df = df.dropna(subset=[
        "enh_chrom",
        "enh_start",
        "enh_end",
        "enh_id",
        "te_type",
        "overlap_bp"
    ]).copy()

    df = df[df["te_type"] != "."].copy()
    df = df[df["overlap_bp"] > 0].copy()

    df["te_bin"] = df["te_type"].map(map_te_bin)
    df = df[df["te_bin"].notna()].copy()

    return df


def unique_enhancers(df):
    enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
    enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
    enh["win_end"] = enh["enh_end"] + FLANK_BP
    return enh


# =========================================================
# BP-BASED COUNTS AND COMPOSITIONS
# =========================================================
def observed_overlap_bp_by_category(df):
    """
    Observed TE contribution inside enhancers, measured by overlap bp.
    """
    return (
        df.groupby("te_bin")["overlap_bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def composition_from_counts(counts):
    counts = pd.Series(counts, index=cat_order, dtype=float)
    total = counts.sum()

    if total <= 0:
        return pd.Series(np.nan, index=cat_order)

    return counts / total * 100


def genome_bp_by_category(bg):
    """
    Genome-wide TE background, measured by TE bp.
    """
    return (
        bg.groupby("te_bin")["bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def neighborhood_bp_per_enhancer(enhancers, bg):
    """
    For each enhancer window, calculate TE bp from each class
    overlapping the ±FLANK_BP neighborhood.

    BP-based analysis:
        No TE insertion deduplication is done.
    """
    rows = []

    for chrom, enh_chr in enhancers.groupby("enh_chrom"):
        te_chr = bg[bg["chrom"] == chrom].copy()

        te_st = te_chr["start"].to_numpy()
        te_en = te_chr["end"].to_numpy()
        te_bin = te_chr["te_bin"].to_numpy()

        for _, row in enh_chr.iterrows():
            ws = row["win_start"]
            we = row["win_end"]

            mask = (te_en > ws) & (te_st < we)

            sub_st = te_st[mask]
            sub_en = te_en[mask]
            sub_bin = te_bin[mask]

            if len(sub_st) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            ov = np.minimum(sub_en, we) - np.maximum(sub_st, ws)
            valid = ov > 0

            ov = ov[valid]
            sub_bin = sub_bin[valid]

            if len(ov) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            tmp = pd.DataFrame({
                "te_bin": sub_bin,
                "bp": ov
            })

            sums = (
                tmp.groupby("te_bin")["bp"]
                .sum()
                .reindex(cat_order, fill_value=0)
                .astype(float)
            )

            sums.name = row["enh_id"]
            rows.append(sums)

    out = pd.DataFrame(rows)
    out.index = enhancers["enh_id"].tolist()

    return out.reindex(columns=cat_order, fill_value=0)


def pooled_local_bp_by_category(enhancers, bg):
    """
    Pooled local neighborhood bp across all enhancer windows.

    This is the denominator used for:
        - stacked neighborhood bars
        - log2 Observed/Neighborhood
        - per-class tests
        - GOF tests
    """
    per_enh_bp = neighborhood_bp_per_enhancer(enhancers, bg)
    return per_enh_bp.sum(axis=0).reindex(cat_order, fill_value=0).astype(float)


# =========================================================
# STATS
# =========================================================
def two_prop_ztest(count1, n1, count2, n2):
    """
    Two-proportion z-test using bp counts.

    Returns:
        z, p, log10_p

    log10_p is included because p can underflow to 0.0 for large BP counts.
    """
    if n1 <= 0 or n2 <= 0:
        return np.nan, np.nan, np.nan

    p1 = count1 / n1
    p2 = count2 / n2

    p_pool = (count1 + count2) / (n1 + n2)
    se = np.sqrt(p_pool * (1 - p_pool) * ((1 / n1) + (1 / n2)))

    if se == 0 or np.isnan(se):
        return np.nan, np.nan, np.nan

    z = (p1 - p2) / se

    p = 2 * norm.sf(abs(z))
    log10_p = (np.log(2) + norm.logsf(abs(z))) / np.log(10)

    return z, p, log10_p


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(len(pvals), np.nan)

    valid = np.isfinite(pvals)

    if valid.sum() == 0:
        return out

    p = pvals[valid]
    order = np.argsort(p)
    ranked = p[order]

    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    tmp = np.empty(len(q))
    tmp[order] = q
    out[valid] = tmp

    return out


def p_to_star(p):
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def category_bp_tests(obs_bp, local_bp, group_name):
    """
    Per-class Observed vs Neighborhood BP tests.
    """
    obs_bp = pd.Series(obs_bp, index=cat_order, dtype=float)
    local_bp = pd.Series(local_bp, index=cat_order, dtype=float)

    obs_total = obs_bp.sum()
    local_total = local_bp.sum()

    rows = []

    for cat in cat_order:
        z, p, log10_p = two_prop_ztest(
            obs_bp[cat],
            obs_total,
            local_bp[cat],
            local_total
        )

        obs_prop = obs_bp[cat] / obs_total if obs_total > 0 else np.nan
        local_prop = local_bp[cat] / local_total if local_total > 0 else np.nan

        rows.append({
            "group": group_name,
            "te_class": cat,
            "obs_bp": obs_bp[cat],
            "local_bp": local_bp[cat],
            "obs_prop": obs_prop,
            "local_prop": local_prop,
            "log2_observed_over_neighborhood": np.log2((obs_prop + EPS) / (local_prop + EPS)),
            "z_value": z,
            "p_value": p,
            "log10_p": log10_p,
            "direction": "enriched" if obs_prop > local_prop else "depleted"
        })

    out = pd.DataFrame(rows)
    out["p_adj"] = bh_fdr(out["p_value"].values) if USE_FDR else out["p_value"]
    out["stars"] = out["p_adj"].map(p_to_star)

    return out


def run_gtest_and_chisq(obs_bp, ref_bp, label):
    """
    Overall BP-composition goodness-of-fit test.

    obs_bp:
        observed enhancer-overlap bp by TE class

    ref_bp:
        neighborhood bp or genome bp by TE class
    """
    obs_bp = pd.Series(obs_bp, index=cat_order, dtype=float)
    ref_bp = pd.Series(ref_bp, index=cat_order, dtype=float)

    obs_total = obs_bp.sum()
    ref_total = ref_bp.sum()

    ref_props = ref_bp / ref_total
    expected = ref_props * obs_total

    # Make scipy happy about exact total matching.
    expected = expected * (obs_total / expected.sum())

    keep = ~((obs_bp == 0) & (expected == 0))

    obs_use = obs_bp[keep]
    exp_use = expected[keep]

    if np.any((exp_use == 0) & (obs_use > 0)):
        summary = pd.DataFrame([{
            "comparison": label,
            "obs_total_bp": obs_total,
            "ref_total_bp": ref_total,
            "n_categories_tested": int(keep.sum()),
            "chi2_stat": np.nan,
            "chi2_p": np.nan,
            "chi2_log10_p": np.nan,
            "g_stat": np.nan,
            "g_p": np.nan,
            "g_log10_p": np.nan,
            "note": "Invalid because at least one category has expected=0 but observed>0"
        }])

        diagnostic = pd.DataFrame({
            "comparison": label,
            "te_class": cat_order,
            "observed_bp": obs_bp.values,
            "expected_bp": expected.values,
            "obs_prop": obs_bp.values / obs_total,
            "ref_prop": ref_props.values,
            "chi_contribution": np.nan
        })

        return summary, diagnostic

    chi_stat, chi_p = chisquare(
        f_obs=obs_use.values,
        f_exp=exp_use.values
    )

    g_stat, g_p = power_divergence(
        f_obs=obs_use.values,
        f_exp=exp_use.values,
        lambda_="log-likelihood"
    )

    df = len(obs_use) - 1

    chi_log10_p = chi2_dist.logsf(chi_stat, df) / np.log(10)
    g_log10_p = chi2_dist.logsf(g_stat, df) / np.log(10)

    diagnostic = pd.DataFrame({
        "comparison": label,
        "te_class": obs_use.index,
        "observed_bp": obs_use.values,
        "expected_bp": exp_use.values,
        "obs_prop": obs_use.values / obs_total,
        "ref_prop": ref_props[obs_use.index].values,
        "chi_contribution": ((obs_use.values - exp_use.values) ** 2) / exp_use.values
    })

    summary = pd.DataFrame([{
        "comparison": label,
        "obs_total_bp": obs_total,
        "ref_total_bp": ref_total,
        "n_categories_tested": int(keep.sum()),
        "chi2_stat": chi_stat,
        "chi2_p": chi_p,
        "chi2_log10_p": chi_log10_p,
        "g_stat": g_stat,
        "g_p": g_p,
        "g_log10_p": g_log10_p,
        "note": ""
    }])

    return summary, diagnostic

def composition_effect_size(obs_counts, ref_counts, label):
    obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
    ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)

    obs_total = obs_counts.sum()
    ref_total = ref_counts.sum()

    obs_props = obs_counts / obs_total
    ref_props = ref_counts / ref_total

    keep = ref_props > 0
    obs_props = obs_props[keep]
    ref_props = ref_props[keep]

    prop_chi2 = (((obs_props - ref_props) ** 2) / ref_props).sum()
    cohens_w = np.sqrt(prop_chi2)

    return pd.DataFrame([{
        "comparison": label,
        "obs_total": obs_total,
        "ref_total": ref_total,
        "proportional_chi2": prop_chi2,
        "cohens_w": cohens_w,
        "note": "Effect size from proportions only; not a chi-square p-value"
    }])
# =========================================================
# PLOTTING HELPERS
# =========================================================
def plot_stacked_bar(ax, values, xpos, width=0.62):
    bottom = 0

    for val, color in zip(values, te_colors):
        ax.bar(
            xpos,
            val,
            width=width,
            bottom=bottom,
            color=color,
            edgecolor="black"
        )

        if val > 0:
            if val < 5:
                ax.text(
                    xpos + 0.34,
                    bottom + val / 2,
                    f"{val:.1f}%",
                    ha="left",
                    va="center",
                    fontsize=15
                )
            else:
                ax.text(
                    xpos,
                    bottom + val / 2,
                    f"{val:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=15
                )

        bottom += val


def make_grouped_te_legend(ax):
    box_w = 34
    box_text_sep = 6

    def legend_row(color, label):
        da = DrawingArea(box_w, 14, 0, 0)

        rect = patches.Rectangle(
            (0, 2),
            30,
            10,
            facecolor=color,
            edgecolor="black",
            linewidth=1.2
        )

        da.add_artist(rect)

        txt = TextArea(
            label,
            textprops=dict(size=12, family="Arial")
        )

        return HPacker(
            children=[da, txt],
            align="center",
            pad=0,
            sep=box_text_sep
        )

    def header_row(label):
        spacer = DrawingArea(box_w + box_text_sep, 1, 0, 0)

        txt = TextArea(
            label,
            textprops=dict(size=14, weight="bold", family="Arial")
        )

        return HPacker(
            children=[spacer, txt],
            align="center",
            pad=0,
            sep=0
        )

    legend_box = VPacker(
        children=[
            header_row("RNA"),
            legend_row(te_colors[0], "LTR"),
            legend_row(te_colors[1], "LINE"),
            legend_row(te_colors[2], "SINE"),
            header_row("DNA"),
            legend_row(te_colors[3], "DNA-TIR"),
        ],
        align="left",
        pad=0,
        sep=4
    )

    anchored_box = AnchoredOffsetbox(
        loc="upper left",
        child=legend_box,
        pad=0.3,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        bbox_transform=ax.transAxes,
        borderpad=0
    )

    ax.add_artist(anchored_box)


def style_axes(ax):
    ax.grid(False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)

    ax.tick_params(
        axis="y",
        which="both",
        left=True,
        right=False,
        length=4,
        width=1,
        color="black",
        labelcolor="black"
    )

    ax.tick_params(
        axis="x",
        which="both",
        bottom=True,
        top=False,
        length=4,
        width=1,
        color="black",
        labelcolor="black"
    )


def make_stacked_plot(summary):
    fig, ax = plt.subplots(figsize=(6.5, 3.2))

    all_values = [
        summary["shadow_observed_pct"].values,
        summary["shadow_neighborhood_pct"].values,
        summary["single_observed_pct"].values,
        summary["single_neighborhood_pct"].values,
        summary["genome_pct"].values
    ]

    all_labels = [
        "shadows",
        "shadow\nneighborhood",
        "singles",
        "single\nneighborhood",
        "genome"
    ]

    x_positions = [0, 1, 2.4, 3.4, 4.8]

    for xpos, vals in zip(x_positions, all_values):
        plot_stacked_bar(ax, vals, xpos, width=0.62)

    ax.set_ylabel("TE bp %")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(all_labels, fontsize=15)
    ax.set_ylim(0, 100)

    make_grouped_te_legend(ax)
    style_axes(ax)

    plt.tight_layout()

    fig.savefig(
        outdir / "mouse_stacked_TE_observed_neighborhood_genome_BP_based.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.1
    )

    fig.savefig(
        outdir / "mouse_stacked_TE_observed_neighborhood_genome_BP_based.pdf",
        bbox_inches="tight",
        pad_inches=0.1
    )

    plt.show()


def make_log2_plot(stats):
    plot_df = stats.pivot(
        index="te_class",
        columns="group",
        values="log2_observed_over_neighborhood"
    ).reindex(cat_order)

    labels = cat_order

    shadow_vals = plot_df["shadow"].to_numpy()
    single_vals = plot_df["single"].to_numpy()

    x = np.arange(len(cat_order))
    bar_width = 0.46

    fig, ax = plt.subplots(figsize=(5.3, 3.2))

    bars1 = ax.bar(
        x - bar_width / 2,
        shadow_vals,
        width=bar_width,
        color=[te_color_map[lbl] for lbl in labels],
        edgecolor="black",
        hatch="//"
    )

    bars2 = ax.bar(
        x + bar_width / 2,
        single_vals,
        width=bar_width,
        color=[te_color_map[lbl] for lbl in labels],
        edgecolor="black"
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=15)

    ax.set_ylabel(
        r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$',
        fontsize=15
    )

    ax.axhline(0, color="black", linewidth=1)

    finite_vals = np.concatenate([
        shadow_vals[np.isfinite(shadow_vals)],
        single_vals[np.isfinite(single_vals)]
    ])

    if len(finite_vals) > 0:
        ymin = np.floor(finite_vals.min()) - 0.4
        ymax = np.ceil(finite_vals.max()) + 0.9
        ax.set_ylim(ymin, ymax)

    # numeric labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()

            if np.isfinite(height):
                offset = 0.1 if height >= 0 else -0.2

                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    height + offset,
                    f"{height:.2f}",
                    ha="center",
                    va="bottom" if height >= 0 else "top",
                    fontsize=16
                )

    # stars from FDR-adjusted per-class BP tests
    for i, cat in enumerate(cat_order):
        sh_row = stats[
            (stats["group"] == "shadow") &
            (stats["te_class"] == cat)
        ].iloc[0]

        si_row = stats[
            (stats["group"] == "single") &
            (stats["te_class"] == cat)
        ].iloc[0]

        if sh_row["stars"] != "ns":
            y = shadow_vals[i] + (0.32 if shadow_vals[i] >= 0 else -0.55)

            ax.text(
                x[i] - bar_width / 2,
                y,
                sh_row["stars"],
                ha="center",
                va="bottom" if shadow_vals[i] >= 0 else "top",
                fontsize=15,
                fontweight="bold"
            )

        if si_row["stars"] != "ns":
            y = single_vals[i] + (0.32 if single_vals[i] >= 0 else -0.55)

            ax.text(
                x[i] + bar_width / 2,
                y,
                si_row["stars"],
                ha="center",
                va="bottom" if single_vals[i] >= 0 else "top",
                fontsize=15,
                fontweight="bold"
            )

    shadow_patch = mpatches.Patch(
        facecolor="white",
        edgecolor="black",
        hatch="//",
        label="shadows"
    )

    single_patch = mpatches.Patch(
        facecolor="white",
        edgecolor="black",
        label="singles"
    )

    ax.legend(
        handles=[shadow_patch, single_patch],
        fontsize=15
    )

    style_axes(ax)

    plt.tight_layout()

    fig.savefig(
        outdir / "mouse_log2_observed_over_neighborhood_BP_based.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.1
    )

    fig.savefig(
        outdir / "mouse_log2_observed_over_neighborhood_BP_based.pdf",
        bbox_inches="tight",
        pad_inches=0.1
    )

    plt.show()


# =========================================================
# MAIN
# =========================================================
def main():
    shadow_df = load_mouse_overlap(shadow_file)
    single_df = load_mouse_overlap(single_file)
    genome_df = load_background(genome_file)

    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    print("Unique shadow enhancers:", len(shadow_enh))
    print("Unique single enhancers:", len(single_enh))

    print("\nMapped observed overlap rows:")
    print("Shadow rows:", len(shadow_df))
    print("Single rows:", len(single_df))

    print("\nMapped genome TE rows:")
    print("Genome rows:", len(genome_df))

    # Raw BP counts
    shadow_obs_bp = observed_overlap_bp_by_category(shadow_df)
    single_obs_bp = observed_overlap_bp_by_category(single_df)

    shadow_neighborhood_bp = pooled_local_bp_by_category(shadow_enh, genome_df)
    single_neighborhood_bp = pooled_local_bp_by_category(single_enh, genome_df)

    genome_bp = genome_bp_by_category(genome_df)

    # Percent compositions
    shadow_obs_pct = composition_from_counts(shadow_obs_bp)
    single_obs_pct = composition_from_counts(single_obs_bp)

    shadow_neighborhood_pct = composition_from_counts(shadow_neighborhood_bp)
    single_neighborhood_pct = composition_from_counts(single_neighborhood_bp)

    genome_pct = composition_from_counts(genome_bp)

    # Per-class BP tests
    shadow_stats = category_bp_tests(
        shadow_obs_bp,
        shadow_neighborhood_bp,
        group_name="shadow"
    )

    single_stats = category_bp_tests(
        single_obs_bp,
        single_neighborhood_bp,
        group_name="single"
    )

    stats = pd.concat([shadow_stats, single_stats], ignore_index=True)

    # Overall GOF tests
    gof_rows = []
    gof_diag_rows = []

    comparisons = [
        (
            shadow_obs_bp,
            shadow_neighborhood_bp,
            "Shadow observed bp vs Shadow neighborhood bp"
        ),
        (
            shadow_obs_bp,
            genome_bp,
            "Shadow observed bp vs Genome bp"
        ),
        (
            single_obs_bp,
            single_neighborhood_bp,
            "Single observed bp vs Single neighborhood bp"
        ),
        (
            single_obs_bp,
            genome_bp,
            "Single observed bp vs Genome bp"
        ),
    ]

    for obs_bp, ref_bp, label in comparisons:
        s, d = run_gtest_and_chisq(obs_bp, ref_bp, label)
        gof_rows.append(s)
        gof_diag_rows.append(d)

    def composition_effect_size(obs_counts, ref_counts, label):
        """
        Percent/proportion-based composition effect size.

        This is NOT a chi-square p-value.
        It is chi-square divided by N, also called Cohen's w^2.
        Cohen's w = sqrt(proportional_chi2).
        """
        obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
        ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)

        obs_total = obs_counts.sum()
        ref_total = ref_counts.sum()

        obs_props = obs_counts / obs_total
        ref_props = ref_counts / ref_total

        keep = ref_props > 0
        obs_props = obs_props[keep]
        ref_props = ref_props[keep]

        proportional_chi2 = (((obs_props - ref_props) ** 2) / ref_props).sum()
        cohens_w = np.sqrt(proportional_chi2)

        return pd.DataFrame([{
            "comparison": label,
            "obs_total": obs_total,
            "ref_total": ref_total,
            "proportional_chi2": proportional_chi2,
            "cohens_w": cohens_w,
            "note": "Effect size from proportions only; not a chi-square p-value"
        }])
    gof_stats = pd.concat(gof_rows, ignore_index=True)
    gof_diagnostics = pd.concat(gof_diag_rows, ignore_index=True)
    # Percent/proportion-based composition effect sizes
    effect_rows = []

    for obs_bp, ref_bp, label in comparisons:
        effect_rows.append(
            composition_effect_size(obs_bp, ref_bp, label)
        )

    effect_sizes = pd.concat(effect_rows, ignore_index=True)
    summary = pd.DataFrame({
        "shadow_observed_pct": shadow_obs_pct,
        "shadow_neighborhood_pct": shadow_neighborhood_pct,
        "single_observed_pct": single_obs_pct,
        "single_neighborhood_pct": single_neighborhood_pct,
        "genome_pct": genome_pct,

        "shadow_observed_bp": shadow_obs_bp,
        "shadow_neighborhood_bp": shadow_neighborhood_bp,
        "single_observed_bp": single_obs_bp,
        "single_neighborhood_bp": single_neighborhood_bp,
        "genome_bp": genome_bp
    }).reindex(cat_order)

    print("\nBP-based composition summary:")
    print(summary.round(4))

    print("\nObserved vs Neighborhood BP per-class tests:")
    print(stats.round(6))

    print("\nOverall BP composition GOF tests:")
    print(gof_stats.round(6))

    print("\nGOF diagnostics by TE class:")
    print(gof_diagnostics.round(6))
    print("\nBP composition effect sizes:")
    print(effect_sizes.round(6))
    # Save tables
    summary.to_csv(
        outdir / "mouse_stacked_TE_BP_based_summary.tsv",
        sep="\t"
    )

    stats.to_csv(
        outdir / "mouse_TEclass_observed_vs_neighborhood_BP_ztests.tsv",
        sep="\t",
        index=False
    )

    gof_stats.to_csv(
        outdir / "mouse_BP_based_gtest_chisq_stats.tsv",
        sep="\t",
        index=False
    )

    gof_diagnostics.to_csv(
        outdir / "mouse_BP_based_gof_diagnostics_by_TEclass.tsv",
        sep="\t",
        index=False
    )
    effect_sizes.to_csv(
        outdir / "mouse_BP_based_composition_effect_sizes.tsv",
        sep="\t",
        index=False
    )
    # Plots
    make_stacked_plot(summary)
    make_log2_plot(stats)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3



# =========================================================
# FILE INPUTS — MOUSE
# =========================================================
genome_file = str(INPUT_DIR / "FINAL_TE_genomefile_merged_dedup.bed")
shadow_file = str(INPUT_DIR / "Shadows_mouse_TEcooption_lastcol_Final.bed")
single_file = str(INPUT_DIR / "singlesmouse_TEcooption.bed")

outdir = OUTPUT_FIGURES


# =========================================================
# PARAMETERS
# =========================================================
FLANK_BP = 20000
USE_FDR = True
EPS = 1e-9

cat_order = ["LTR", "LINE", "SINE", "DNA-TIR"]
tissue_order = ["forebrain", "heart", "limb"]
te_colors = ["#c6dbef", "#6baed6", "#1f78b4", "#f16913"]
tissue_colors = {
    "forebrain": "#4C78A8",
    "heart": "#E45756",
    "limb": "#54A24B",
    "pooled": "black",
}

PANEL_HEIGHT = 3.0
STACKED_FIGSIZE = (6.5, 3.2)
TISSUE_SPLIT_FIGSIZE = (8.3, PANEL_HEIGHT)
ENRICHMENT_YLIM = (-3, 3.2)

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# =========================================================
# HELPERS
# =========================================================
def map_te_bin(te_type):
    te_type = str(te_type)

    if te_type.startswith("LTR"):
        return "LTR"
    if te_type.startswith("LINE"):
        return "LINE"
    if te_type.startswith("SINE"):
        return "SINE"
    if te_type.startswith((
        "DNA/hAT",
        "DNA/TcMar",
        "DNA/PiggyBac",
        "DNA/MULE",
        "DNA/MuDR",
    )):
        return "DNA-TIR"
    return None


def assign_tissue_from_enh_id(enh_id):
    s = str(enh_id).lower()
    if "forebrain" in s or "fb" in s:
        return "forebrain"
    if "heart" in s:
        return "heart"
    if "limb" in s:
        return "limb"
    return "other"


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(len(pvals), np.nan)

    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return out

    p = pvals[valid]
    order = np.argsort(p)
    ranked = p[order]

    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    tmp = np.empty(len(q))
    tmp[order] = q
    out[valid] = tmp
    return out


def p_to_star(p):
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def style_axes(ax):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(
        axis="both",
        which="both",
        left=True,
        right=False,
        bottom=True,
        top=False,
        length=4,
        width=1.0,
        color="black",
        labelcolor="black",
    )


def save_figure(fig, png_path=None, pdf_path=None, left=0.10, right=0.82, bottom=0.22, top=0.92):
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    if png_path:
        fig.savefig(png_path, dpi=600)
    if pdf_path:
        fig.savefig(pdf_path)
    plt.close(fig)


# =========================================================
# LOADING
# =========================================================
def load_background(path):
    bg = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "te_name", "te_type"],
    )

    bg["start"] = pd.to_numeric(bg["start"], errors="coerce")
    bg["end"] = pd.to_numeric(bg["end"], errors="coerce")

    bg = bg.dropna(subset=["chrom", "start", "end", "te_type"]).copy()
    bg["bp"] = bg["end"] - bg["start"]
    bg = bg[bg["bp"] > 0].copy()

    bg["te_bin"] = bg["te_type"].map(map_te_bin)
    bg = bg[bg["te_bin"].notna()].copy()
    return bg


def load_mouse_overlap(path):
    raw = pd.read_csv(path, sep="\t", header=None)

    df = pd.DataFrame({
        "enh_chrom": raw.iloc[:, 0],
        "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "enh_id": raw.iloc[:, 3].astype(str),
        "te_type": raw.iloc[:, -2],
        "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
    })

    df = df.dropna(subset=[
        "enh_chrom",
        "enh_start",
        "enh_end",
        "enh_id",
        "te_type",
        "overlap_bp",
    ]).copy()

    df = df[df["te_type"] != "."].copy()
    df = df[df["overlap_bp"] > 0].copy()

    df["te_bin"] = df["te_type"].map(map_te_bin)
    df = df[df["te_bin"].notna()].copy()
    df["tissue"] = df["enh_id"].map(assign_tissue_from_enh_id)
    return df


def unique_enhancers(df):
    enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
    enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
    enh["win_end"] = enh["enh_end"] + FLANK_BP
    enh["tissue"] = enh["enh_id"].map(assign_tissue_from_enh_id)
    return enh


# =========================================================
# BP-BASED COMPOSITIONS
# =========================================================
def observed_overlap_bp_by_category(df):
    return (
        df.groupby("te_bin")["overlap_bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def observed_overlap_composition(df):
    x = observed_overlap_bp_by_category(df)
    return x / x.sum() * 100 if x.sum() > 0 else x


def composition_from_counts(counts):
    counts = pd.Series(counts, index=cat_order, dtype=float)
    total = counts.sum()
    if total <= 0:
        return pd.Series(np.nan, index=cat_order)
    return counts / total * 100


def genome_bp_by_category(bg):
    return (
        bg.groupby("te_bin")["bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def genome_composition(bg):
    x = genome_bp_by_category(bg)
    return x / x.sum() * 100 if x.sum() > 0 else x


def neighborhood_bp_per_enhancer(enhancers, bg):
    rows = []

    for chrom, enh_chr in enhancers.groupby("enh_chrom"):
        te_chr = bg[bg["chrom"] == chrom].copy()

        te_st = te_chr["start"].to_numpy()
        te_en = te_chr["end"].to_numpy()
        te_bin = te_chr["te_bin"].to_numpy()

        for _, row in enh_chr.iterrows():
            ws = row["win_start"]
            we = row["win_end"]

            mask = (te_en > ws) & (te_st < we)
            sub_st = te_st[mask]
            sub_en = te_en[mask]
            sub_bin = te_bin[mask]

            if len(sub_st) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            ov = np.minimum(sub_en, we) - np.maximum(sub_st, ws)
            valid = ov > 0

            ov = ov[valid]
            sub_bin = sub_bin[valid]

            if len(ov) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            tmp = pd.DataFrame({"te_bin": sub_bin, "bp": ov})
            sums = (
                tmp.groupby("te_bin")["bp"]
                .sum()
                .reindex(cat_order, fill_value=0)
                .astype(float)
            )
            sums.name = row["enh_id"]
            rows.append(sums)

    out = pd.DataFrame(rows)
    out.index = enhancers["enh_id"].tolist()
    return out.reindex(columns=cat_order, fill_value=0)


def mean_local_bp_composition(enhancers, bg):
    per_enh_bp = neighborhood_bp_per_enhancer(enhancers, bg)
    per_enh_pct = per_enh_bp.div(per_enh_bp.sum(axis=1).replace(0, np.nan), axis=0) * 100
    return per_enh_pct.mean(axis=0).reindex(cat_order, fill_value=0)


def pooled_local_bp_by_category(enhancers, bg):
    per_enh_bp = neighborhood_bp_per_enhancer(enhancers, bg)
    return per_enh_bp.sum(axis=0).reindex(cat_order, fill_value=0).astype(float)


# =========================================================
# STATS
# =========================================================
def two_prop_ztest(count1, n1, count2, n2):
    if n1 <= 0 or n2 <= 0:
        return np.nan

    p1 = count1 / n1
    p2 = count2 / n2
    p_pool = (count1 + count2) / (n1 + n2)
    se = np.sqrt(p_pool * (1 - p_pool) * ((1 / n1) + (1 / n2)))

    if se == 0 or np.isnan(se):
        return np.nan

    z = (p1 - p2) / se
    return 2 * norm.sf(abs(z))


def category_bp_tests(obs_bp, local_bp, group_name):
    obs_total = obs_bp.sum()
    local_total = local_bp.sum()
    rows = []

    for cat in cat_order:
        p = two_prop_ztest(obs_bp[cat], obs_total, local_bp[cat], local_total)
        obs_prop = obs_bp[cat] / obs_total if obs_total > 0 else np.nan
        local_prop = local_bp[cat] / local_total if local_total > 0 else np.nan

        rows.append({
            "group": group_name,
            "te_class": cat,
            "obs_bp": obs_bp[cat],
            "local_bp": local_bp[cat],
            "obs_prop": obs_prop,
            "local_prop": local_prop,
            "log2_obs_over_local": np.log2((obs_prop + EPS) / (local_prop + EPS)),
            "p_value": p,
        })

    out = pd.DataFrame(rows)
    out["p_adj"] = bh_fdr(out["p_value"].values) if USE_FDR else out["p_value"]
    out["stars"] = out["p_adj"].map(p_to_star)
    return out


def per_class_bp_enrichment(obs_bp, local_bp, group_name):
    obs_bp = pd.Series(obs_bp, index=cat_order, dtype=float)
    local_bp = pd.Series(local_bp, index=cat_order, dtype=float)

    obs_total = obs_bp.sum()
    local_total = local_bp.sum()
    rows = []

    for cat in cat_order:
        obs_prop = obs_bp[cat] / obs_total if obs_total > 0 else np.nan
        local_prop = local_bp[cat] / local_total if local_total > 0 else np.nan
        p = two_prop_ztest(obs_bp[cat], obs_total, local_bp[cat], local_total)

        rows.append({
            "group": group_name,
            "te_class": cat,
            "obs_bp": obs_bp[cat],
            "obs_total_bp": obs_total,
            "local_bp": local_bp[cat],
            "local_total_bp": local_total,
            "obs_prop": obs_prop,
            "local_prop": local_prop,
            "log2_obs_over_local": np.log2((obs_prop + EPS) / (local_prop + EPS)),
            "p_value": p,
        })

    out = pd.DataFrame(rows)
    out["q_value"] = bh_fdr(out["p_value"].values)
    out["label"] = out["q_value"].map(p_to_star)
    return out


# =========================================================
# PLOTTING
# =========================================================
def plot_stacked_bar(ax, values, xpos, width=0.62):
    bottom = 0
    for val, color in zip(values, te_colors):
        ax.bar(xpos, val, width=width, bottom=bottom, color=color, edgecolor="black")
        if val > 0:
            if val < 5:
                ax.text(xpos + 0.34, bottom + val / 2, f"{val:.1f}%", ha="left", va="center", fontsize=15)
            else:
                ax.text(xpos, bottom + val / 2, f"{val:.1f}%", ha="center", va="center", fontsize=15)
        bottom += val


def make_grouped_te_legend(ax):
    box_w = 34
    box_text_sep = 6

    def legend_row(color, label):
        da = DrawingArea(box_w, 14, 0, 0)
        rect = patches.Rectangle((0, 2), 30, 10, facecolor=color, edgecolor="black", linewidth=1.2)
        da.add_artist(rect)
        txt = TextArea(label, textprops=dict(size=12, family="Arial"))
        return HPacker(children=[da, txt], align="center", pad=0, sep=box_text_sep)

    def header_row(label):
        spacer = DrawingArea(box_w + box_text_sep, 1, 0, 0)
        txt = TextArea(label, textprops=dict(size=14, weight="bold", family="Arial"))
        return HPacker(children=[spacer, txt], align="center", pad=0, sep=0)

    legend_box = VPacker(
        children=[
            header_row("RNA"),
            legend_row(te_colors[0], "LTR"),
            legend_row(te_colors[1], "LINE"),
            legend_row(te_colors[2], "SINE"),
            header_row("DNA"),
            legend_row(te_colors[3], "DNA-TIR"),
        ],
        align="left",
        pad=0,
        sep=4,
    )

    anchored_box = AnchoredOffsetbox(
        loc="upper left",
        child=legend_box,
        pad=0.3,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    ax.add_artist(anchored_box)


def make_plot(summary):
    fig, ax = plt.subplots(figsize=STACKED_FIGSIZE)
    ax.grid(False)

    all_values = [
        summary["shadow_observed_pct"].values,
        summary["shadow_local_pct"].values,
        summary["single_observed_pct"].values,
        summary["single_local_pct"].values,
        summary["genome_pct"].values,
    ]
    all_labels = ["shadows", "shadow\nneighborhood", "singles", "single\nneighborhood", "genome"]
    x_positions = [0, 1, 2.4, 3.4, 4.8]

    for xpos, vals in zip(x_positions, all_values):
        plot_stacked_bar(ax, vals, xpos, width=0.62)

    ax.set_ylabel("TE bp %")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(all_labels, fontsize=16)
    ax.set_ylim(0, 100)
    make_grouped_te_legend(ax)
    style_axes(ax)

    save_figure(
        fig,
        png_path=outdir / "mouse_stacked_TE_observed_local_genome_BP_based.png",
        pdf_path=outdir / "mouse_stacked_TE_observed_local_genome_BP_based.pdf",
        left=0.10,
        right=0.80,
        bottom=0.24,
        top=0.92,
    )


def plot_tissue_dotplot(enrichment_df, metric="log2_obs_over_local", outfile_png=None, outfile_pdf=None):
    fig, axes = plt.subplots(1, 2, figsize=TISSUE_SPLIT_FIGSIZE, sharey=True)

    group_order = ["shadow", "single"]
    group_titles = {"shadow": "Shadows", "single": "Singles"}
    yvals = np.arange(len(cat_order))
    tissue_yoffset = {"forebrain": -0.18, "heart": 0.00, "limb": 0.18, "pooled": 0.00}

    for ax, group in zip(axes, group_order):
        sub = enrichment_df[enrichment_df["group"] == group].copy()
        ax.axvline(0, color="black", linewidth=1)

        for tissue in tissue_order:
            tissue_sub = sub[sub["tissue"] == tissue].copy()
            for i, te in enumerate(cat_order):
                row = tissue_sub[tissue_sub["te_class"] == te]
                if row.empty:
                    continue

                x = float(row.iloc[0][metric])
                y = i + tissue_yoffset[tissue]
                ax.scatter(
                    x, y, s=55, color=tissue_colors[tissue],
                    edgecolor="black", linewidth=0.8, zorder=3
                )

                label = str(row.iloc[0]["label"])
                if label != "ns":
                    ax.text(x + 0.08, y, label, ha="left", va="center", fontsize=15, fontweight="bold")

        pooled_sub = sub[sub["tissue"] == "pooled"].copy()
        for i, te in enumerate(cat_order):
            row = pooled_sub[pooled_sub["te_class"] == te]
            if row.empty:
                continue

            x = float(row.iloc[0][metric])
            y = i + tissue_yoffset["pooled"]
            ax.scatter(
                x, y, s=75, color=tissue_colors["pooled"],
                edgecolor="black", linewidth=0.8, marker="D", zorder=4
            )

        ax.set_title(group_titles[group], fontsize=15)
        ax.set_xlim(*ENRICHMENT_YLIM)
        ax.set_yticks(yvals)
        ax.set_yticklabels(cat_order, fontsize=15)
        ax.invert_yaxis()
        style_axes(ax)

    axes[0].set_ylabel("TE class", fontsize=15)
    fig.supxlabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=16, y=0.08)

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["forebrain"], markeredgecolor="black", markersize=7, label="Forebrain"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["heart"], markeredgecolor="black", markersize=7, label="Heart"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["limb"], markeredgecolor="black", markersize=7, label="Limb"),
        Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Pooled"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=4, frameon=False, fontsize=16)

    save_figure(fig, png_path=outfile_png, pdf_path=outfile_pdf, left=0.10, right=0.98, bottom=0.24, top=0.82)


# =========================================================
# TISSUE ANALYSIS
# =========================================================
def run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df):
    shadow_df = shadow_df_all[shadow_df_all["tissue"] == tissue].copy()
    single_df = single_df_all[single_df_all["tissue"] == tissue].copy()

    if shadow_df.empty or single_df.empty:
        print(f"Skipping {tissue}: missing data in one group")
        return None

    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    shadow_obs_bp = observed_overlap_bp_by_category(shadow_df)
    single_obs_bp = observed_overlap_bp_by_category(single_df)
    shadow_local_bp = pooled_local_bp_by_category(shadow_enh, genome_df)
    single_local_bp = pooled_local_bp_by_category(single_enh, genome_df)

    shadow_stats = per_class_bp_enrichment(shadow_obs_bp, shadow_local_bp, group_name="shadow").copy()
    shadow_stats["tissue"] = tissue

    single_stats = per_class_bp_enrichment(single_obs_bp, single_local_bp, group_name="single").copy()
    single_stats["tissue"] = tissue

    return pd.concat([shadow_stats, single_stats], ignore_index=True)


def compute_tissue_enrichment_tables(shadow_df_all, single_df_all, genome_df):
    rows = []

    for tissue in tissue_order:
        tissue_stats = run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df)
        if tissue_stats is not None:
            rows.append(tissue_stats)

    shadow_enh_all = unique_enhancers(shadow_df_all)
    single_enh_all = unique_enhancers(single_df_all)
    shadow_obs_bp_all = observed_overlap_bp_by_category(shadow_df_all)
    single_obs_bp_all = observed_overlap_bp_by_category(single_df_all)
    shadow_local_bp_all = pooled_local_bp_by_category(shadow_enh_all, genome_df)
    single_local_bp_all = pooled_local_bp_by_category(single_enh_all, genome_df)

    shadow_stats_all = per_class_bp_enrichment(shadow_obs_bp_all, shadow_local_bp_all, group_name="shadow").copy()
    shadow_stats_all["tissue"] = "pooled"

    single_stats_all = per_class_bp_enrichment(single_obs_bp_all, single_local_bp_all, group_name="single").copy()
    single_stats_all["tissue"] = "pooled"

    rows.append(shadow_stats_all)
    rows.append(single_stats_all)

    out = pd.concat(rows, ignore_index=True)
    out["te_class"] = pd.Categorical(out["te_class"], categories=cat_order, ordered=True)
    out["tissue"] = pd.Categorical(out["tissue"], categories=tissue_order + ["pooled"], ordered=True)
    return out.sort_values(["group", "te_class", "tissue"])


# =========================================================
# MAIN
# =========================================================
def main():
    shadow_df = load_mouse_overlap(shadow_file)
    single_df = load_mouse_overlap(single_file)
    genome_df = load_background(genome_file)

    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    print("Unique shadow enhancers:", len(shadow_enh))
    print("Unique single enhancers:", len(single_enh))

    shadow_obs_pct = observed_overlap_composition(shadow_df)
    single_obs_pct = observed_overlap_composition(single_df)
    shadow_local_pct = mean_local_bp_composition(shadow_enh, genome_df)
    single_local_pct = mean_local_bp_composition(single_enh, genome_df)
    genome_pct = genome_composition(genome_df)

    shadow_obs_bp = observed_overlap_bp_by_category(shadow_df)
    single_obs_bp = observed_overlap_bp_by_category(single_df)
    shadow_local_bp = pooled_local_bp_by_category(shadow_enh, genome_df)
    single_local_bp = pooled_local_bp_by_category(single_enh, genome_df)
    genome_bp = genome_bp_by_category(genome_df)

    shadow_stats = category_bp_tests(shadow_obs_bp, shadow_local_bp, group_name="shadow")
    single_stats = category_bp_tests(single_obs_bp, single_local_bp, group_name="single")
    stats = pd.concat([shadow_stats, single_stats], ignore_index=True)

    summary = pd.DataFrame({
        "shadow_observed_pct": shadow_obs_pct,
        "shadow_local_pct": shadow_local_pct,
        "single_observed_pct": single_obs_pct,
        "single_local_pct": single_local_pct,
        "genome_pct": genome_pct,
        "shadow_observed_bp": shadow_obs_bp,
        "shadow_local_bp": shadow_local_bp,
        "single_observed_bp": single_obs_bp,
        "single_local_bp": single_local_bp,
        "genome_bp": genome_bp,
    }).reindex(cat_order)

    print("\nBP-based composition summary:")
    print(summary.round(4))
    print("\nObserved vs local BP tests:")
    print(stats.round(6))

    summary.to_csv(OUTPUT_FILES / "mouse_stacked_TE_BP_based_summary.tsv", sep="\t")
    stats.to_csv(OUTPUT_FILES / "mouse_TEclass_observed_vs_local_BP_ztests.tsv", sep="\t", index=False)

    tissue_enrichment = compute_tissue_enrichment_tables(shadow_df, single_df, genome_df)
    tissue_enrichment.to_csv(
        OUTPUT_FILES / "mouse_TEclass_observed_vs_local_BP_by_tissue.tsv",
        sep="\t",
        index=False,
    )

    make_plot(summary)
    plot_tissue_dotplot(
        tissue_enrichment,
        metric="log2_obs_over_local",
        outfile_png=OUTPUT_FIGURES / "mouse_TEclass_observed_vs_local_BP_dotplot_by_tissue.png",
        outfile_pdf=OUTPUT_FIGURES / "mouse_TEclass_observed_vs_local_BP_dotplot_by_tissue.pdf",
    )


if __name__ == "__main__":
    main()
#!/usr/bin/env python3


# =========================================================
# FILE INPUTS — MOUSE
# =========================================================
genome_file = str(INPUT_DIR / "FINAL_TE_genomefile_merged_dedup.bed")
shadow_file = str(INPUT_DIR / "Shadows_mouse_TEcooption_lastcol_Final.bed")
single_file = str(INPUT_DIR / "singlesmouse_TEcooption.bed")

outdir = OUTPUT_FIGURES


# =========================================================
# PARAMETERS
# =========================================================
FLANK_BP = 20000
USE_FDR = True
EPS = 1e-9

cat_order = ["LTR", "LINE", "SINE", "DNA-TIR"]
tissue_order = ["forebrain", "heart", "limb"]
te_colors = ["#c6dbef", "#6baed6", "#1f78b4", "#f16913"]
tissue_colors = {
    "forebrain": "#4C78A8",
    "heart": "#E45756",
    "limb": "#54A24B",
    "pooled": "black",
}

PANEL_HEIGHT = 3.0
STACKED_FIGSIZE = (6.5, 3.2)
TISSUE_SPLIT_FIGSIZE = (8.3, PANEL_HEIGHT)
COMBINED_DOTPLOT_FIGSIZE = (6.2, 3)
ENRICHMENT_YLIM = (-3, 3.2)

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# =========================================================
# HELPERS
# =========================================================
def map_te_bin(te_type):
    te_type = str(te_type)

    if te_type.startswith("LTR"):
        return "LTR"
    if te_type.startswith("LINE"):
        return "LINE"
    if te_type.startswith("SINE"):
        return "SINE"
    if te_type.startswith((
        "DNA/hAT",
        "DNA/TcMar",
        "DNA/PiggyBac",
        "DNA/MULE",
        "DNA/MuDR",
    )):
        return "DNA-TIR"
    return None


def assign_tissue_from_enh_id(enh_id):
    s = str(enh_id).lower()
    if "forebrain" in s or "fb" in s:
        return "forebrain"
    if "heart" in s:
        return "heart"
    if "limb" in s:
        return "limb"
    return "other"


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(len(pvals), np.nan)

    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return out

    p = pvals[valid]
    order = np.argsort(p)
    ranked = p[order]

    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    tmp = np.empty(len(q))
    tmp[order] = q
    out[valid] = tmp
    return out


def p_to_star(p):
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def style_axes(ax):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(
        axis="both",
        which="both",
        left=True,
        right=False,
        bottom=True,
        top=False,
        length=4,
        width=1.0,
        color="black",
        labelcolor="black",
    )


def save_figure(fig, png_path=None, pdf_path=None, left=0.10, right=0.82, bottom=0.22, top=0.92):
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    if png_path:
        fig.savefig(png_path, dpi=600)
    if pdf_path:
        fig.savefig(pdf_path)
    plt.close(fig)


# =========================================================
# LOADING
# =========================================================
def load_background(path):
    bg = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "te_name", "te_type"],
    )

    bg["start"] = pd.to_numeric(bg["start"], errors="coerce")
    bg["end"] = pd.to_numeric(bg["end"], errors="coerce")

    bg = bg.dropna(subset=["chrom", "start", "end", "te_type"]).copy()
    bg["bp"] = bg["end"] - bg["start"]
    bg = bg[bg["bp"] > 0].copy()

    bg["te_bin"] = bg["te_type"].map(map_te_bin)
    bg = bg[bg["te_bin"].notna()].copy()
    return bg


def load_mouse_overlap(path):
    raw = pd.read_csv(path, sep="\t", header=None)

    df = pd.DataFrame({
        "enh_chrom": raw.iloc[:, 0],
        "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "enh_id": raw.iloc[:, 3].astype(str),
        "te_type": raw.iloc[:, -2],
        "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
    })

    df = df.dropna(subset=[
        "enh_chrom",
        "enh_start",
        "enh_end",
        "enh_id",
        "te_type",
        "overlap_bp",
    ]).copy()

    df = df[df["te_type"] != "."].copy()
    df = df[df["overlap_bp"] > 0].copy()

    df["te_bin"] = df["te_type"].map(map_te_bin)
    df = df[df["te_bin"].notna()].copy()
    df["tissue"] = df["enh_id"].map(assign_tissue_from_enh_id)
    return df


def unique_enhancers(df):
    enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
    enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
    enh["win_end"] = enh["enh_end"] + FLANK_BP
    enh["tissue"] = enh["enh_id"].map(assign_tissue_from_enh_id)
    return enh


# =========================================================
# BP-BASED COMPOSITIONS
# =========================================================
def observed_overlap_bp_by_category(df):
    return (
        df.groupby("te_bin")["overlap_bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def observed_overlap_composition(df):
    x = observed_overlap_bp_by_category(df)
    return x / x.sum() * 100 if x.sum() > 0 else x


def composition_from_counts(counts):
    counts = pd.Series(counts, index=cat_order, dtype=float)
    total = counts.sum()
    if total <= 0:
        return pd.Series(np.nan, index=cat_order)
    return counts / total * 100


def genome_bp_by_category(bg):
    return (
        bg.groupby("te_bin")["bp"]
        .sum()
        .reindex(cat_order, fill_value=0)
        .astype(float)
    )


def genome_composition(bg):
    x = genome_bp_by_category(bg)
    return x / x.sum() * 100 if x.sum() > 0 else x


def neighborhood_bp_per_enhancer(enhancers, bg):
    rows = []

    for chrom, enh_chr in enhancers.groupby("enh_chrom"):
        te_chr = bg[bg["chrom"] == chrom].copy()

        te_st = te_chr["start"].to_numpy()
        te_en = te_chr["end"].to_numpy()
        te_bin = te_chr["te_bin"].to_numpy()

        for _, row in enh_chr.iterrows():
            ws = row["win_start"]
            we = row["win_end"]

            mask = (te_en > ws) & (te_st < we)
            sub_st = te_st[mask]
            sub_en = te_en[mask]
            sub_bin = te_bin[mask]

            if len(sub_st) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            ov = np.minimum(sub_en, we) - np.maximum(sub_st, ws)
            valid = ov > 0

            ov = ov[valid]
            sub_bin = sub_bin[valid]

            if len(ov) == 0:
                rows.append(pd.Series(0.0, index=cat_order, name=row["enh_id"]))
                continue

            tmp = pd.DataFrame({"te_bin": sub_bin, "bp": ov})
            sums = (
                tmp.groupby("te_bin")["bp"]
                .sum()
                .reindex(cat_order, fill_value=0)
                .astype(float)
            )
            sums.name = row["enh_id"]
            rows.append(sums)

    out = pd.DataFrame(rows)
    out.index = enhancers["enh_id"].tolist()
    return out.reindex(columns=cat_order, fill_value=0)


def mean_local_bp_composition(enhancers, bg):
    per_enh_bp = neighborhood_bp_per_enhancer(enhancers, bg)
    per_enh_pct = per_enh_bp.div(per_enh_bp.sum(axis=1).replace(0, np.nan), axis=0) * 100
    return per_enh_pct.mean(axis=0).reindex(cat_order, fill_value=0)


def pooled_local_bp_by_category(enhancers, bg):
    per_enh_bp = neighborhood_bp_per_enhancer(enhancers, bg)
    return per_enh_bp.sum(axis=0).reindex(cat_order, fill_value=0).astype(float)


# =========================================================
# STATS
# =========================================================
def two_prop_ztest(count1, n1, count2, n2):
    if n1 <= 0 or n2 <= 0:
        return np.nan

    p1 = count1 / n1
    p2 = count2 / n2
    p_pool = (count1 + count2) / (n1 + n2)
    se = np.sqrt(p_pool * (1 - p_pool) * ((1 / n1) + (1 / n2)))

    if se == 0 or np.isnan(se):
        return np.nan

    z = (p1 - p2) / se
    return 2 * norm.sf(abs(z))


def category_bp_tests(obs_bp, local_bp, group_name):
    obs_total = obs_bp.sum()
    local_total = local_bp.sum()
    rows = []

    for cat in cat_order:
        p = two_prop_ztest(obs_bp[cat], obs_total, local_bp[cat], local_total)
        obs_prop = obs_bp[cat] / obs_total if obs_total > 0 else np.nan
        local_prop = local_bp[cat] / local_total if local_total > 0 else np.nan

        rows.append({
            "group": group_name,
            "te_class": cat,
            "obs_bp": obs_bp[cat],
            "local_bp": local_bp[cat],
            "obs_prop": obs_prop,
            "local_prop": local_prop,
            "log2_obs_over_local": np.log2((obs_prop + EPS) / (local_prop + EPS)),
            "p_value": p,
        })

    out = pd.DataFrame(rows)
    out["p_adj"] = bh_fdr(out["p_value"].values) if USE_FDR else out["p_value"]
    out["stars"] = out["p_adj"].map(p_to_star)
    return out


def per_class_bp_enrichment(obs_bp, local_bp, group_name):
    obs_bp = pd.Series(obs_bp, index=cat_order, dtype=float)
    local_bp = pd.Series(local_bp, index=cat_order, dtype=float)

    obs_total = obs_bp.sum()
    local_total = local_bp.sum()
    rows = []

    for cat in cat_order:
        obs_prop = obs_bp[cat] / obs_total if obs_total > 0 else np.nan
        local_prop = local_bp[cat] / local_total if local_total > 0 else np.nan
        p = two_prop_ztest(obs_bp[cat], obs_total, local_bp[cat], local_total)

        rows.append({
            "group": group_name,
            "te_class": cat,
            "obs_bp": obs_bp[cat],
            "obs_total_bp": obs_total,
            "local_bp": local_bp[cat],
            "local_total_bp": local_total,
            "obs_prop": obs_prop,
            "local_prop": local_prop,
            "log2_obs_over_local": np.log2((obs_prop + EPS) / (local_prop + EPS)),
            "p_value": p,
        })

    out = pd.DataFrame(rows)
    out["q_value"] = bh_fdr(out["p_value"].values)
    out["label"] = out["q_value"].map(p_to_star)
    return out


# =========================================================
# PLOTTING
# =========================================================
def plot_stacked_bar(ax, values, xpos, width=0.62):
    bottom = 0
    for val, color in zip(values, te_colors):
        ax.bar(xpos, val, width=width, bottom=bottom, color=color, edgecolor="black")
        if val > 0:
            if val < 5:
                ax.text(xpos + 0.34, bottom + val / 2, f"{val:.1f}%", ha="left", va="center", fontsize=16)
            else:
                ax.text(xpos, bottom + val / 2, f"{val:.1f}%", ha="center", va="center", fontsize=16)
        bottom += val


def make_grouped_te_legend(ax):
    box_w = 34
    box_text_sep = 6

    def legend_row(color, label):
        da = DrawingArea(box_w, 14, 0, 0)
        rect = patches.Rectangle((0, 2), 30, 10, facecolor=color, edgecolor="black", linewidth=1.2)
        da.add_artist(rect)
        txt = TextArea(label, textprops=dict(size=12, family="Arial"))
        return HPacker(children=[da, txt], align="center", pad=0, sep=box_text_sep)

    def header_row(label):
        spacer = DrawingArea(box_w + box_text_sep, 1, 0, 0)
        txt = TextArea(label, textprops=dict(size=14, weight="bold", family="Arial"))
        return HPacker(children=[spacer, txt], align="center", pad=0, sep=0)

    legend_box = VPacker(
        children=[
            header_row("RNA"),
            legend_row(te_colors[0], "LTR"),
            legend_row(te_colors[1], "LINE"),
            legend_row(te_colors[2], "SINE"),
            header_row("DNA"),
            legend_row(te_colors[3], "DNA-TIR"),
        ],
        align="left",
        pad=0,
        sep=4,
    )

    anchored_box = AnchoredOffsetbox(
        loc="upper left",
        child=legend_box,
        pad=0.3,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    ax.add_artist(anchored_box)


def make_plot(summary):
    fig, ax = plt.subplots(figsize=STACKED_FIGSIZE)
    ax.grid(False)

    all_values = [
        summary["shadow_observed_pct"].values,
        summary["shadow_local_pct"].values,
        summary["single_observed_pct"].values,
        summary["single_local_pct"].values,
        summary["genome_pct"].values,
    ]
    all_labels = ["shadows", "shadow\nneighborhood", "singles", "single\nneighborhood", "genome"]
    x_positions = [0, 1, 2.4, 3.4, 4.8]

    for xpos, vals in zip(x_positions, all_values):
        plot_stacked_bar(ax, vals, xpos, width=0.62)

    ax.set_ylabel("TE bp %")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(all_labels, fontsize=12)
    ax.set_ylim(0, 100)
    make_grouped_te_legend(ax)
    style_axes(ax)

    save_figure(
        fig,
        png_path=outdir / "mouse_stacked_TE_observed_local_genome_BP_based.png",
        pdf_path=outdir / "mouse_stacked_TE_observed_local_genome_BP_based.pdf",
        left=0.10,
        right=0.80,
        bottom=0.24,
        top=0.92,
    )


def plot_tissue_dotplot(enrichment_df, metric="log2_obs_over_local", outfile_png=None, outfile_pdf=None):
    fig, axes = plt.subplots(1, 2, figsize=TISSUE_SPLIT_FIGSIZE, sharey=True)

    group_order = ["shadow", "single"]
    group_titles = {"shadow": "Shadows", "single": "Singles"}
    yvals = np.arange(len(cat_order))
    tissue_yoffset = {"forebrain": -0.18, "heart": 0.00, "limb": 0.18, "pooled": 0.00}

    for ax, group in zip(axes, group_order):
        sub = enrichment_df[enrichment_df["group"] == group].copy()
        ax.axvline(0, color="black", linewidth=1)

        for tissue in tissue_order:
            tissue_sub = sub[sub["tissue"] == tissue].copy()
            for i, te in enumerate(cat_order):
                row = tissue_sub[tissue_sub["te_class"] == te]
                if row.empty:
                    continue

                x = float(row.iloc[0][metric])
                y = i + tissue_yoffset[tissue]
                ax.scatter(
                    x, y, s=55, color=tissue_colors[tissue],
                    edgecolor="black", linewidth=0.8, zorder=3
                )

                label = str(row.iloc[0]["label"])
                if label != "ns":
                    ax.text(x + 0.08, y, label, ha="left", va="center", fontsize=9, fontweight="bold")

        pooled_sub = sub[sub["tissue"] == "pooled"].copy()
        for i, te in enumerate(cat_order):
            row = pooled_sub[pooled_sub["te_class"] == te]
            if row.empty:
                continue

            x = float(row.iloc[0][metric])
            y = i + tissue_yoffset["pooled"]
            ax.scatter(
                x, y, s=75, color=tissue_colors["pooled"],
                edgecolor="black", linewidth=0.8, marker="D", zorder=4
            )

        ax.set_title(group_titles[group], fontsize=14)
        ax.set_xlim(*ENRICHMENT_YLIM)
        ax.set_yticks(yvals)
        ax.set_yticklabels(cat_order, fontsize=12)
        ax.invert_yaxis()
        style_axes(ax)

    axes[0].set_ylabel("TE class", fontsize=12)
    fig.supxlabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=13, y=0.08)

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["forebrain"], markeredgecolor="black", markersize=7, label="Forebrain"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["heart"], markeredgecolor="black", markersize=7, label="Heart"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["limb"], markeredgecolor="black", markersize=7, label="Limb"),
        Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Pooled"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=4, frameon=False, fontsize=11)

    save_figure(fig, png_path=outfile_png, pdf_path=outfile_pdf, left=0.10, right=0.98, bottom=0.24, top=0.82)


def plot_combined_vertical_dotplot_clean(
    enrichment_df,
    metric="log2_obs_over_local",
    outfile_png=None,
    outfile_pdf=None,
    show_pooled=True,
    show_stars=True,
):
    fig, ax = plt.subplots(figsize=COMBINED_DOTPLOT_FIGSIZE)

    x_base = np.arange(len(cat_order)) * 1.45
    tissue_xoffset = {
        "forebrain": -0.30,
        "heart": 0.00,
        "limb": 0.30,
        "pooled": 0.55,
    }
    group_xoffset = {
        "shadow": -0.05,
        "single": 0.05,
    }
    star_offset = 0.16

    ax.axhline(0, color="black", linewidth=1.1, zorder=1)
    for x_sep in (x_base[:-1] + x_base[1:]) / 2:
        ax.axvline(x_sep, color="0.90", linewidth=0.8, zorder=0)

    for tissue in tissue_order:
        color = tissue_colors[tissue]

        for te_i, te in enumerate(cat_order):
            sh = enrichment_df[
                (enrichment_df["tissue"] == tissue) &
                (enrichment_df["group"] == "shadow") &
                (enrichment_df["te_class"] == te)
            ]
            si = enrichment_df[
                (enrichment_df["tissue"] == tissue) &
                (enrichment_df["group"] == "single") &
                (enrichment_df["te_class"] == te)
            ]

            if sh.empty or si.empty:
                continue

            x_sh = x_base[te_i] + tissue_xoffset[tissue] + group_xoffset["shadow"]
            x_si = x_base[te_i] + tissue_xoffset[tissue] + group_xoffset["single"]
            y_sh = float(sh.iloc[0][metric])
            y_si = float(si.iloc[0][metric])

            ax.scatter(
                x_sh,
                y_sh,
                s=58,
                facecolor=color,
                edgecolor="black",
                linewidth=0.9,
                zorder=4,
            )
            ax.scatter(
                x_si,
                y_si,
                s=58,
                facecolor="white",
                edgecolor=color,
                linewidth=1.8,
                zorder=4,
            )

            if show_stars:
                for x_pt, y_pt, row in [
                    (x_sh, y_sh, sh.iloc[0]),
                    (x_si, y_si, si.iloc[0]),
                ]:
                    label = str(row["label"])
                    if label != "ns":
                        ax.text(
                            x_pt,
                            y_pt + star_offset,
                            label,
                            ha="center",
                            va="bottom",
                            fontsize=8,
                            fontweight="bold",
                            zorder=6,
                            clip_on=False,
                        )

    if show_pooled:
        for te_i, te in enumerate(cat_order):
            for group in ["shadow", "single"]:
                row = enrichment_df[
                    (enrichment_df["tissue"] == "pooled") &
                    (enrichment_df["group"] == group) &
                    (enrichment_df["te_class"] == te)
                ]
                if row.empty:
                    continue

                x = x_base[te_i] + tissue_xoffset["pooled"] + group_xoffset[group]
                y = float(row.iloc[0][metric])

                if group == "shadow":
                    ax.scatter(x, y, s=84, marker="D", facecolor="black", edgecolor="black", linewidth=1.0, zorder=5)
                else:
                    ax.scatter(x, y, s=84, marker="D", facecolor="white", edgecolor="black", linewidth=1.8, zorder=5)

                if show_stars:
                    label = str(row.iloc[0]["label"])
                    if label != "ns":
                        ax.text(
                            x,
                            y + star_offset,
                            label,
                            ha="center",
                            va="bottom",
                            fontsize=8,
                            fontweight="bold",
                            zorder=6,
                            clip_on=False,
                        )

    ax.set_xticks(x_base)
    ax.set_xticklabels(cat_order, fontsize=13)
    ax.set_xlim(x_base[0] - 0.65, x_base[-1] + 0.90)
    ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=13)
    ax.set_ylim(*ENRICHMENT_YLIM)
    style_axes(ax)

    tissue_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["forebrain"], markeredgecolor="black", markersize=7, label="Forebrain"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["heart"], markeredgecolor="black", markersize=7, label="Heart"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["limb"], markeredgecolor="black", markersize=7, label="Limb"),
    ]
    if show_pooled:
        tissue_handles.append(
            Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Pooled")
        )

    class_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Shadows"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="white", markeredgecolor="black", markersize=7, label="Singles"),
    ]

    leg1 = ax.legend(
        handles=tissue_handles,
        title="Tissue",
        frameon=False,
        fontsize=14,
        title_fontsize=11,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
    )
    ax.add_artist(leg1)

    ax.legend(
        handles=class_handles,
        title="Enhancer class",
        frameon=False,
        fontsize=14,
        title_fontsize=11,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.50),
    )

    save_figure(fig, png_path=outfile_png, pdf_path=outfile_pdf, left=0.10, right=0.77, bottom=0.24, top=0.92)


# =========================================================
# TISSUE ANALYSIS
# =========================================================
def run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df):
    shadow_df = shadow_df_all[shadow_df_all["tissue"] == tissue].copy()
    single_df = single_df_all[single_df_all["tissue"] == tissue].copy()

    if shadow_df.empty or single_df.empty:
        print(f"Skipping {tissue}: missing data in one group")
        return None

    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    shadow_obs_bp = observed_overlap_bp_by_category(shadow_df)
    single_obs_bp = observed_overlap_bp_by_category(single_df)
    shadow_local_bp = pooled_local_bp_by_category(shadow_enh, genome_df)
    single_local_bp = pooled_local_bp_by_category(single_enh, genome_df)

    shadow_stats = per_class_bp_enrichment(shadow_obs_bp, shadow_local_bp, group_name="shadow").copy()
    shadow_stats["tissue"] = tissue

    single_stats = per_class_bp_enrichment(single_obs_bp, single_local_bp, group_name="single").copy()
    single_stats["tissue"] = tissue

    return pd.concat([shadow_stats, single_stats], ignore_index=True)


def compute_tissue_enrichment_tables(shadow_df_all, single_df_all, genome_df):
    rows = []

    for tissue in tissue_order:
        tissue_stats = run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df)
        if tissue_stats is not None:
            rows.append(tissue_stats)

    shadow_enh_all = unique_enhancers(shadow_df_all)
    single_enh_all = unique_enhancers(single_df_all)
    shadow_obs_bp_all = observed_overlap_bp_by_category(shadow_df_all)
    single_obs_bp_all = observed_overlap_bp_by_category(single_df_all)
    shadow_local_bp_all = pooled_local_bp_by_category(shadow_enh_all, genome_df)
    single_local_bp_all = pooled_local_bp_by_category(single_enh_all, genome_df)

    shadow_stats_all = per_class_bp_enrichment(shadow_obs_bp_all, shadow_local_bp_all, group_name="shadow").copy()
    shadow_stats_all["tissue"] = "pooled"

    single_stats_all = per_class_bp_enrichment(single_obs_bp_all, single_local_bp_all, group_name="single").copy()
    single_stats_all["tissue"] = "pooled"

    rows.append(shadow_stats_all)
    rows.append(single_stats_all)

    out = pd.concat(rows, ignore_index=True)
    out["te_class"] = pd.Categorical(out["te_class"], categories=cat_order, ordered=True)
    out["tissue"] = pd.Categorical(out["tissue"], categories=tissue_order + ["pooled"], ordered=True)
    return out.sort_values(["group", "te_class", "tissue"])


# =========================================================
# MAIN
# =========================================================
def main():
    shadow_df = load_mouse_overlap(shadow_file)
    single_df = load_mouse_overlap(single_file)
    genome_df = load_background(genome_file)

    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    print("Unique shadow enhancers:", len(shadow_enh))
    print("Unique single enhancers:", len(single_enh))

    shadow_obs_pct = observed_overlap_composition(shadow_df)
    single_obs_pct = observed_overlap_composition(single_df)
    shadow_local_pct = mean_local_bp_composition(shadow_enh, genome_df)
    single_local_pct = mean_local_bp_composition(single_enh, genome_df)
    genome_pct = genome_composition(genome_df)

    shadow_obs_bp = observed_overlap_bp_by_category(shadow_df)
    single_obs_bp = observed_overlap_bp_by_category(single_df)
    shadow_local_bp = pooled_local_bp_by_category(shadow_enh, genome_df)
    single_local_bp = pooled_local_bp_by_category(single_enh, genome_df)
    genome_bp = genome_bp_by_category(genome_df)

    shadow_stats = category_bp_tests(shadow_obs_bp, shadow_local_bp, group_name="shadow")
    single_stats = category_bp_tests(single_obs_bp, single_local_bp, group_name="single")
    stats = pd.concat([shadow_stats, single_stats], ignore_index=True)

    summary = pd.DataFrame({
        "shadow_observed_pct": shadow_obs_pct,
        "shadow_local_pct": shadow_local_pct,
        "single_observed_pct": single_obs_pct,
        "single_local_pct": single_local_pct,
        "genome_pct": genome_pct,
        "shadow_observed_bp": shadow_obs_bp,
        "shadow_local_bp": shadow_local_bp,
        "single_observed_bp": single_obs_bp,
        "single_local_bp": single_local_bp,
        "genome_bp": genome_bp,
    }).reindex(cat_order)

    print("\nBP-based composition summary:")
    print(summary.round(4))
    print("\nObserved vs local BP tests:")
    print(stats.round(6))

    summary.to_csv(OUTPUT_FILES / "mouse_stacked_TE_BP_based_summary.tsv", sep="\t")
    stats.to_csv(OUTPUT_FILES / "mouse_TEclass_observed_vs_local_BP_ztests.tsv", sep="\t", index=False)

    tissue_enrichment = compute_tissue_enrichment_tables(shadow_df, single_df, genome_df)
    tissue_enrichment.to_csv(
        OUTPUT_FILES / "mouse_TEclass_observed_vs_local_BP_by_tissue.tsv",
        sep="\t",
        index=False,
    )

    make_plot(summary)
    plot_combined_vertical_dotplot_clean(
        tissue_enrichment,
        metric="log2_obs_over_local",
        outfile_png=OUTPUT_FIGURES / "combined_vertical_dotplot_TEclass_enrichment_clean_BP.png",
        outfile_pdf=OUTPUT_FIGURES / "combined_vertical_dotplot_TEclass_enrichment_clean_BP.pdf",
        show_pooled=True,
        show_stars=True,
    )


if __name__ == "__main__":
    main()
