from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch

# ----------------------------
# file inputs
# ----------------------------
ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"
OUTPUT = ROOT / "output"
PLOT_FIGURES = ROOT / "plot_figures"
OUTPUT.mkdir(parents=True, exist_ok=True)
PLOT_FIGURES.mkdir(parents=True, exist_ok=True)

# legacy external path examples removed
single_file = str(INPUT / "Full_unfiltered_overlap_singles_merged.bed")
shadow_file = str(INPUT / "Final_overlap_shadows_TEtype.bed")
genome_file = str(INPUT / "FINAL_TE_genomefile_merged_dedup.bed")

# legacy external path examples removed
# ----------------------------
# parameters
# ----------------------------
MIN_OVERLAP_BP = 50
FLANK_BP = 50000

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "DejaVu Sans"
})

cat_order = ["LTR", "LINE", "SINE", "DNA/RC", "DNA/other"]
te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#fdd49e', '#f16913']

# ----------------------------
# helpers
# ----------------------------
def map_te_bin(te_type):
    te_type = str(te_type)
    if te_type.startswith("LTR"):
        return "LTR"
    if te_type.startswith("LINE"):
        return "LINE"
    if te_type.startswith("SINE"):
        return "SINE"
    if te_type == "DNA/RC":
        return "DNA/RC"
    if te_type.startswith("DNA"):
        return "DNA/other"
    return None

def load_background(path):
    bg = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "te_type"]
    )
    bg["start"] = pd.to_numeric(bg["start"], errors="coerce")
    bg["end"] = pd.to_numeric(bg["end"], errors="coerce")
    bg = bg.dropna(subset=["start", "end", "te_type"]).copy()
    bg["bp"] = bg["end"] - bg["start"]
    bg = bg[bg["bp"] > 0].copy()
    bg["te_bin"] = bg["te_type"].map(map_te_bin)
    bg = bg[bg["te_bin"].notna()].copy()
    return bg

def load_single_overlap(path):
    # 9-column file
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["enh_chrom", "enh_start", "enh_end", "enh_id",
               "te_chrom", "te_start", "te_end", "te_type", "overlap_bp"]
    )
    df["enh_start"] = pd.to_numeric(df["enh_start"], errors="coerce")
    df["enh_end"] = pd.to_numeric(df["enh_end"], errors="coerce")
    df["overlap_bp"] = pd.to_numeric(df["overlap_bp"], errors="coerce")
    df = df.dropna(subset=["enh_start", "enh_end", "te_type", "overlap_bp"]).copy()
    df = df[(df["overlap_bp"] >= MIN_OVERLAP_BP) & (df["te_type"] != ".")].copy()
    df["te_bin"] = df["te_type"].map(map_te_bin)
    df = df[df["te_bin"].notna()].copy()
    return df

def load_shadow_overlap(path):
    # enhancer in cols 0-3, TE type in col 8, overlap in last col
    raw = pd.read_csv(path, sep="\t", header=None)
    df = pd.DataFrame({
        "enh_chrom": raw.iloc[:, 0],
        "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "enh_id": raw.iloc[:, 3].astype(str),
        "te_type": raw.iloc[:, 8],
        "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
    })
    df = df.dropna(subset=["enh_start", "enh_end", "te_type", "overlap_bp"]).copy()
    df = df[(df["overlap_bp"] >= MIN_OVERLAP_BP) & (df["te_type"] != ".")].copy()
    df["te_bin"] = df["te_type"].map(map_te_bin)
    df = df[df["te_bin"].notna()].copy()
    return df

def unique_enhancers(df):
    enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
    enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
    enh["win_end"] = enh["enh_end"] + FLANK_BP
    return enh

# ----------------------------
# FIXED observed composition:
# count each enhancer once per TE category
# ----------------------------
def observed_enhancer_category_counts(df):
    pairs = df[["enh_id", "te_bin"]].drop_duplicates().copy()
    counts = pairs["te_bin"].value_counts().reindex(cat_order, fill_value=0).astype(float)
    return counts

def observed_enhancer_category_composition(df):
    counts = observed_enhancer_category_counts(df)
    total = counts.sum()
    return (counts / total * 100) if total > 0 else counts

def genome_composition(bg):
    x = bg.groupby("te_bin")["bp"].sum().reindex(cat_order, fill_value=0)
    return x / x.sum() * 100

# ----------------------------
# FIXED neighborhood:
# for each enhancer, a TE category is present once or absent
# repeated same-category elements do not count multiple times
# ----------------------------
def neighborhood_presence_per_enhancer(enhancers, bg):
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

            counts = pd.Series(0.0, index=cat_order, name=row["enh_id"])

            if len(sub_st) == 0:
                rows.append(counts)
                continue

            ov = np.minimum(sub_en, we) - np.maximum(sub_st, ws)
            valid = ov > 0
            present_bins = pd.Series(sub_bin[valid]).drop_duplicates()

            present_bins = present_bins[present_bins.isin(cat_order)]
            counts.loc[present_bins] = 1.0
            rows.append(counts)

    out = pd.DataFrame(rows)
    out.index = enhancers["enh_id"].tolist()
    return out

def mean_local_composition(enhancers, bg):
    per_enh = neighborhood_presence_per_enhancer(enhancers, bg)
    mean_counts = per_enh.mean(axis=0).reindex(cat_order, fill_value=0)
    total = mean_counts.sum()
    return (mean_counts / total * 100) if total > 0 else mean_counts

def plot_stacked_bar(ax, values, xpos, width=0.6):
    bottom = 0
    for val, color in zip(values, te_colors):
        ax.bar(xpos, val, width=width, bottom=bottom, color=color, edgecolor="black")
        if val > 0:
            if val < 5:
                ax.text(xpos + 0.33, bottom + val / 2, f"{val:.1f}%", ha="left", va="center", fontsize=10)
            else:
                ax.text(xpos, bottom + val / 2, f"{val:.1f}%", ha="center", va="center", fontsize=10)
        bottom += val

# ----------------------------
# load data
# ----------------------------
shadow_df = load_shadow_overlap(shadow_file)
single_df = load_single_overlap(single_file)
genome_df = load_background(genome_file)

shadow_enh = unique_enhancers(shadow_df)
single_enh = unique_enhancers(single_df)

# ----------------------------
# compute compositions
# ----------------------------
shadow_obs = observed_enhancer_category_composition(shadow_df)
shadow_local = mean_local_composition(shadow_enh, genome_df)
single_obs = observed_enhancer_category_composition(single_df)
single_local = mean_local_composition(single_enh, genome_df)
genome_comp = genome_composition(genome_df)

summary = pd.DataFrame({
    "shadow_observed": shadow_obs,
    "shadow_local_100bp": shadow_local,
    "single_observed": single_obs,
    "single_local_100bp": single_local,
    "genome": genome_comp
}).reindex(cat_order)

print("Unique shadow enhancers:", len(shadow_enh))
print("Unique single enhancers:", len(single_enh))
print("\nObserved raw counts (enhancer counted once per TE category):")
print(pd.DataFrame({
    "shadow_counts": observed_enhancer_category_counts(shadow_df).astype(int),
    "single_counts": observed_enhancer_category_counts(single_df).astype(int)
}).reindex(cat_order))

print("\nComposition table (%):")
print(summary.round(2))

summary.to_csv(OUTPUT / "stackedplot_local100bp_summary.tsv", sep="\t")

# ----------------------------
# plot
# ----------------------------
fig, ax = plt.subplots(figsize=(8.5, 4.2))
ax.grid(False)

all_values = [
    shadow_obs.values,
    shadow_local.values,
    single_obs.values,
    single_local.values,
    genome_comp.values
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

ax.set_ylabel("percentage")
ax.set_xticks(x_positions)
ax.set_xticklabels(all_labels, fontsize=12)
ax.set_ylim(0, 100)

handles = [
    Patch(facecolor=te_colors[0], edgecolor='black', label='LTR'),
    Patch(facecolor=te_colors[1], edgecolor='black', label='LINE'),
    Patch(facecolor=te_colors[2], edgecolor='black', label='SINE'),
    Patch(facecolor=te_colors[3], edgecolor='black', label='DNA/RC'),
    Patch(facecolor=te_colors[4], edgecolor='black', label='DNA/other')
]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1), frameon=True)

for spine in ax.spines.values():
    spine.set_visible(True)

plt.title("Observed TE-category usage and local TE neighborhood composition (±100 bp)")
plt.tight_layout()
plt.savefig(PLOT_FIGURES / "stacked_TE_observed_local_genome_100bp.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.show()
