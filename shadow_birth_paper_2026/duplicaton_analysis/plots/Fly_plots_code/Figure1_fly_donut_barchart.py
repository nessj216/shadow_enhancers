
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

PAPER_ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUT_DIR = HERE / "input"
OUTPUT_DIR = PAPER_ROOT / "plots" / "output_pngs" / "fly"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ====== MASTER FONT SETTINGS ======
mpl.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# ====== Input files ======
BREAKDOWN_FILE = str(INPUT_DIR / "FINAL_breakdown_single_double_flanks_UPDATED_1e4_flankmethod.csv")
TOTAL_BED_FILE = str(PAPER_ROOT / "duplicaton_analysis" / "enhancer_hits_per_shadow_bin" / "Fly" / "data" / "011925_all_shadowsets_DM6.bed")
csv_file = BREAKDOWN_FILE
OUT_PNG = str(OUTPUT_DIR / "pie_left_nested_donut.png")

FINAL_BREAKDOWN_CSV = BREAKDOWN_FILE
SHADOW_BED = Path(TOTAL_BED_FILE)
# ====== Appearance ======
C_ANY    = "#1F4E99"
C_NONE   = "#E5E7E9"
C_ENH    = "#F39C12"
C_FLANK  = "#56B4E9"
TRANSPARENT = (1, 1, 1, 0)

SHOW_TEXT = True
DPI = 600


# ------------------------------------------------------------------
# Read files and compute counts
# ------------------------------------------------------------------

# final breakdown file expected columns:
# gene_name, pair, source, flank_hit_class
df = pd.read_csv(BREAKDOWN_FILE)

# normalize column names just in case
df.columns = [c.strip() for c in df.columns]

required_cols = {"gene_name", "source"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in breakdown file: {missing}")

# BED file: total universe comes from column 4 (0-based index 3)
bed = pd.read_csv(TOTAL_BED_FILE, sep="\t", header=None)

if bed.shape[1] < 4:
    raise ValueError("BED file does not have at least 4 columns.")

# Unique gene sets from BED
total_gene_sets = set(
    bed.iloc[:, 3]
    .dropna()
    .astype(str)
    .str.strip()
)

TOTAL = len(total_gene_sets)

# Unique gene sets in each source group
enh_gene_sets = set(
    df.loc[df["source"] == "enhancer_hit", "gene_name"]
    .dropna()
    .astype(str)
    .str.strip()
)

flank_gene_sets = set(
    df.loc[df["source"] == "flank", "gene_name"]
    .dropna()
    .astype(str)
    .str.strip()
)

# Restrict everything to the BED universe
enh_gene_sets_in_bed = enh_gene_sets & total_gene_sets
flank_gene_sets_in_bed = flank_gene_sets & total_gene_sets

any_gene_sets_in_bed = enh_gene_sets_in_bed | flank_gene_sets_in_bed

flank_only_gene_sets_in_bed = flank_gene_sets_in_bed - enh_gene_sets_in_bed

TOTAL = len(total_gene_sets)
ANY = len(any_gene_sets_in_bed)
NONE = len(total_gene_sets - any_gene_sets_in_bed)

ENH = len(enh_gene_sets_in_bed)
FLANK = len(flank_only_gene_sets_in_bed)

print(f"TOTAL BED gene sets: {TOTAL}")
print(f"Any duplicated gene sets in BED: {ANY}")
print(f"No-hit BED gene sets: {NONE}")
print(f"Enhancer-hit gene sets in BED: {ENH}")
print(f"Flank-only gene sets in BED: {FLANK}")

if ANY != ENH + FLANK:
    print("Warning: ANY != ENH + FLANK")


# ------------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------------

def full_pie(ax, values, colors, center_pct=None, center_subtitle=None, startangle=90):
    wedges, _ = ax.pie(
        values,
        colors=colors,
        startangle=startangle,
        wedgeprops=dict(edgecolor="black", linewidth=0.3)
    )
    ax.set(aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])

    if SHOW_TEXT and center_pct is not None:
        ax.text(
            0, 0.08, f"{center_pct:.1f}%",
            ha="center", va="center",
            fontweight="bold"
        )
        if center_subtitle:
            ax.text(
                0, -0.15, center_subtitle,
                ha="center", va="center"
            )

    for spine in ax.spines.values():
        spine.set_visible(False)
    return wedges


def nested_donut_breakdown_of_any(ax, startangle=90):
    outer_radius = 1.0
    outer_width  = 0.32
    inner_radius = 0.68
    inner_width  = 0.28

    # Outer ring: any hit vs no hit
    ax.pie(
        [ANY, NONE],
        colors=[C_ANY, C_NONE],
        startangle=startangle,
        radius=outer_radius,
        wedgeprops=dict(width=outer_width, edgecolor="black", linewidth=0.3)
    )

    # Inner ring: enhancer-hit vs flank-only vs transparent filler
    ax.pie(
        [ENH, FLANK, NONE],
        colors=[C_ENH, C_FLANK, TRANSPARENT],
        startangle=startangle,
        radius=inner_radius,
        wedgeprops=dict(width=inner_width, edgecolor="black", linewidth=0.3)
    )

    ax.set(aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


# ------------------------------------------------------------------
# Build figure
# ------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(5, 3), dpi=DPI)

pct_any = 100 * ANY / TOTAL if TOTAL > 0 else 0

full_pie(
    axes[0],
    values=[ANY, NONE],
    colors=[C_ANY, C_NONE],
    center_pct=pct_any,
    center_subtitle="hit sets"
)

nested_donut_breakdown_of_any(axes[1])

if SHOW_TEXT:
    right_handles = [
        Line2D(
            [0], [0],
            marker='o', color='w',
            markerfacecolor=C_ENH,
            markersize=8,
            label=f"Enhancer cohort (n={ENH})"
        ),
        Line2D(
            [0], [0],
            marker='o', color='w',
            markerfacecolor=C_FLANK,
            markersize=8,
            label=f"Flank-only cohort (n={FLANK})"
        ),
    ]
    axes[1].legend(
        handles=right_handles,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(0.90, 0.5)
    )

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
plt.show()





import matplotlib as mpl
mpl.rcParams["font.family"] = "Arial"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- Input file ----
# legacy external path example removed
df = pd.read_csv(csv_file)

# ---- Clean columns ----
df["source"] = df["source"].astype(str).str.strip().str.lower()
df["flank_hit_class"] = df["flank_hit_class"].fillna("").astype(str).str.strip().str.lower()

# ---- Count rows for enhancer-hit pairs ----
# Rules:
# - source == enhancer_hit  -> goes under "Enhancer-hit pairs"
# - flank_hit_class == single -> counted as 1 flank hit
# - flank_hit_class == double -> counted as 2 flank hits
# - otherwise -> counted as 0 flank hits

enh_df = df[df["source"] == "enhancer_hit"].copy()

enh_double = (enh_df["flank_hit_class"] == "double").sum()
enh_single = (enh_df["flank_hit_class"] == "single").sum()
enh_none = (~enh_df["flank_hit_class"].isin(["single", "double"])).sum()
enh_total = enh_double + enh_single + enh_none

# ---- Count rows for flank cohort ----
# Rule:
# - source == flank -> goes under "Flank cohort"
flank_df = df[df["source"] == "flank"].copy()
flank_total = len(flank_df)

# Optional: if you want to confirm what kinds are present in flank rows
flank_single = (flank_df["flank_hit_class"] == "single").sum()
flank_double = (flank_df["flank_hit_class"] == "double").sum()
flank_none = (~flank_df["flank_hit_class"].isin(["single", "double"])).sum()

print("Enhancer-hit pairs")
print(f"  2 flank hits: {enh_double}")
print(f"  1 flank hit : {enh_single}")
print(f"  0 flank hits: {enh_none}")
print(f"  total       : {enh_total}")
print()
print("Flank cohort")
print(f"  total       : {flank_total}")
print(f"    single    : {flank_single}")
print(f"    double    : {flank_double}")
print(f"    none      : {flank_none}")

# ---- Colors ----
C_DOUBLE = "#D35400"   # dark orange
C_SINGLE = "#F8C471"   # light orange
C_NONE   = "#E5E7E9"   # light gray

# ---- Plot ----
fig, ax = plt.subplots(figsize=(2, 2))

labels = ["Enhancer-hit pairs", "Flank cohort"]
x = np.arange(len(labels))

# ========== Enhancer bar (stacked vertical): DOUBLE -> SINGLE -> ZERO ==========
bottom = 0
ax.bar(
    x[0], enh_double, bottom=bottom,
    color=C_DOUBLE, edgecolor="white", linewidth=1.2
)
bottom += enh_double

ax.bar(
    x[0], enh_single, bottom=bottom,
    color=C_SINGLE, edgecolor="white", linewidth=1.2
)
bottom += enh_single

ax.bar(
    x[0], enh_none, bottom=bottom,
    color=C_NONE, edgecolor="white", linewidth=1.2
)

# ========== Flank bar ==========
# Per your rule: anything in source == flank goes in the flank cohort
ax.bar(
    x[1], flank_total, bottom=0,
    color=C_DOUBLE, edgecolor="white", linewidth=1.2
)

# ---- Labels inside segments ----
def label_segment(xpos, start, height, text):
    if height > 0:
        ax.text(
            xpos, start + height / 2, text,
            ha="center", va="center", fontsize=8, color="black"
        )

# Enhancer labels
start = 0
label_segment(x[0], start, enh_double, f"n={enh_double}")
start += enh_double
label_segment(x[0], start, enh_single, f"n={enh_single}")
start += enh_single
label_segment(x[0], start, enh_none, f"n={enh_none}")

# Flank label
label_segment(x[1], 0, flank_total, f"n={flank_total}")

# ---- Remove axes entirely ----
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# ---- Cohort labels below bars ----
max_total = max(enh_total, flank_total)
ax.text(x[0], -max_total * 0.06, labels[0], ha="center", va="top", fontsize=8)
ax.text(x[1], -max_total * 0.06, labels[1], ha="center", va="top", fontsize=8)

# ---- Limits ----
ax.set_xlim(-0.6, len(labels) - 0.4)
ax.set_ylim(-max_total * 0.14, max_total * 1.02)

# ---- Save ----
plt.savefig(OUTPUT_DIR / "stacked_bars_enhancer_vs_flank_pairs_vertical_reordered_2.png", dpi=600)
plt.savefig(OUTPUT_DIR / "stacked_bars_enhancer_vs_flank_pairs_vertical_reordered_2.pdf", dpi=600, bbox_inches="tight")
plt.show()


import matplotlib as mpl
mpl.rcParams["font.family"] = "Arial"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- Input ----

df = pd.read_csv(csv_file)

# ---- Clean ----
df["source"] = df["source"].astype(str).str.strip().str.lower()
df["flank_hit_class"] = (
    df["flank_hit_class"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.lower()
)

# ---- Apply your rules ----
# enhancer_hit rows -> enhancer-hit pairs
enh_df = df[df["source"] == "enhancer_hit"].copy()

enh_double = (enh_df["flank_hit_class"] == "double").sum()
enh_single = (enh_df["flank_hit_class"] == "single").sum()
enh_none = ~enh_df["flank_hit_class"].isin(["single", "double"])
enh_none = enh_none.sum()

enh_total = enh_double + enh_single + enh_none

# flank rows -> flank cohort
flank_df = df[df["source"] == "flank"].copy()
flank_total = len(flank_df)

print("Enhancer-hit pairs")
print("2 flank hits:", enh_double)
print("1 flank hit :", enh_single)
print("0 flank hits:", enh_none)
print("total       :", enh_total)

print("\nFlank cohort")
print("total       :", flank_total)

# ---- Colors ----
C_DOUBLE = "#D35400"   # dark orange
C_SINGLE = "#F8C471"   # light orange
C_NONE   = "#E5E7E9"   # light gray

# ---- Plot ----
fig, ax = plt.subplots(figsize=(2.2, 2.2))

labels = ["Enhancer-hit pairs", "Flank cohort"]
x = np.arange(len(labels))

# Enhancer stacked bar
bottom = 0
ax.bar(x[0], enh_double, bottom=bottom, color=C_DOUBLE, edgecolor="white", linewidth=1.2)
bottom += enh_double

ax.bar(x[0], enh_single, bottom=bottom, color=C_SINGLE, edgecolor="white", linewidth=1.2)
bottom += enh_single

ax.bar(x[0], enh_none, bottom=bottom, color=C_NONE, edgecolor="white", linewidth=1.2)

# Flank cohort bar
ax.bar(x[1], flank_total, bottom=0, color=C_DOUBLE, edgecolor="white", linewidth=1.2)

# ---- Labels inside segments ----
def label_segment(xpos, start, height, text):
    if height > 0:
        ax.text(
            xpos,
            start + height / 2,
            text,
            ha="center",
            va="center",
            fontsize=8,
            color="black"
        )

start = 0
label_segment(x[0], start, enh_double, f"n={enh_double}")
start += enh_double

label_segment(x[0], start, enh_single, f"n={enh_single}")
start += enh_single

label_segment(x[0], start, enh_none, f"n={enh_none}")

#label_segment(x[1], 0, flank_total, f"n={flank_total}")
plot_flank_height = max(flank_total, max_total * 0.03)
ax.bar(x[1], plot_flank_height, bottom=0, color=C_DOUBLE, edgecolor="white", linewidth=1.2)
label_segment(x[1], 0, plot_flank_height, f"n={flank_total}")
# ---- Remove axes ----
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# ---- Labels below bars ----
max_total = max(enh_total, flank_total)
ax.text(x[0], -max_total * 0.06, labels[0], ha="center", va="top", fontsize=8)
ax.text(x[1], -max_total * 0.06, labels[1], ha="center", va="top", fontsize=8)

ax.set_xlim(-0.6, len(labels) - 0.4)
ax.set_ylim(-max_total * 0.14, max_total * 1.02)


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# ====== Font / figure settings ======
mpl.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# ====== Input file ======
# legacy external path example removed

# ====== Read data ======
df = pd.read_csv(csv_file)
df.columns = [c.strip() for c in df.columns]

# Expecting columns like:
# gene_name, pair, source, flank_hit_class
required_cols = {"gene_name", "pair"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Drop missing gene_name/pair rows
df = df.dropna(subset=["gene_name", "pair"]).copy()

# Normalize text
df["gene_name"] = df["gene_name"].astype(str).str.strip()
df["pair"] = df["pair"].astype(str).str.strip()

# Keep only unique gene_name/pair combinations
pair_per_gene = df[["gene_name", "pair"]].drop_duplicates()

# Count number of unique pairs per gene set
pair_counts = (
    pair_per_gene.groupby("gene_name")["pair"]
    .nunique()
    .reset_index(name="n_unique_pairs")
)

# ====== Bin gene sets by number of unique pairs ======
n_1 = (pair_counts["n_unique_pairs"] == 1).sum()
n_23 = pair_counts["n_unique_pairs"].between(2, 3, inclusive="both").sum()
n_ge4 = (pair_counts["n_unique_pairs"] >= 4).sum()

TOTAL = n_1 + n_23 + n_ge4

print(f"Gene sets with 1 pair: {n_1}")
print(f"Gene sets with 2–3 pairs: {n_23}")
print(f"Gene sets with ≥4 pairs: {n_ge4}")
print(f"Total gene sets: {TOTAL}")

# ====== Appearance ======
C1   = "#F1B6DA"  # light magenta
C23  = "#C51B7D"  # medium magenta
CGE4 = "#7A0177"  # dark magenta

DPI = 600

# ====== Build figure ======
fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=DPI)

values = [n_1, n_23, n_ge4]
labels = [
    f"1 pair/set\n{n_1} sets ({100*n_1/TOTAL:.1f}%)",
    f"2–3 pairs/set\n{n_23} sets ({100*n_23/TOTAL:.1f}%)",
    f"≥4 pairs/set\n{n_ge4} sets ({100*n_ge4/TOTAL:.1f}%)",
]
colors = [C1, C23, CGE4]

ax.pie(
    values,
    labels=labels,
    colors=colors,
    startangle=90,
    counterclock=False,
    labeldistance=1.15,
    wedgeprops=dict(width=0.35, edgecolor="black", linewidth=0.3)
)

ax.set(aspect="equal")
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

# ====== Save outputs ======
plt.savefig(
    OUTPUT_DIR / "donut_unique_pairs_per_geneset.png",
    dpi=DPI,
    bbox_inches="tight"
)
plt.savefig(
    OUTPUT_DIR / "donut_unique_pairs_per_geneset.pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# =============================================================================
# PLOT 3: MOUSE BAR PLOT
# Counts unique enhancers in each shadow-set-size bin that are implicated in
# a duplication event. An enhancer can contribute to multiple bins if it
# appears in sets from different bins.
# =============================================================================

from pathlib import Path
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


OUT_PNG = str(OUTPUT_DIR / "enhancers_per_bin.png")
OUT_PDF = str(OUTPUT_DIR / "enhancers_per_bin.pdf")

# -----------------------------
# style
# -----------------------------
DPI = 600
BASE_FS = 12
PURPLE = "#cdb4db"

mpl.rcParams.update({
    "font.size": BASE_FS,
    "axes.titlesize": BASE_FS,
    "axes.labelsize": BASE_FS,
    "xtick.labelsize": BASE_FS,
    "ytick.labelsize": BASE_FS,
    "legend.fontsize": BASE_FS,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# -----------------------------
# helpers
# -----------------------------
def to_bucket(n):
    if n == 2:
        return "2 shadows/set"
    elif n == 3:
        return "3 shadows/set"
    elif n >= 4:
        return "≥4 shadows/set"
    return None

def parse_pair_field(pair_str):
    pair_str = str(pair_str).strip()
    if not pair_str or "_" not in pair_str:
        return None
    a, b = pair_str.split("_", 1)
    return a.strip(), b.strip()

# -----------------------------
# load bed
# enhancer ID = chr:start-end
# set ID = column 4
# -----------------------------
bed = pd.read_csv(SHADOW_BED, sep="\t", header=None)

if bed.shape[1] < 4:
    raise ValueError("BED file must have at least 4 columns.")

bed = bed.iloc[:, :4].copy()
bed.columns = ["chrom", "start", "end", "set_name"]

bed["enhancer_id"] = (
    bed["chrom"].astype(str).str.strip() + ":" +
    bed["start"].astype(str).str.strip() + "-" +
    bed["end"].astype(str).str.strip()
)
bed["set_name"] = bed["set_name"].astype(str).str.strip()

# count shadows per set
set_sizes = bed["set_name"].value_counts().to_dict()
set_to_bucket = {
    set_name: to_bucket(size)
    for set_name, size in set_sizes.items()
    if to_bucket(size) is not None
}

bed["bucket"] = bed["set_name"].map(set_to_bucket)
bed = bed.dropna(subset=["bucket"])

bucket_order = ["2 shadows/set", "3 shadows/set", "≥4 shadows/set"]

# -----------------------------
# load final csv
# pair column contains enhancer1_enhancer2
# both enhancers must be counted
# -----------------------------
final = pd.read_csv(FINAL_BREAKDOWN_CSV)
final.columns = [c.strip() for c in final.columns]

if "pair" not in final.columns:
    raise ValueError("CSV must contain a 'pair' column.")

dup_enhancers = set()
for pair_val in final["pair"].dropna():
    parsed = parse_pair_field(pair_val)
    if parsed is None:
        continue
    enh1, enh2 = parsed
    dup_enhancers.add(enh1)
    dup_enhancers.add(enh2)

# -----------------------------
# collect unique enhancers per bin
# denominator = all unique enhancers in that bin
# numerator   = unique enhancers in that bin found in dup_enhancers
# -----------------------------
bucket_to_all = {b: set() for b in bucket_order}
bucket_to_hit = {b: set() for b in bucket_order}

for _, row in bed.iterrows():
    bucket = row["bucket"]
    enhancer_id = row["enhancer_id"]

    bucket_to_all[bucket].add(enhancer_id)

    if enhancer_id in dup_enhancers:
        bucket_to_hit[bucket].add(enhancer_id)

totals = [len(bucket_to_all[b]) for b in bucket_order]
hits = [len(bucket_to_hit[b]) for b in bucket_order]
proportions = [(h / t) if t > 0 else np.nan for h, t in zip(hits, totals)]

print("Proportion of unique enhancers implicated in duplication, by set-size bin:")
for b, h, t, p in zip(bucket_order, hits, totals, proportions):
    print(f"{b}: {h}/{t} = {p:.4f}")

# -----------------------------
# plot
# -----------------------------
fig, ax = plt.subplots(figsize=(3.5, 3.4), dpi=DPI)

x = np.arange(len(bucket_order))
ax.bar(x, proportions, color=PURPLE, edgecolor="black", width=0.78)

ax.set_xticks(x)
ax.set_xticklabels([f"{b}\n(n={t})" for b, t in zip(bucket_order, totals)],
                   rotation=45, ha="right")
ax.set_ylabel("Proportion of hit enhancers")

ymax = max([p for p in proportions if pd.notna(p)] + [0])
ax.set_ylim(0, ymax + 0.08)

# for xi, yi, h, t in zip(x, proportions, hits, totals):
#     if pd.notna(yi):
#         ax.text(xi, yi + 0.01, f"{h}/{t}", ha="center", va="bottom")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", direction="out", length=4, width=1.0)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
plt.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
plt.show()



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})
# ---- Input file ----
csv_file = BREAKDOWN_FILE
df = pd.read_csv(csv_file)

# ---- Clean columns ----
df["source"] = df["source"].astype(str).str.strip().str.lower()
df["flank_hit_class"] = df["flank_hit_class"].fillna("").astype(str).str.strip().str.lower()

# ---- Count rows for enhancer-hit pairs ----
enh_df = df[df["source"] == "enhancer_hit"].copy()

enh_double = (enh_df["flank_hit_class"] == "double").sum()
enh_single = (enh_df["flank_hit_class"] == "single").sum()
enh_none = (~enh_df["flank_hit_class"].isin(["single", "double"])).sum()
enh_total = enh_double + enh_single + enh_none

# ---- Count rows for flank cohort ----
flank_df = df[df["source"] == "flank"].copy()
flank_total = len(flank_df)

flank_single = (flank_df["flank_hit_class"] == "single").sum()
flank_double = (flank_df["flank_hit_class"] == "double").sum()
flank_none = (~flank_df["flank_hit_class"].isin(["single", "double"])).sum()

print("Enhancer-hit pairs")
print(f"  2 flank hits: {enh_double}")
print(f"  1 flank hit : {enh_single}")
print(f"  0 flank hits: {enh_none}")
print(f"  total       : {enh_total}")
print()
print("Flank cohort")
print(f"  total       : {flank_total}")
print(f"    single    : {flank_single}")
print(f"    double    : {flank_double}")
print(f"    none      : {flank_none}")

# ---- Colors ----
C_DOUBLE = "#D35400"   # dark orange
C_SINGLE = "#F8C471"   # light orange
C_NONE   = "#E5E7E9"   # light gray

# ---- Plot ----
fig, ax = plt.subplots(figsize=(2, 2))

labels = ["Enhancer-hit pairs", "Flank cohort"]
x = np.arange(len(labels))

# ========== Enhancer bar (stacked vertical): DOUBLE -> SINGLE -> ZERO ==========
bottom = 0
ax.bar(
    x[0], enh_double, bottom=bottom,
    color=C_DOUBLE, edgecolor="white", linewidth=1.2
)
bottom += enh_double

ax.bar(
    x[0], enh_single, bottom=bottom,
    color=C_SINGLE, edgecolor="white", linewidth=1.2
)
bottom += enh_single

ax.bar(
    x[0], enh_none, bottom=bottom,
    color=C_NONE, edgecolor="white", linewidth=1.2
)

# ========== Flank bar ==========
ax.bar(
    x[1], flank_total, bottom=0,
    color=C_DOUBLE, edgecolor="white", linewidth=1.2
)

# ---- Labels inside segments: print percent instead of n ----
def label_segment_percent(xpos, start, height, total, color="black"):
    if height > 0 and total > 0:
        pct = 100 * height / total
        ax.text(
            xpos, start + height / 2, f"{pct:.1f}%",
            ha="center", va="center", fontsize=7, color=color
        )

# Enhancer labels
start = 0
label_segment_percent(x[0], start, enh_double, enh_total)
start += enh_double
label_segment_percent(x[0], start, enh_single, enh_total, color="black")
start += enh_single
label_segment_percent(x[0], start, enh_none, enh_total, color="black")

# Flank label
label_segment_percent(x[1], 0, flank_total, flank_total)

# ---- Remove axes entirely ----
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# ---- Cohort labels below bars ----
max_total = max(enh_total, flank_total)
ax.text(x[0], -max_total * 0.06, labels[0], ha="center", va="top", fontsize=8)
ax.text(x[1], -max_total * 0.06, labels[1], ha="center", va="top", fontsize=8)

# ---- Limits ----
ax.set_xlim(-0.6, len(labels) - 0.4)
ax.set_ylim(-max_total * 0.14, max_total * 1.02)

# ---- Save ----
plt.savefig(OUTPUT_DIR / "stacked_bars_enhancer_vs_flank_pairs_vertical_percent.png", dpi=600, bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "stacked_bars_enhancer_vs_flank_pairs_vertical_percent.pdf", dpi=600, bbox_inches="tight")
plt.show()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- Input file ----
csv_file = BREAKDOWN_FILE
df = pd.read_csv(csv_file)

# ---- Clean columns ----
df["source"] = df["source"].astype(str).str.strip().str.lower()
df["flank_hit_class"] = df["flank_hit_class"].fillna("").astype(str).str.strip().str.lower()
df["gene_name"] = df["gene_name"].astype(str).str.strip()

# =========================================================
# BIN GENE NAMES FOR ENHANCER_HIT ROWS
# Each gene_name goes into exactly one bin:
#   double > single > none
# =========================================================
enh_df = df[df["source"] == "enhancer_hit"].copy()

gene_bins = []
for gene, g in enh_df.groupby("gene_name"):
    has_double = (g["flank_hit_class"] == "double").any()
    has_single = (g["flank_hit_class"] == "single").any()

    if has_double:
        gene_bin = "double"
    elif has_single:
        gene_bin = "single"
    else:
        gene_bin = "none"

    gene_bins.append((gene, gene_bin))

gene_bin_df = pd.DataFrame(gene_bins, columns=["gene_name", "gene_bin"])

enh_double = (gene_bin_df["gene_bin"] == "double").sum()
enh_single = (gene_bin_df["gene_bin"] == "single").sum()
enh_none = (gene_bin_df["gene_bin"] == "none").sum()
enh_total = len(gene_bin_df)

# ---- Flank cohort: unique gene_names in source == flank ----
flank_df = df[df["source"] == "flank"].copy()
flank_total = flank_df["gene_name"].nunique()

print("Enhancer-hit gene bins")
print(f"  2 flank hits bin: {enh_double}")
print(f"  1 flank hit bin : {enh_single}")
print(f"  0 flank hits bin: {enh_none}")
print(f"  total genes     : {enh_total}")
print()
print("Flank cohort")
print(f"  unique genes    : {flank_total}")

# ---- Colors ----
C_DOUBLE = "#D35400"   # dark orange
C_SINGLE = "#F8C471"   # light orange
C_NONE   = "#E5E7E9"   # light gray

# ---- Plot ----
fig, ax = plt.subplots(figsize=(2, 2))

labels = ["Enhancer-hit genes", "Flank cohort"]
x = np.arange(len(labels))

# ========== Enhancer bar (stacked vertical): DOUBLE -> SINGLE -> NONE ==========
bottom = 0
ax.bar(
    x[0], enh_double, bottom=bottom,
    color=C_DOUBLE, edgecolor="white", linewidth=1.2
)
bottom += enh_double

ax.bar(
    x[0], enh_single, bottom=bottom,
    color=C_SINGLE, edgecolor="white", linewidth=1.2
)
bottom += enh_single

ax.bar(
    x[0], enh_none, bottom=bottom,
    color=C_NONE, edgecolor="white", linewidth=1.2
)

# ========== Flank bar ==========
ax.bar(
    x[1], flank_total, bottom=0,
    color=C_DOUBLE, edgecolor="white", linewidth=1.2
)

# ---- Labels inside segments: print % of gene_names ----
def label_segment_percent(xpos, start, height, total, color="black"):
    if height > 0 and total > 0:
        pct = 100 * height / total
        ax.text(
            xpos, start + height / 2, f"{pct:.1f}%",
            ha="center", va="center", fontsize=7, color=color
        )

# Enhancer labels
start = 0
label_segment_percent(x[0], start, enh_double, enh_total)
start += enh_double
label_segment_percent(x[0], start, enh_single, enh_total, color="black")
start += enh_single
label_segment_percent(x[0], start, enh_none, enh_total, color="black")

# Flank label
label_segment_percent(x[1], 0, flank_total, flank_total)

# ---- Remove axes entirely ----
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# ---- Cohort labels below bars ----
max_total = max(enh_total, flank_total)
ax.text(x[0], -max_total * 0.06, labels[0], ha="center", va="top", fontsize=8)
ax.text(x[1], -max_total * 0.06, labels[1], ha="center", va="top", fontsize=8)

# ---- Limits ----
ax.set_xlim(-0.6, len(labels) - 0.4)
ax.set_ylim(-max_total * 0.14, max_total * 1.02)

# ---- Save ----
plt.savefig(OUTPUT_DIR / "fly_stacked_bars_gene_binned_percent.png", dpi=600, bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "fly_stacked_bars_gene_binned_percent.pdf", dpi=600, bbox_inches="tight")
plt.show()
