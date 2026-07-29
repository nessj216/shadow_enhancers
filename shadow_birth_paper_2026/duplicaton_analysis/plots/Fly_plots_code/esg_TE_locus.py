import pandas as pd
import matplotlib.pyplot as plt
import re
from pathlib import Path

# ----------------------------
# input
# ----------------------------
PAPER_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = Path(__file__).resolve().parent / "input"
OUTPUT_DIR = PAPER_ROOT / "plots" / "output_pngs" / "fly"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

csv_path = str(INPUT_DIR / "esg_TE_raw_splice_sizes_across_lines.csv")

# ----------------------------
# load data
# ----------------------------
df = pd.read_csv(csv_path)

# shorten genome names
def shorten_name(s):
    if s == "Dmel_ISO1":
        return s
    m = re.search(r"ASM(\d+)v1", str(s))
    if m:
        return "ASM" + m.group(1)
    m2 = re.search(r"GCA_(\d+\.\d+)", str(s))
    if m2:
        return "GCA_" + m2.group(1)
    return str(s)

df["short_sample"] = df["sample"].astype(str).apply(shorten_name)
# more intuitive labels for x-axis
sample_label_map = {
    "ASM2014149": "STO-022",
    "ASM2014150": "TEN-015",
    "ASM2014151": "TOM-008",
    "ASM2014157": "RAL-855",
    "ASM2014159": "SLA-001",
    "ASM2014162": "RAL-426",
    "ASM2014170": "RAL-091",
    "ASM2014174": "RAL-059",
    "ASM2014198": "GIM-012",
    "ASM2014204": "COR-018",

    # fallback if short_sample is GCA instead of ASM
    "GCA_020141495.1": "STO-022",
    "GCA_020141505.1": "TEN-015",
    "GCA_020141515.1": "TOM-008",
    "GCA_020141575.1": "RAL-855",
    "GCA_020141595.1": "SLA-001",
    "GCA_020141625.1": "RAL-426",
    "GCA_020141705.1": "RAL-091",
    "GCA_020141745.1": "RAL-059",
    "GCA_020141985.1": "GIM-012",
    "GCA_020142045.1": "COR-018",
}

df["plot_label"] = df["short_sample"].replace(sample_label_map)
# ----------------------------
# plot style
# ----------------------------
plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"]
})

ref_rows = df[df["sample"].isin(["DM6", "Dmel_ISO1"])]

if len(ref_rows) == 0:
    raise ValueError("Could not find DM6 or Dmel_ISO1 reference row.")

dmel_size = float(ref_rows["raw_splice_gap_bp"].iloc[0])

# remove reference from plotted samples
plot_df = df[~df["sample"].isin(["DM6", "Dmel_ISO1"])].copy()
# ----------------------------
# make plot
# ----------------------------
# ----------------------------
# make plot: lollipop/dot plot
# ----------------------------

plot_df["frac_recovered"] = plot_df["raw_splice_gap_bp"] / dmel_size

plot_df["recovery_status"] = plot_df["frac_recovered"].apply(
    lambda x: "full-length" if x >= 0.9 else "partial recovery"
)

# optional: sort by recovered size
plot_df = plot_df.sort_values("raw_splice_gap_bp", ascending=False).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(5, 2), dpi=300)

x = range(len(plot_df))

# lollipop stems
# ax.vlines(
#     x,
#     ymin=0,
#     ymax=plot_df["raw_splice_gap_bp"],
#     linewidth=0.8,
#     alpha=0.6
# )

colors = {
    "full-length": "#1f77b4",
    "partial recovery": "#999999"
}

for status, sub in plot_df.groupby("recovery_status"):
    ax.scatter(
        sub.index,
        sub["raw_splice_gap_bp"],
        s=45,
        color=colors[status],
        label=status,
        zorder=3
    )

# ISO1/DM6 reference length shown only as dashed line
ax.axhline(
    dmel_size,
    linestyle="--",
    linewidth=1,
    color="black",
    label="ISO1/DM6 TE size"
)

# optional 90–110% recovery band
ax.axhspan(
    dmel_size * 0.9,
    dmel_size * 1.1,
    color="gray",
    alpha=0.15,
    linewidth=0
)

ax.set_ylabel(r"Recovered $\mathit{esg}$ locus"+"\n "+"TE size (bp)")
ax.set_xlabel("")
#ax.set_title("esg TE recovery across lines")

ax.set_xticks(list(x))
ax.set_xticklabels(plot_df["plot_label"], rotation=45, ha="right")

ax.set_ylim(0, dmel_size * 1.18)

ax.legend(
    frameon=False,
    fontsize=9,
    loc="upper left",
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)
plt.tight_layout()
#plt.tight_layout(rect=[0, 0, 0.78, 1])
plt.savefig(
    OUTPUT_DIR / "esglocus.png",
    dpi=600,
    bbox_inches="tight"
)
plt.show()


