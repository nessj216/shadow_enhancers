from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"
PLOT_FIGURES = ROOT / "plot_figures"
PLOT_FIGURES.mkdir(parents=True, exist_ok=True)

SHADOW_BED = str(INPUT / "011925_all_shadowsets_DM6.bed")
TE_BED = str(INPUT / "Final_overlap_shadows_TEtype.bed")
OUTDIR = str(PLOT_FIGURES)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial"
})
def to_bucket(n: int) -> str | None:
    if n == 2:
        return "2 shadows/set"
    if n == 3:
        return "3 shadows/set"
    if n >= 4:
        return "≥4 shadows/set"
    return None


def bonferroni(pvals):
    m = len(pvals)
    return [min(p * m, 1.0) for p in pvals]


def plot_te_enh_hit_by_shadow_category(
    shadow_bed: str = SHADOW_BED,
    te_bed: str = TE_BED,
    outbase: str = str(PLOT_FIGURES / "TE_enhancer_hits_by_shadow_category_2"),
) -> None:
    # 1) Read files
    shadow = pd.read_csv(
        shadow_bed,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "set_name"],
    )

    # robust parser for ragged TE file:
    # enhancer info = first 4 columns
    # overlap = last column
    te_rows = []
    with open(te_bed, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                print(f"Skipping line {line_num}: fewer than 5 columns")
                continue

            te_rows.append({
                "chr": parts[0],
                "start": pd.to_numeric(parts[1], errors="coerce"),
                "end": pd.to_numeric(parts[2], errors="coerce"),
                "set_name": str(parts[3]),
                "overlap": pd.to_numeric(parts[-1], errors="coerce"),  # LAST column
            })

    te = pd.DataFrame(te_rows)

    te = te.dropna(subset=["chr", "start", "end", "set_name"])
    te = te[te["overlap"].fillna(0) > 0].copy()    # 2) Build enhancer IDs
    shadow["enhancer_id"] = (
        shadow["chr"].astype(str)
        + ":"
        + shadow["start"].astype(str)
        + "-"
        + shadow["end"].astype(str)
    )

    te["enhancer_id"] = (
        te["chr"].astype(str)
        + ":"
        + te["start"].astype(str)
        + "-"
        + te["end"].astype(str)
    )

    # 3) Count total shadows per set, assign each set to a bucket
    set_sizes = shadow.groupby("set_name").size().rename("shadow_count").reset_index()
    set_sizes["bucket"] = set_sizes["shadow_count"].map(to_bucket)

    shadow = shadow.merge(set_sizes[["set_name", "shadow_count", "bucket"]], on="set_name", how="left")
    te = te.merge(set_sizes[["set_name", "shadow_count", "bucket"]], on="set_name", how="left")

    bucket_order = ["2 shadows/set", "3 shadows/set", "≥4 shadows/set"]

    # 4) TOTAL enhancers per bucket:
    # count BED rows directly (same enhancer can count in different bins if it appears there)
    totals = (
        shadow.groupby("bucket", observed=True)
        .size()
        .reindex(bucket_order, fill_value=0)
        .to_dict()
    )

    # 5) HIT enhancers per bucket:
    # deduplicate duplicate TE+ rows only within the same shadow bin
    # so same enhancer in same bin counted once, but same enhancer in different bins can count in each bin
    te_hits_unique = (
        te.loc[te["bucket"].notna(), ["enhancer_id", "bucket"]]
        .drop_duplicates()
    )

    hits = (
        te_hits_unique.groupby("bucket", observed=True)
        .size()
        .reindex(bucket_order, fill_value=0)
        .to_dict()
    )

    proportions = [hits[b] / totals[b] if totals[b] else np.nan for b in bucket_order]

    # 6) Fisher exact tests + Bonferroni
    hits_row = [hits[b] for b in bucket_order]
    nonhits_row = [max(totals[b] - hits[b], 0) for b in bucket_order]

    comps = [(0, 1), (0, 2), (1, 2)]
    p_raw = []
    for i, j in comps:
        table = [
            [hits_row[i], nonhits_row[i]],
            [hits_row[j], nonhits_row[j]],
        ]
        _, p = fisher_exact(table, alternative="two-sided")
        p_raw.append(p)
    p_adj = bonferroni(p_raw)

    # 7) Plot
    colors = list(plt.cm.Purples([0.6, 0.6, 0.6]))

    fig, ax = plt.subplots(figsize=(4, 4))
    x = np.arange(len(bucket_order))

    ax.bar(
        x,
        proportions,
        color=colors,
        alpha=0.90,
        edgecolor="black",
        linewidth=1.2,
    )

    ax.set_xticks(x)

    ax.set_xticklabels([
        f"2 shadows/set\n(n={totals['2 shadows/set']})",
        f"3 shadows/set\n(n={totals['3 shadows/set']})",
        f"≥4 shadows/set\n(n={totals['≥4 shadows/set']})",
    ], rotation=45, ha="right")
    ax.set_ylabel("Proportion of TE+ enhancers")
    #ax.set_xlabel("Shadow enhancers per set", labelpad=26)

    ax.yaxis.grid(False)
    ax.xaxis.grid(False)

    # for i, b in enumerate(bucket_order):
    #     ax.text(
    #         i,
    #         -0.12,
    #         f"(n={totals[b]})",
    #         ha="center",
    #         va="top",
    #         transform=ax.get_xaxis_transform()
    #     )
    #     ax.text(
    #         i,
    #         proportions[i] + 0.004,
    #         f"{hits[b]}/{totals[b]} = {proportions[i]*100:.2f}%",
    #         ha="center",
    #         va="bottom"
    #     )

    y_max = np.nanmax(proportions)
    ax.set_ylim(0, y_max + 0.02)
    # keep only left and bottom axes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    # make axes black
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)

    # y-axis tick marks on, black
    ax.tick_params(axis="y", which="both", left=True, right=False,
                   length=4, width=1, color="black", labelcolor="black")

    plt.tight_layout()
    plt.savefig(f"{outbase}.png", dpi=300, bbox_inches="tight")
    plt.show()
    #plt.close(fig)

    # 8) Helpful prints
    print("\nTE+ enhancer proportions by shadow-set bin")
    for b in bucket_order:
        prop = hits[b] / totals[b] if totals[b] else float("nan")
        print(f"{b}: hits={hits[b]} totals={totals[b]} prop={prop:.6f}")

    print("\nPairwise Fisher exact tests (raw p, Bonferroni-adjusted p)")
    for (i, j), raw_p, adj_p in zip(comps, p_raw, p_adj):
        print(f"{bucket_order[i]} vs {bucket_order[j]}: raw={raw_p:.6g}, bonf={adj_p:.6g}")


plot_te_enh_hit_by_shadow_category()
