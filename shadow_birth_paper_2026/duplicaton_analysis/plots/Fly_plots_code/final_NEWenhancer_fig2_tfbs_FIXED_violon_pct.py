"""
final_fig2_tfbs_FIXED_violin.py

Fixes:
- ONE save_figure() only (no later redefinitions)
- No bbox_inches="tight", no constrained_layout, no tight_layout
- Layout presets so axes/plot area heights are consistent across figures
- Two-panel k-mer vs TFBS plot uses violin plots (with jitter points) instead of boxplots
"""

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple
import collections
import csv

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.stats import mannwhitneyu, fisher_exact

PAPER_ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUT_DIR = HERE / "input"
OUTPUT_DIR = PAPER_ROOT / "plots" / "output_pngs" / "fly" / "figures_export"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# INPUTS (edit paths)
# =============================================================================
KMER_CSV      = str(PAPER_ROOT / "duplicaton_analysis" / "TFBS_similarity_analysis" / "final_Alignment_with_kmer.no_blacklist__filtered.csv")
# legacy external path example removed
TFBS_SIM_FILE = str(PAPER_ROOT / "duplicaton_analysis" / "TFBS_similarity_analysis" / "tfbs_similarity_three_groups_name_plus_jaccardcounts.tsv")
# Barplot inputs
SHADOW_BED = PAPER_ROOT / "duplicaton_analysis" / "enhancer_hits_per_shadow_bin" / "Fly" / "data" / "011925_all_shadowsets_DM6.bed"
HIT_PAIRS_CSV = INPUT_DIR / "Combined_pairs__gene_first__single-column_pair_.csv"

# Residual scatter "allowed pairs" filter
ALLOWED_PAIRS_FILE = str(INPUT_DIR / "FINAL_breakdown_single_double_flanks_UPDATED_1e4_flankmethod.csv")
ALLOWED_PAIRS_COL = "pair"
FINAL_PAIRS_CSV=ALLOWED_PAIRS_FILE
FIG2_COUNT_MODE = "unique_within_bin"
# legacy external path example removed
#csv_path = DUP_EXTRA_CSV
# =============================================================================
# OUTPUT + STYLE
# =============================================================================
OUTDIR = OUTPUT_DIR
OUTDIR.mkdir(parents=True, exist_ok=True)

DPI = 600
plt.rcParams.update({
    "savefig.dpi": DPI,
    "figure.dpi": DPI,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})

sns.set_style("white")
sns.set_context("paper", rc={
    "axes.titlesize": 16,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
})

# Keep these as you had them
FIGSIZE_TWO_PANEL = (7.0, 3.5)
FIGSIZE_SINGLE = (7.0, 3.5)
np.random.seed(0)

def add_n_under_xticks(ax, df, group_col, order, y=-0.16, fontsize=11):
    """
    Writes '(n=...)' under each x tick. y is in axis coords (0 = axis line).
    More negative -> lower.
    """
    for i, g in enumerate(order):
        n = int(df.loc[df[group_col] == g].shape[0])
        ax.text(
            i, y, f"(n={n})",
            ha="center", va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=fontsize
        )

# =============================================================================
# LAYOUT PRESETS (controls axes/plot-area height)
# - Use NO_TITLE when you do not set ax.set_title(...)
# - Use WITH_TITLE when you DO set titles
# - Use COLORBAR when you have a colorbar (reserves right margin)
# =============================================================================
LAYOUT = {
    "WITH_TITLE": dict(left=0.10, right=0.98, bottom=0.22, top=0.88, wspace=0.35),
    "NO_TITLE":   dict(left=0.10, right=0.95, bottom=0.22, top=0.97, wspace=0.35),
    "COLORBAR":   dict(left=0.10, right=0.88, bottom=0.22, top=0.94, wspace=0.35),
}

def apply_layout(fig, key: str, two_panel: bool = False) -> None:
    kw = LAYOUT[key].copy()
    wspace = kw.pop("wspace", None)
    if two_panel and wspace is not None:
        fig.subplots_adjust(**kw, wspace=wspace)
    else:
        fig.subplots_adjust(**kw)

def save_figure(fig, basename: str, layout: str = "NO_TITLE", two_panel: bool = False) -> None:
    """
    Save with fixed margins (no bbox_inches='tight') so plotted axes heights are consistent.
    """
    apply_layout(fig, layout, two_panel=two_panel)

    png = OUTDIR / f"{basename}.png"
    pdf = OUTDIR / f"{basename}.pdf"
    fig.savefig(png, dpi=DPI)
    fig.savefig(pdf, dpi=DPI)
    plt.close(fig)
    print(f"Saved: {png}\nSaved: {pdf}")

def enforce_ticks(ax):
    ax.tick_params(
        axis="both", which="both",
        direction="out", length=4, width=1.0,
        bottom=True, top=False, left=True, right=False
    )

# =============================================================================
# STATS + ANNOTATION HELPERS
# =============================================================================
def p_to_stars(p):
    if pd.isna(p): return "ns"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"

def add_sig_bar_two_group(ax, p, x1=0, x2=1):
    stars = p_to_stars(p)
    if stars == "ns":
        return
    ymin, ymax = ax.get_ylim()
    y = ymax - (ymax - ymin) * 0.10
    h = (ymax - ymin) * 0.02
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], color="black", linewidth=1)
    ax.text((x1+x2)/2, y + h*1.02, stars, ha="center", va="bottom", fontsize=11)

def bonferroni(pvals: List[float]) -> List[float]:
    m = len(pvals)
    return [min(p * m, 1.0) if not pd.isna(p) else np.nan for p in pvals]

# =============================================================================
# CORE DRAW: 2-group violin + jitter points (replaces boxplots)
# =============================================================================
def draw_two_group_violin_jitter(
    ax,
    df,
    group_col,
    value_col,
    order=("duplication_hits", "neither"),
    labels=("Duplication hits", "no-hit shadows"),
    palette=None,
    ylabel="",
    violin_width=0.80,
    jitter_sd=0.06,
    point_size=20,

):

    if palette is None:
        palette = {"duplication_hits": "#0B3C5D", "neither": "whitesmoke"}

    # MWU
    x = df.loc[df[group_col] == order[0], value_col].dropna().values
    y = df.loc[df[group_col] == order[1], value_col].dropna().values
    p = np.nan
    if len(x) and len(y):
        _, p = mannwhitneyu(x, y, alternative="two-sided")

    # violin
    sns.violinplot(
        data=df, x=group_col, y=value_col,
        order=list(order), palette=palette,
        cut=0, inner="box",
        width=violin_width,
        linewidth=1.0,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_xticklabels(list(labels))
    ax.set_ylabel(ylabel)
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    enforce_ticks(ax)

    add_sig_bar_two_group(ax, p, x1=0, x2=1)
    return p, len(x), len(y)

# =============================================================================
# FIGURE 1: Two-panel k-mer (left) + TFBS (right) with VIOLINS
# =============================================================================
def plot_kmer_tfbs_two_panel_violin(
    kmer_csv=KMER_CSV,
    tfbs_tsv=TFBS_SIM_FILE,
    outbase="fig_kmer_tfbs_dupHits_vs_neither_two_panel_VIOLIN"
):
    kmer = pd.read_csv(kmer_csv)
    tfbs = pd.read_csv(tfbs_tsv, sep="\t")

    # collapse groups
    kmer = kmer[kmer["cohort"].isin(["enhancer_hit", "flank_hit", "neither"])].copy()
    kmer["group2"] = np.where(kmer["cohort"].isin(["enhancer_hit", "flank_hit"]),
                              "duplication_hits", "neither")

    tfbs = tfbs[tfbs["group"].isin(["enhancer_hit", "flank_hit", "neither"])].copy()
    tfbs["group2"] = np.where(tfbs["group"].isin(["enhancer_hit", "flank_hit"]),
                              "duplication_hits", "neither")

    # convert to percent so both panels share the same 0–100% y-scale
    kmer["cosine_pct"] = pd.to_numeric(kmer["cosine"], errors="coerce") * 100.0
    tfbs["jaccard_pct"] = pd.to_numeric(tfbs["mesoderm_jaccard_count"], errors="coerce") * 100.0


    order   = ("duplication_hits", "neither")
    labels  = ("Duplication hits", "no-hit shadows")
    palette = {"duplication_hits":  "#1F4E99", "neither": "whitesmoke"}

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE_TWO_PANEL)

    p_kmer, n_k1, n_k2 = draw_two_group_violin_jitter(
        axes[0], kmer,
        group_col="group2", value_col="cosine_pct",
        order=order, labels=labels, palette=palette,
        ylabel="k-mer cosine similarity (%)",
        violin_width=0.80,
    )
    axes[0].set_ylim(0, 110)
    axes[0].set_yticks([0, 25, 50, 75, 100])


    p_tfbs, n_t1, n_t2 = draw_two_group_violin_jitter(
        axes[1], tfbs,
        group_col="group2", value_col="jaccard_pct",
        order=order, labels=labels, palette=palette,
        ylabel="TFBS Jaccard similarity (%)",
        violin_width=0.80,
    )

    axes[1].set_ylim(0, 110)
    axes[1].set_yticks([0, 25, 50, 75, 100])

    # No titles -> use NO_TITLE layout so axes can be taller
    save_figure(fig, outbase, layout="NO_TITLE", two_panel=True)

    print(f"k-mer MWU p={p_kmer:.3g} (n={n_k1} vs {n_k2})")
    print(f"TFBS  MWU p={p_tfbs:.3g} (n={n_t1} vs {n_t2})")
# =============================================================================
# FIGURE 2: Enhancer hits by shadow category (bar plot)
# =============================================================================
def iter_bed_rows(path: Path) -> Iterable[Tuple[str, str, str, str, List[str]]]:
    with path.open() as sf:
        for line in sf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            chrom, start, end, set_name = cols[0], cols[1], cols[2], cols[3]
            yield chrom, start, end, set_name, cols

def to_bucket(n: int) -> Optional[str]:
    if n == 2:
        return "2 shadows/set"
    if n == 3:
        return "3 shadows/set"
    if n >= 4:
        return "≥4 shadows/set"
    return None

def parse_ids(s: str) -> List[str]:
    s = str(s).strip()
    if s.endswith(".txt"):
        s = s[:-4]
    return [p for p in s.split("_") if p]

def detect_header(first_row: List[str]) -> bool:
    if not first_row:
        return True
    col0 = first_row[0].strip().lower()
    if col0 in {"set", "set_name", "shadow_set"}:
        return True
    if len(first_row) < 2:
        return True
    ids = parse_ids(first_row[1])
    return len(ids) < 2
import csv, collections, re
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Set, Tuple, Iterable, List, Optional
from scipy.stats import fisher_exact

# Uses your existing:
# - SHADOW_BED (Path)
# - save_figure, enforce_ticks
# - FIGSIZE_TWO_PANEL
# - bonferroni
# - to_bucket(), iter_bed_rows()   (from your snippet)
#
# IMPORTANT: set this to your FINAL file path:


def parse_final_pair(s: str) -> Optional[Tuple[str, str]]:
    """
    FINAL file 'pair' strings look like:
      chr3R:16814262-16814590_chr3R:16904976-16905213
    """
    if s is None:
        return None
    s = str(s).strip()
    if not s:
        return None

    # split on "_" most of the time
    if "_" in s:
        parts = [p.strip() for p in s.split("_") if p.strip()]
    elif "/" in s:
        parts = [p.strip() for p in s.split("/") if p.strip()]
    else:
        parts = re.split(r"\s+", s)

    # fallback: explicitly pull out chr:start-end tokens
    if len(parts) != 2:
        toks = re.findall(r"(chr[^_/\s]+:\d+-\d+)", s)
        if len(toks) >= 2:
            parts = toks[:2]

    if len(parts) != 2:
        return None
    return parts[0], parts[1]

import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict
from scipy.stats import fisher_exact

def plot_enh_hit_by_shadow_category(
    outbase: str = "fig_size_bin_purple",
    count_mode: str = FIG2_COUNT_MODE,
) -> None:
    """
    count_mode:
      - "unique_within_bin": count each enhancer once within a shadow-size bin,
        even if the exact interval appears in multiple sets in that bin.
      - "membership": count every BED row / enhancer membership separately.
    """
    bucket_order = ["2 shadows/set", "3 shadows/set", "≥4 shadows/set"]

    shadow_counts: Dict[str, int] = collections.Counter()
    for chrom, start, end, set_name, cols in iter_bed_rows(SHADOW_BED):
        shadow_counts[set_name] += 1

    shadow_buckets: Dict[str, str] = {}
    for set_name, c in shadow_counts.items():
        b = to_bucket(c)
        if b:
            shadow_buckets[set_name] = b

    final = pd.read_csv(FINAL_PAIRS_CSV)
    dup_enhancers: Set[str] = set()
    bad_pairs = 0

    for val in final.get("pair", pd.Series([], dtype=str)).dropna().astype(str):
        parsed = parse_final_pair(val)
        if not parsed:
            bad_pairs += 1
            continue
        a, b = parsed
        dup_enhancers.add(a)
        dup_enhancers.add(b)

    totals_membership: Dict[str, int] = {b: 0 for b in bucket_order}
    hits_membership: Dict[str, int] = {b: 0 for b in bucket_order}
    bucket_to_all: Dict[str, Set[str]] = {b: set() for b in bucket_order}
    bucket_to_hit: Dict[str, Set[str]] = {b: set() for b in bucket_order}

    for chrom, start, end, set_name, cols in iter_bed_rows(SHADOW_BED):
        bucket = shadow_buckets.get(set_name)
        if not bucket:
            continue

        enhancer_id = f"{chrom}:{start}-{end}"
        totals_membership[bucket] += 1
        bucket_to_all[bucket].add(enhancer_id)

        if enhancer_id in dup_enhancers:
            hits_membership[bucket] += 1
            bucket_to_hit[bucket].add(enhancer_id)

    if count_mode == "unique_within_bin":
        totals = {b: len(bucket_to_all[b]) for b in bucket_order}
        hits = {b: len(bucket_to_hit[b]) for b in bucket_order}
    elif count_mode == "membership":
        totals = totals_membership
        hits = hits_membership
    else:
        raise ValueError("count_mode must be 'unique_within_bin' or 'membership'")

    proportions = [hits[b] / totals[b] if totals[b] else np.nan for b in bucket_order]

    comps = [(0, 1), (0, 2), (1, 2)]
    p_raw = []
    for i, j in comps:
        bi = bucket_order[i]
        bj = bucket_order[j]
        table = [
            [hits[bi], max(totals[bi] - hits[bi], 0)],
            [hits[bj], max(totals[bj] - hits[bj], 0)],
        ]
        _, p = fisher_exact(table, alternative="two-sided")
        p_raw.append(p)
    p_adj = bonferroni(p_raw)

    fig, axes = plt.subplots(
        1, 2,
        figsize=FIGSIZE_TWO_PANEL,
        gridspec_kw={"width_ratios": [1.65, 0.35]}
    )
    ax = axes[0]
    axes[1].axis("off")

    x = np.arange(len(bucket_order))
    colors = list(plt.cm.Purples([0.55, 0.65, 0.78]))
    ax.bar(x, proportions, color=colors, alpha=0.95, edgecolor="white", linewidth=1.2)

    ax.set_xticks(x)
    ax.set_xticklabels(bucket_order)
    ax.set_ylabel("Proportion of hit enhancers")
    ax.set_xlabel("Shadow enhancers per set", labelpad=26)
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    enforce_ticks(ax)

    for i, b in enumerate(bucket_order):
        ax.text(
            i, -0.12,
            f"(n={totals[b]})",
            ha="center", va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=11
        )

    y_max = np.nanmax(proportions)
    ax.set_ylim(0, y_max + 0.08)

    save_figure(fig, outbase, layout="NO_TITLE", two_panel=True)

    repeated_within_bin = {
        b: totals_membership[b] - len(bucket_to_all[b]) for b in bucket_order
    }
    print(f"[Figure 2] count_mode={count_mode}")
    print(f"[Figure 2] bad_rows_in_final_pair_file={bad_pairs}")
    for b in bucket_order:
        print(
            f"{b}: hits={hits[b]} totals={totals[b]} prop={hits[b]/totals[b]:.3f} | "
            f"membership_total={totals_membership[b]} unique_total={len(bucket_to_all[b])} "
            f"repeated_within_bin={repeated_within_bin[b]}"
        )
    print("Raw p:", p_raw)
    print("Bonferroni p:", p_adj)

# =============================================================================
# FIGURE 3: TFBS residual-z scatter with gray dead-zone (fixed height, reserves colorbar)
# =============================================================================
def _round_sig(x, sig=2):
    if x == 0 or not np.isfinite(x):
        return x
    return float(f"{x:.{sig}g}")

def plot_tfbs_residual_z_scatter_deadzone(
    tfbs_sim_file=TFBS_SIM_FILE,
    align_file=KMER_CSV,
    outbase="fig_tfbs_resid_z_scatter_deadzone",
    deg=2,
    deadzone=1.0,
    allowed_pairs_file=ALLOWED_PAIRS_FILE,
    allowed_pairs_col=ALLOWED_PAIRS_COL,
):
    S_POINTS = 30

    def norm_label(x: str) -> str:
        s = str(x).strip()
        return s.replace(":", "_").replace("|", "_").replace(" ", "")

    def pair_key(a: str, b: str) -> str:
        a2, b2 = norm_label(a), norm_label(b)
        return "__".join(sorted([a2, b2]))

    def parse_pair_field(s: str):
        s = str(s).strip()
        if not s:
            return None
        parts = s.split("_")
        if len(parts) != 2:
            return None
        return parts[0], parts[1]

    pairs_df = pd.read_csv(allowed_pairs_file)
    if allowed_pairs_col not in pairs_df.columns:
        raise ValueError(
            f"Expected column '{allowed_pairs_col}' in {allowed_pairs_file}, "
            f"found: {list(pairs_df.columns)}"
        )

    allowed_pair_keys = set()
    bad = 0
    for s in pairs_df[allowed_pairs_col].dropna():
        ab = parse_pair_field(s)
        if ab is None:
            bad += 1
            continue
        a, b = ab
        allowed_pair_keys.add(pair_key(a, b))

    if len(allowed_pair_keys) == 0:
        raise ValueError(
            f"No valid pairs parsed from {allowed_pairs_file} column '{allowed_pairs_col}'. "
            f"(bad_rows={bad})"
        )

    sim_all = pd.read_csv(tfbs_sim_file, sep="\t")
    align = pd.read_csv(align_file)

    sim_all["pair_key_u"] = [pair_key(a, b) for a, b in zip(sim_all["enhancer_a"], sim_all["enhancer_b"])]
    align["pair_key_u"]   = [pair_key(a, b) for a, b in zip(align["nameA"], align["nameB"])]

    sim_all_f = sim_all[sim_all["pair_key_u"].isin(allowed_pair_keys)].copy()
    align_f   = align[align["pair_key_u"].isin(allowed_pair_keys)].copy()

    print(
        f"[allowed pairs] parsed={len(allowed_pair_keys)} | "
        f"tfbs_pairs_present={sim_all_f['pair_key_u'].nunique()} | "
        f"align_pairs_present={align_f['pair_key_u'].nunique()} | "
        f"bad_rows_in_pairs_file={bad}"
    )

    sim_pair = (sim_all_f.groupby("pair_key_u", as_index=False)
                        .agg(tfbs_pct=("jaccard", "mean")))
    sim_pair["tfbs_pct"] *= 100.0

    align_f["_pid"] = pd.to_numeric(align_f["percent_identity"], errors="coerce")
    align_pair = (align_f.groupby("pair_key_u", as_index=False)
                          .agg(seq_pid=("_pid", "mean")))

    df = (sim_pair.merge(align_pair, on="pair_key_u", how="inner")
                 .dropna(subset=["tfbs_pct", "seq_pid"])
                 .copy())

    if df.empty:
        raise ValueError("After filtering + merge, df is empty (pair naming mismatch).")

    x = df["seq_pid"].to_numpy()
    y = df["tfbs_pct"].to_numpy()

    coef = np.polyfit(x, y, deg=deg)
    yhat = np.polyval(coef, x)
    resid = y - yhat

    med = np.median(resid)
    mad = np.median(np.abs(resid - med))
    robust_sd = 1.4826 * mad if mad > 0 else np.std(resid)
    df["tfbs_resid_z"] = resid / robust_sd if robust_sd and np.isfinite(robust_sd) else np.nan

    z = df["tfbs_resid_z"].to_numpy()
    zmin = _round_sig(float(np.nanmin(z)), sig=2)
    zmax = _round_sig(float(np.nanmax(z)), sig=2)

    base = mpl.cm.get_cmap("viridis")
    N = 256
    colors = []
    for i in range(N):
        t = i / (N - 1)
        zz = zmin + t * (zmax - zmin)
        if -deadzone <= zz <= deadzone:
            colors.append((0.70, 0.70, 0.70, 1.0))
        elif zz < -deadzone:
            u = (zz - zmin) / (-deadzone - zmin)
            colors.append(base(0.0 + u * 0.45))
        else:
            u = (zz - deadzone) / (zmax - deadzone)
            colors.append(base(0.55 + u * 0.45))

    cmap_dead = mpl.colors.ListedColormap(colors)
    norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)

    inlier = np.abs(df["tfbs_resid_z"]) <= deadzone
    point_colors = cmap_dead(norm(df["tfbs_resid_z"].to_numpy()))
    point_colors[inlier.to_numpy()] = (0.70, 0.70, 0.70, 0.55)

    FIGSIZE_EVOL = (5.0, 4)  # narrower width, same height as before

    fig, ax = plt.subplots(figsize=FIGSIZE_EVOL)
    ax.scatter(df["seq_pid"], df["tfbs_pct"],
               s=S_POINTS, c=point_colors, edgecolors="none")

    xx = np.linspace(df["seq_pid"].min(), df["seq_pid"].max(), 300)
    yy = np.polyval(coef, xx)
    ax.plot(xx, yy, linewidth=1.5)

    ax.set_xlabel("global sequence identity (%)")
    ax.set_ylabel("TFBS Jaccard similarity (%)")
    enforce_ticks(ax)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_dead)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("TFBS similarity z-score")
    cbar.set_ticks([zmin, -2, -1, 0, 1, 2, zmax])

    # Has colorbar + no title -> COLORBAR preset (reserves right margin)
    save_figure(fig, outbase, layout="COLORBAR", two_panel=False)
    df.to_csv(OUTDIR / f"{outbase}_table.tsv", sep="\t", index=False)
    print(f"zmin={zmin}, zmax={zmax}, n={len(df)}")

# =============================================================================
# RUN
# =============================================================================
if __name__ == "__main__":
    plot_kmer_tfbs_two_panel_violin()
    plot_enh_hit_by_shadow_category(outbase="fig_size_bin_new")
    plot_tfbs_residual_z_scatter_deadzone(outbase="fig_tfbs_resid_z_scatter_deadzone_v4")

import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "Arial",   # use Arial everywhere
    "font.size": 8,           # default text size
    "legend.fontsize": 8,     # legend text size
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 600,
    "savefig.dpi": 600,
})

# import pandas as pd
# import matplotlib.pyplot as plt
#
# # ====== Input file ======
#
#
# # Read data
# df = pd.read_csv(csv_path)
#
# # Make sure column is numeric (in case it was read as text)
# df["dup_enhancers_extra_copies"] = pd.to_numeric(
#     df["dup_enhancers_extra_copies"], errors="coerce"
# )
#
# # Optional: drop rows with missing values in the column
# df = df.dropna(subset=["dup_enhancers_extra_copies"])
#
# # ====== Bin counts from dup_enhancers_extra_copies ======
# n_1   = (df["dup_enhancers_extra_copies"] == 1).sum()
# n_23  = df["dup_enhancers_extra_copies"].between(2, 3, inclusive="both").sum()
# n_ge4 = (df["dup_enhancers_extra_copies"] >= 4).sum()
# TOTAL = n_1 + n_23 + n_ge4
#
# # ====== Appearance (increasing shades of magenta) ======
# C1   = "#F1B6DA"  # light magenta
# C23  = "#C51B7D"  # medium magenta
# CGE4 = "#7A0177"  # dark magenta
#
# DPI = 600
#
# # ====== Build figure ======
# fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=DPI)
#
# values = [n_1, n_23, n_ge4]
# labels = [
#     f"1 hit/set\n{n_1} sets ({100*n_1/TOTAL:.1f}%)",
#     f"2–3 hits/set\n{n_23} sets ({100*n_23/TOTAL:.1f}%)",
#     f"≥4 hits/set\n{n_ge4} sets ({100*n_ge4/TOTAL:.1f}%)",
# ]
# colors = [C1, C23, CGE4]
#
# # Donut plot (labels outside)
# ax.pie(
#     values,
#     labels=labels,
#     colors=colors,
#     startangle=90,
#     counterclock=False,
#     labeldistance=1.15,
#     wedgeprops=dict(width=0.35, edgecolor="white")
# )
#
# ax.set(aspect="equal")
# ax.set_xticks([]); ax.set_yticks([])
# for spine in ax.spines.values():
#     spine.set_visible(False)
#
# plt.tight_layout()
#
# # ====== Save outputs (edit paths as needed) ======
# legacy external save example removed
#             dpi=DPI, bbox_inches="tight")
# legacy external save example removed
#             dpi=300, bbox_inches="tight")
#
# plt.show()

# # Print counts to confirm
# print(f"Total sets: {TOTAL}")
# print(f"1 hit/set: {n_1}")
# print(f"2–3 hits/set: {n_23}")
# print(f"≥4 hits/set: {n_ge4}")
