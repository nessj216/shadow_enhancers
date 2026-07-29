
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from scipy.stats import chisquare, power_divergence, norm
from pathlib import Path

# ----------------------------
# file inputs
# ----------------------------
ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"
OUTPUT = ROOT / "output"
PLOT_FIGURES = ROOT / "plot_figures"
OUTPUT.mkdir(parents=True, exist_ok=True)
PLOT_FIGURES.mkdir(parents=True, exist_ok=True)

genome_file = str(INPUT / "FINAL_TE_genomefile_merged_dedup.bed")

# legacy external path examples removed
single_file = str(INPUT / 'FINAL_overlap_singles_dm6_TEtype.bed')
shadow_file = str(INPUT / 'Full_unfiltered_overlap_shadows_merged.bed')
# legacy external path examples removed
# ----------------------------
# parameters
# ----------------------------
MIN_OVERLAP_BP = 50
FLANK_BP = 30000
EPS = 1e-9

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial"
})

cat_order = ["LTR", "LINE", "SINE", "Helitron/RC", "TIR"]
te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#fdd49e', '#f16913']
PLOT_H = 3.2

AX_LEFT = 0.14
AX_BOTTOM = 0.18
AX_WIDTH = 0.68
AX_HEIGHT = 0.72
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
        return "Helitron/RC"
    if te_type.startswith("DNA"):
        return "TIR"
    return None

def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
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

def plot_stacked_bar(ax, values, xpos, width=0.6):
    bottom = 0
    for val, color in zip(values, te_colors):
        ax.bar(xpos, val, width=width, bottom=bottom, color=color, edgecolor="black")
        if val > 0:
            if val < 6:
                ax.text(xpos + 0.33, bottom + val/2, f"{val:.1f}%", ha="left", va="center", fontsize=10)
            else:
                ax.text(xpos, bottom + val/2, f"{val:.1f}%", ha="center", va="center", fontsize=10)
        bottom += val

def one_proportion_ztest(obs_pos, n, p_ref):
    """
    Two-sided one-sample proportion z-test.
    H0: observed proportion == p_ref
    """
    obs_pos = float(obs_pos)
    n = float(n)
    p_ref = float(p_ref)

    if n <= 0 or pd.isna(p_ref) or p_ref < 0 or p_ref > 1:
        return np.nan, np.nan

    if p_ref in (0, 1):
        return np.nan, np.nan

    p_obs = obs_pos / n
    se = np.sqrt(p_ref * (1 - p_ref) / n)

    if se == 0 or np.isnan(se):
        return np.nan, np.nan

    z = (p_obs - p_ref) / se
    p = 2 * norm.sf(abs(z))
    return z, p

# ----------------------------
# loading
# ----------------------------
def load_background(path):
    bg = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "te_type"])
    bg["start"] = pd.to_numeric(bg["start"], errors="coerce")
    bg["end"] = pd.to_numeric(bg["end"], errors="coerce")
    bg = bg.dropna(subset=["start", "end", "te_type"]).copy()
    bg = bg[bg["end"] > bg["start"]].copy()
    bg["te_bin"] = bg["te_type"].map(map_te_bin)
    bg = bg[bg["te_bin"].notna()].copy()
    return bg

def load_single_overlap(path):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["enh_chrom","enh_start","enh_end","enh_id",
                            "te_chrom","te_start","te_end","te_type","overlap_bp"])
    df["enh_start"] = pd.to_numeric(df["enh_start"], errors="coerce")
    df["enh_end"] = pd.to_numeric(df["enh_end"], errors="coerce")
    df["overlap_bp"] = pd.to_numeric(df["overlap_bp"], errors="coerce")
    df = df.dropna(subset=["enh_start", "enh_end", "te_type", "overlap_bp"]).copy()
    df = df[(df["overlap_bp"] >= MIN_OVERLAP_BP) & (df["te_type"] != ".")].copy()
    df["te_bin"] = df["te_type"].map(map_te_bin)
    return df[df["te_bin"].notna()].copy()

def load_shadow_overlap(path):
    raw = pd.read_csv(path, sep="\t", header=None)
    df = pd.DataFrame({
        "enh_chrom": raw.iloc[:, 0],
        "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "enh_id": raw.iloc[:, 3].astype(str),
        "te_type": raw.iloc[:, 7],
        "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
    })
    df = df.dropna(subset=["enh_start", "enh_end", "te_type", "overlap_bp"]).copy()
    df = df[(df["overlap_bp"] >= MIN_OVERLAP_BP) & (df["te_type"] != ".")].copy()
    df["te_bin"] = df["te_type"].map(map_te_bin)
    return df[df["te_bin"].notna()].copy()

def unique_enhancers(df):
    enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
    enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
    enh["win_end"] = enh["enh_end"] + FLANK_BP
    return enh

# ----------------------------
# observed composition
# enhancer-based category presence
# ----------------------------
def observed_presence_matrix(df, enhancers):
    pairs = df[["enh_id", "te_bin"]].drop_duplicates().copy()
    pairs["value"] = 1
    mat = pairs.pivot_table(index="enh_id", columns="te_bin", values="value",
                            aggfunc="max", fill_value=0)
    mat = mat.reindex(index=enhancers["enh_id"], columns=cat_order, fill_value=0)
    return mat.astype(int)

def observed_overlap_insert_counts(df, enhancers):
    mat = observed_presence_matrix(df, enhancers)
    counts = mat.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)
    return counts

def observed_overlap_insert_composition(df, enhancers):
    counts = observed_overlap_insert_counts(df, enhancers).astype(float)
    return counts / counts.sum() * 100

# ----------------------------
# genome composition
# ----------------------------
def genome_insert_counts(bg):
    return bg["te_bin"].value_counts().reindex(cat_order, fill_value=0).astype(int)

def genome_insert_composition(bg):
    counts = genome_insert_counts(bg).astype(float)
    return counts / counts.sum() * 100

# ----------------------------
# neighborhood composition
# count inserts per enhancer window
# ----------------------------
def neighborhood_insert_counts_per_enhancer(enhancers, bg):
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
                rows.append(pd.Series(0, index=cat_order, name=row["enh_id"]))
                continue

            ov = np.minimum(sub_en, we) - np.maximum(sub_st, ws)
            valid = ov > 0
            sub_bin = sub_bin[valid]

            if len(sub_bin) == 0:
                rows.append(pd.Series(0, index=cat_order, name=row["enh_id"]))
                continue

            counts = pd.Series(sub_bin).value_counts().reindex(cat_order, fill_value=0).astype(int)
            counts.name = row["enh_id"]
            rows.append(counts)

    out = pd.DataFrame(rows)
    out.index = enhancers["enh_id"].tolist()
    out = out.reindex(columns=cat_order, fill_value=0)
    return out.astype(int)

def neighborhood_total_insert_counts(enhancers, bg):
    per_enh = neighborhood_insert_counts_per_enhancer(enhancers, bg)
    total_counts = per_enh.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)
    return total_counts

def mean_local_insert_composition(enhancers, bg):
    per_enh = neighborhood_insert_counts_per_enhancer(enhancers, bg)
    mean_counts = per_enh.mean(axis=0).reindex(cat_order, fill_value=0)
    return mean_counts / mean_counts.sum() * 100

# ----------------------------
# stats
# ----------------------------
def run_gtest_and_chisq(obs_counts, ref_counts, label):
    obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
    ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)

    obs_total = obs_counts.sum()
    ref_props = ref_counts / ref_counts.sum()
    expected = ref_props * obs_total

    keep = ~((obs_counts == 0) & (expected == 0))
    obs_use = obs_counts[keep]
    exp_use = expected[keep]

    if np.any((exp_use == 0) & (obs_use > 0)):
        return pd.DataFrame([{
            "comparison": label,
            "obs_total": int(obs_total),
            "n_categories_tested": int(keep.sum()),
            "chi2_stat": np.nan,
            "chi2_p": np.nan,
            "g_stat": np.nan,
            "g_p": np.nan,
            "note": "Invalid because at least one category has expected=0 but observed>0"
        }])

    chi_stat, chi_p = chisquare(f_obs=obs_use.values, f_exp=exp_use.values)
    g_stat, g_p = power_divergence(f_obs=obs_use.values, f_exp=exp_use.values, lambda_="log-likelihood")

    return pd.DataFrame([{
        "comparison": label,
        "obs_total": int(obs_total),
        "n_categories_tested": int(keep.sum()),
        "chi2_stat": chi_stat,
        "chi2_p": chi_p,
        "g_stat": g_stat,
        "g_p": g_p,
        "note": ""
    }])

def ztest_te_enrichment(obs_presence_mat, ref_counts, group_name):
    """
    For each TE class:
      observed = enhancer presence/absence in the enhancer set
      expected = reference probability from neighborhood composition

    Uses a one-sample proportion z-test.
    """
    n_enh = obs_presence_mat.shape[0]
    obs_pos = obs_presence_mat.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)
    obs_neg = n_enh - obs_pos

    ref_props = pd.Series(ref_counts, index=cat_order, dtype=float)
    ref_props = ref_props / ref_props.sum()

    rows = []
    for cat in cat_order:
        p_ref = float(ref_props[cat])

        zval, pval = one_proportion_ztest(
            obs_pos=obs_pos[cat],
            n=n_enh,
            p_ref=p_ref
        )

        obs_prop = obs_pos[cat] / n_enh if n_enh > 0 else np.nan
        log2_or = np.log2((obs_prop + EPS) / (p_ref + EPS))
        direction = "depleted" if obs_prop < p_ref else "enriched"

        rows.append({
            "group": group_name,
            "te_class": cat,
            "n_enhancers": int(n_enh),
            "obs_pos": int(obs_pos[cat]),
            "obs_neg": int(obs_neg[cat]),
            "ref_prop": p_ref,
            "obs_prop": obs_prop,
            "log2_obs_over_ref": log2_or,
            "z_value": zval,
            "p_value": pval,
            "direction": direction
        })

    return pd.DataFrame(rows)

# ----------------------------
# load data
# ----------------------------
shadow_df = load_shadow_overlap(shadow_file)
single_df = load_single_overlap(single_file)
genome_df = load_background(genome_file)

shadow_enh = unique_enhancers(shadow_df)
single_enh = unique_enhancers(single_df)

shadow_obs_mat = observed_presence_matrix(shadow_df, shadow_enh)
single_obs_mat = observed_presence_matrix(single_df, single_enh)

shadow_obs_counts = observed_overlap_insert_counts(shadow_df, shadow_enh)
single_obs_counts = observed_overlap_insert_counts(single_df, single_enh)

shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)
genome_counts = genome_insert_counts(genome_df)

shadow_obs = observed_overlap_insert_composition(shadow_df, shadow_enh)
single_obs = observed_overlap_insert_composition(single_df, single_enh)
shadow_local = mean_local_insert_composition(shadow_enh, genome_df)
single_local = mean_local_insert_composition(single_enh, genome_df)
genome_comp = genome_insert_composition(genome_df)

summary = pd.DataFrame({
    "shadow_observed_insert_based": shadow_obs,
    "shadow_local_insert_based": shadow_local,
    "single_observed_insert_based": single_obs,
    "single_local_insert_based": single_local,
    "genome_insert_based": genome_comp
}).reindex(cat_order)

print("Unique shadow enhancers:", len(shadow_enh))
print("Unique single enhancers:", len(single_enh))
print("\nInsert-based composition table (%):")
print(summary.round(2))

summary.to_csv(OUTPUT / "stackedplot_local_insert_based_summary.tsv", sep="\t")

# ----------------------------
# overall composition tests for stacked plot
# ----------------------------
stacked_stats = pd.concat([
    run_gtest_and_chisq(shadow_obs_counts, shadow_local_counts, "Shadow observed vs Shadow neighborhood"),
    run_gtest_and_chisq(shadow_obs_counts, genome_counts, "Shadow observed vs Genome"),
    run_gtest_and_chisq(single_obs_counts, single_local_counts, "Single observed vs Single neighborhood"),
    run_gtest_and_chisq(single_obs_counts, genome_counts, "Single observed vs Genome")
], ignore_index=True)

stacked_stats.to_csv(OUTPUT / "stacked_plot_gtest_chisq_stats.tsv", sep="\t", index=False)

print("\n=== Stacked plot overall composition tests ===")
print(stacked_stats.round(6))

# ----------------------------
# Z-tests for log2 bars
# compare observed enhancer presence to neighborhood reference
# ----------------------------
shadow_zstats = ztest_te_enrichment(shadow_obs_mat, shadow_local_counts, "Shadow")
single_zstats = ztest_te_enrichment(single_obs_mat, single_local_counts, "Single")

ztest_stats = pd.concat([shadow_zstats, single_zstats], ignore_index=True)
ztest_stats["significant"] = ztest_stats["p_value"] < 0.05
ztest_stats["depleted_significant"] = (ztest_stats["direction"] == "depleted") & (ztest_stats["p_value"] < 0.05)
ztest_stats["enriched_significant"] = (ztest_stats["direction"] == "enriched") & (ztest_stats["p_value"] < 0.05)
ztest_stats["label"] = ztest_stats["p_value"].map(p_to_star)

ztest_stats.to_csv(OUTPUT / "ztest_TEclass_observed_vs_neighborhood_BH.tsv", sep="\t", index=False)

print("\n=== Z-test TE-class tests ===")
print(ztest_stats.round(6))

# ----------------------------
# plot 1
# ----------------------------
fig = plt.figure(figsize=(9, PLOT_H))
ax = fig.add_axes([AX_LEFT, AX_BOTTOM, AX_WIDTH, AX_HEIGHT])

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
    "singles\nneighborhood",
    "genome"
]

x_positions = [0, 1, 2.4, 3.4, 4.8]

for xpos, vals in zip(x_positions, all_values):
    plot_stacked_bar(ax, vals, xpos, width=0.62)

ax.set_ylabel("TE insert %")
ax.set_xticks(x_positions)
ax.set_xticklabels(all_labels, fontsize=12)
ax.set_ylim(0, 150)

handles = [
    Patch(facecolor=te_colors[0], edgecolor='black', label='LTR'),
    Patch(facecolor=te_colors[1], edgecolor='black', label='LINE'),
    Patch(facecolor=te_colors[2], edgecolor='black', label='SINE'),
    Patch(facecolor=te_colors[3], edgecolor='black', label="Helitron/RC"),
    Patch(facecolor=te_colors[4], edgecolor='black', label='TIR')
]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1), frameon=True)

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

# x-axis ticks black too
ax.tick_params(axis="x", which="both", bottom=True, top=False,
               length=4, width=1, color="black", labelcolor="black")

plt.tight_layout()
plt.savefig(PLOT_FIGURES / "stacked_TE_observed_local_genome_insert_based.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.savefig(PLOT_FIGURES / "stacked_TE_observed_local_genome_insert_based.pdf", bbox_inches='tight', pad_inches=0.1)
plt.show()

print("\n=== Stacked plot overall composition tests ===")
for _, row in stacked_stats.iterrows():
    print(
        f"{row['comparison']}: "
        f"chi2 = {row['chi2_stat']:.4f}, p = {row['chi2_p']:.4g}; "
        f"G = {row['g_stat']:.4f}, p = {row['g_p']:.4g}"
    )

# =========================================================
# LOG2 ENRICHMENT PLOT: observed / neighborhood by TE type
# =========================================================
log2_df = pd.DataFrame({
    "Shadow": np.log2((summary["shadow_observed_insert_based"] + EPS) / (summary["shadow_local_insert_based"] + EPS)),
    "Single": np.log2((summary["single_observed_insert_based"] + EPS) / (summary["single_local_insert_based"] + EPS)),
}, index=summary.index).reindex(cat_order)

print("\nlog2(observed / neighborhood) by TE type:")
print(log2_df.round(3))

log2_df.to_csv(OUTPUT / "log2_observed_over_neighborhood_by_TEtype.tsv", sep="\t")

labels = ['LTR', 'LINE', 'SINE', "Helitron/RC", 'TIR']
te_color_map = {
    'LTR': '#c6dbef',
    'LINE': '#6baed6',
    'SINE': '#1f78b4',
    "Helitron/RC": '#fdd49e',
    'TIR': '#f16913'
}

shadow_vals = log2_df["Shadow"].to_numpy()
single_vals = log2_df["Single"].to_numpy()

x = np.arange(len(cat_order))
bar_width = 0.46

fig = plt.figure(figsize=(6, PLOT_H))
ax = fig.add_axes([AX_LEFT, AX_BOTTOM, AX_WIDTH, AX_HEIGHT])

bars1 = ax.bar(
    x - bar_width/2,
    shadow_vals,
    width=bar_width,
    color=[te_color_map[lbl] for lbl in labels],
    edgecolor='black',
    hatch='//'
)

bars2 = ax.bar(
    x + bar_width/2,
    single_vals,
    width=bar_width,
    color=[te_color_map[lbl] for lbl in labels],
    edgecolor='black'
)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)
ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=12)
ax.tick_params(axis='y', labelsize=11)
ax.tick_params(axis='x', labelsize=11)

#ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
#ax.grid(axis='x', which='both', linestyle='', linewidth=0)
ax.axhline(0, color='black', linewidth=1)
ax.grid(False)
finite_vals = np.concatenate([
    shadow_vals[np.isfinite(shadow_vals)],
    single_vals[np.isfinite(single_vals)]
])
if len(finite_vals) > 0:
    ymin = np.floor(finite_vals.min()) - 0.3
    ymax = np.ceil(finite_vals.max()) + 0.8
    ax.set_ylim(ymin, ymax)
ax.set_ylim(-2.5, 3)

for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        if np.isfinite(height):
            offset = 0.1 if height >= 0 else -0.2
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + offset,
                f'{height:.2f}',
                ha='center',
                va='bottom' if height >= 0 else 'top',
                fontsize=11
            )

# stars from one-sample proportion z-test + BH
for i, cat in enumerate(cat_order):
    sh_row = ztest_stats[(ztest_stats["group"] == "Shadow") & (ztest_stats["te_class"] == cat)].iloc[0]
    si_row = ztest_stats[(ztest_stats["group"] == "Single") & (ztest_stats["te_class"] == cat)].iloc[0]

    if sh_row["label"] != "ns":
        y = shadow_vals[i] + (0.32 if shadow_vals[i] >= 0 else -0.6)
        ax.text(x[i] - bar_width/2, y, sh_row["label"],
                ha='center',
                va='bottom' if shadow_vals[i] >= 0 else 'top',
                fontsize=11, fontweight='bold')

    if si_row["label"] != "ns":
        y = single_vals[i] + (0.32 if single_vals[i] >= 0 else -0.6)
        ax.text(x[i] + bar_width/2, y, si_row["label"],
                ha='center',
                va='bottom' if single_vals[i] >= 0 else 'top',
                fontsize=11, fontweight='bold')

shadow_patch = mpatches.Patch(facecolor='white', edgecolor='black', hatch='//', label='shadows')
single_patch = mpatches.Patch(facecolor='white', edgecolor='black', label='singles')
ax.legend(handles=[shadow_patch, single_patch], fontsize=11, loc='upper right')
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

# x-axis ticks black too
ax.tick_params(axis="x", which="both", bottom=True, top=False,
               length=4, width=1, color="black", labelcolor="black")
plt.tight_layout()
plt.savefig(PLOT_FIGURES / "log2_observed_over_neighborhood_by_TEtype.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.savefig(PLOT_FIGURES / "log2_observed_over_neighborhood_by_TEtype.pdf", bbox_inches='tight', pad_inches=0.1)
plt.show()







# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from matplotlib.patches import Patch
# from scipy.stats import chisquare, power_divergence, fisher_exact
#
# # ----------------------------
# # file inputs
# # ----------------------------
# legacy external path examples removed
#
# legacy external path examples removed
#
# # ----------------------------
# # parameters
# # ----------------------------
# MIN_OVERLAP_BP = 50
# FLANK_BP = 30000
# EPS = 1e-9
#
# plt.style.use("seaborn-v0_8-whitegrid")
# plt.rcParams.update({
#     "font.size": 12,
#     "xtick.labelsize": 12,
#     "ytick.labelsize": 12,
#     "font.family": "Arial"
# })
#
# cat_order = ["LTR", "LINE", "SINE", "DNA/RC", "DNA/other"]
# te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#fdd49e', '#f16913']
#
#
# # ----------------------------
# # helpers
# # ----------------------------
# def map_te_bin(te_type):
#     te_type = str(te_type)
#     if te_type.startswith("LTR"):
#         return "LTR"
#     if te_type.startswith("LINE"):
#         return "LINE"
#     if te_type.startswith("SINE"):
#         return "SINE"
#     if te_type == "DNA/RC":
#         return "DNA/RC"
#     if te_type.startswith("DNA"):
#         return "DNA/other"
#     return None
#
# def bh_fdr(pvals):
#     pvals = np.asarray(pvals, dtype=float)
#     n = len(pvals)
#     order = np.argsort(pvals)
#     ranked = pvals[order]
#     q = ranked * n / np.arange(1, n + 1)
#     q = np.minimum.accumulate(q[::-1])[::-1]
#     q = np.clip(q, 0, 1)
#     out = np.empty(n)
#     out[order] = q
#     return out
#
# def p_to_star(p):
#     if p < 0.001:
#         return "***"
#     if p < 0.01:
#         return "**"
#     if p < 0.05:
#         return "*"
#     return "ns"
#
# def plot_stacked_bar(ax, values, xpos, width=0.6):
#     bottom = 0
#     for val, color in zip(values, te_colors):
#         ax.bar(xpos, val, width=width, bottom=bottom, color=color, edgecolor="black")
#         if val > 0:
#             if val < 6:
#                 ax.text(xpos + 0.33, bottom + val/2, f"{val:.1f}%", ha="left", va="center", fontsize=10)
#             else:
#                 ax.text(xpos, bottom + val/2, f"{val:.1f}%", ha="center", va="center", fontsize=10)
#         bottom += val
#
#
# # ----------------------------
# # loading
# # ----------------------------
# def load_background(path):
#     bg = pd.read_csv(path, sep="\t", header=None,
#                      names=["chrom", "start", "end", "te_type"])
#     bg["start"] = pd.to_numeric(bg["start"], errors="coerce")
#     bg["end"] = pd.to_numeric(bg["end"], errors="coerce")
#     bg = bg.dropna(subset=["start", "end", "te_type"]).copy()
#     bg = bg[bg["end"] > bg["start"]].copy()
#     bg["te_bin"] = bg["te_type"].map(map_te_bin)
#     bg = bg[bg["te_bin"].notna()].copy()
#     return bg
#
# def load_single_overlap(path):
#     df = pd.read_csv(path, sep="\t", header=None,
#                      names=["enh_chrom","enh_start","enh_end","enh_id",
#                             "te_chrom","te_start","te_end","te_type","overlap_bp"])
#     df["enh_start"] = pd.to_numeric(df["enh_start"], errors="coerce")
#     df["enh_end"] = pd.to_numeric(df["enh_end"], errors="coerce")
#     df["overlap_bp"] = pd.to_numeric(df["overlap_bp"], errors="coerce")
#     df = df.dropna(subset=["enh_start", "enh_end", "te_type", "overlap_bp"]).copy()
#     df = df[(df["overlap_bp"] >= MIN_OVERLAP_BP) & (df["te_type"] != ".")].copy()
#     df["te_bin"] = df["te_type"].map(map_te_bin)
#     return df[df["te_bin"].notna()].copy()
#
# def load_shadow_overlap(path):
#     raw = pd.read_csv(path, sep="\t", header=None)
#     df = pd.DataFrame({
#         "enh_chrom": raw.iloc[:, 0],
#         "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
#         "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
#         "enh_id": raw.iloc[:, 3].astype(str),
#         "te_type": raw.iloc[:, 7],
#         "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
#     })
#     df = df.dropna(subset=["enh_start", "enh_end", "te_type", "overlap_bp"]).copy()
#     df = df[(df["overlap_bp"] >= MIN_OVERLAP_BP) & (df["te_type"] != ".")].copy()
#     df["te_bin"] = df["te_type"].map(map_te_bin)
#     return df[df["te_bin"].notna()].copy()
#
# def unique_enhancers(df):
#     enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
#     enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
#     enh["win_end"] = enh["enh_end"] + FLANK_BP
#     return enh
#
#
# # ----------------------------
# # observed composition
# # enhancer-based category presence
# # ----------------------------
# def observed_presence_matrix(df, enhancers):
#     pairs = df[["enh_id", "te_bin"]].drop_duplicates().copy()
#     pairs["value"] = 1
#     mat = pairs.pivot_table(index="enh_id", columns="te_bin", values="value",
#                             aggfunc="max", fill_value=0)
#     mat = mat.reindex(index=enhancers["enh_id"], columns=cat_order, fill_value=0)
#     return mat.astype(int)
#
# def observed_overlap_insert_counts(df, enhancers):
#     mat = observed_presence_matrix(df, enhancers)
#     counts = mat.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)
#     return counts
#
# def observed_overlap_insert_composition(df, enhancers):
#     counts = observed_overlap_insert_counts(df, enhancers).astype(float)
#     return counts / counts.sum() * 100
#
#
# # ----------------------------
# # genome composition
# # ----------------------------
# def genome_insert_counts(bg):
#     return bg["te_bin"].value_counts().reindex(cat_order, fill_value=0).astype(int)
#
# def genome_insert_composition(bg):
#     counts = genome_insert_counts(bg).astype(float)
#     return counts / counts.sum() * 100
#
#
# # ----------------------------
# # neighborhood composition
# # count inserts per enhancer window
# # ----------------------------
# def neighborhood_insert_counts_per_enhancer(enhancers, bg):
#     rows = []
#
#     for chrom, enh_chr in enhancers.groupby("enh_chrom"):
#         te_chr = bg[bg["chrom"] == chrom].copy()
#         te_st = te_chr["start"].to_numpy()
#         te_en = te_chr["end"].to_numpy()
#         te_bin = te_chr["te_bin"].to_numpy()
#
#         for _, row in enh_chr.iterrows():
#             ws = row["win_start"]
#             we = row["win_end"]
#
#             mask = (te_en > ws) & (te_st < we)
#             sub_st = te_st[mask]
#             sub_en = te_en[mask]
#             sub_bin = te_bin[mask]
#
#             if len(sub_st) == 0:
#                 rows.append(pd.Series(0, index=cat_order, name=row["enh_id"]))
#                 continue
#
#             ov = np.minimum(sub_en, we) - np.maximum(sub_st, ws)
#             valid = ov > 0
#             sub_bin = sub_bin[valid]
#
#             if len(sub_bin) == 0:
#                 rows.append(pd.Series(0, index=cat_order, name=row["enh_id"]))
#                 continue
#
#             counts = pd.Series(sub_bin).value_counts().reindex(cat_order, fill_value=0).astype(int)
#             counts.name = row["enh_id"]
#             rows.append(counts)
#
#     out = pd.DataFrame(rows)
#     out.index = enhancers["enh_id"].tolist()
#     out = out.reindex(columns=cat_order, fill_value=0)
#     return out.astype(int)
#
# def neighborhood_total_insert_counts(enhancers, bg):
#     per_enh = neighborhood_insert_counts_per_enhancer(enhancers, bg)
#     total_counts = per_enh.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)
#     return total_counts
#
# def mean_local_insert_composition(enhancers, bg):
#     per_enh = neighborhood_insert_counts_per_enhancer(enhancers, bg)
#     mean_counts = per_enh.mean(axis=0).reindex(cat_order, fill_value=0)
#     return mean_counts / mean_counts.sum() * 100
#
#
# # ----------------------------
# # stats
# # ----------------------------
# def run_gtest_and_chisq(obs_counts, ref_counts, label):
#     obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
#     ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)
#
#     obs_total = obs_counts.sum()
#     ref_props = ref_counts / ref_counts.sum()
#     expected = ref_props * obs_total
#
#     # drop categories that are 0 in both observed and expected
#     keep = ~((obs_counts == 0) & (expected == 0))
#     obs_use = obs_counts[keep]
#     exp_use = expected[keep]
#
#     # if expected is still 0 anywhere but observed > 0, test is invalid as-is
#     if np.any((exp_use == 0) & (obs_use > 0)):
#         return pd.DataFrame([{
#             "comparison": label,
#             "obs_total": int(obs_total),
#             "n_categories_tested": int(keep.sum()),
#             "chi2_stat": np.nan,
#             "chi2_p": np.nan,
#             "g_stat": np.nan,
#             "g_p": np.nan,
#             "note": "Invalid because at least one category has expected=0 but observed>0"
#         }])
#
#     chi_stat, chi_p = chisquare(f_obs=obs_use.values, f_exp=exp_use.values)
#     g_stat, g_p = power_divergence(f_obs=obs_use.values, f_exp=exp_use.values, lambda_="log-likelihood")
#
#     return pd.DataFrame([{
#         "comparison": label,
#         "obs_total": int(obs_total),
#         "n_categories_tested": int(keep.sum()),
#         "chi2_stat": chi_stat,
#         "chi2_p": chi_p,
#         "g_stat": g_stat,
#         "g_p": g_p,
#         "note": ""
#     }])
#
# def fisher_te_enrichment(obs_presence_mat, ref_counts, group_name):
#     """
#     For each TE class:
#       observed = enhancer presence/absence in the enhancer set
#       expected = reference probability from neighborhood composition
#     Uses a 2x2 Fisher test on integer pseudo-counts scaled to same total N.
#     """
#     n_enh = obs_presence_mat.shape[0]
#     obs_pos = obs_presence_mat.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)
#     obs_neg = n_enh - obs_pos
#
#     ref_props = pd.Series(ref_counts, index=cat_order, dtype=float)
#     ref_props = ref_props / ref_props.sum()
#
#     rows = []
#     for cat in cat_order:
#         p_ref = ref_props[cat]
#
#         # expected reference counts, scaled to same N as enhancer count
#         ref_pos = int(round(p_ref * n_enh))
#         ref_pos = min(max(ref_pos, 0), n_enh)
#         ref_neg = n_enh - ref_pos
#
#         table = np.array([
#             [int(obs_pos[cat]), int(obs_neg[cat])],
#             [int(ref_pos), int(ref_neg)]
#         ])
#
#         oddsratio, pval = fisher_exact(table, alternative="two-sided")
#
#         obs_prop = obs_pos[cat] / n_enh if n_enh > 0 else np.nan
#         ref_prop_val = ref_pos / n_enh if n_enh > 0 else np.nan
#         log2_or = np.log2((obs_prop + EPS) / (ref_prop_val + EPS))
#
#         direction = "depleted" if obs_prop < ref_prop_val else "enriched"
#
#         rows.append({
#             "group": group_name,
#             "te_class": cat,
#             "n_enhancers": n_enh,
#             "obs_pos": int(obs_pos[cat]),
#             "obs_neg": int(obs_neg[cat]),
#             "ref_pos_scaled": int(ref_pos),
#             "ref_neg_scaled": int(ref_neg),
#             "obs_prop": obs_prop,
#             "ref_prop_scaled": ref_prop_val,
#             "log2_obs_over_ref": log2_or,
#             "oddsratio": oddsratio,
#             "p_value": pval,
#             "direction": direction
#         })
#
#     out = pd.DataFrame(rows)
#     return out
#
#
# # ----------------------------
# # load data
# # ----------------------------
# shadow_df = load_shadow_overlap(shadow_file)
# single_df = load_single_overlap(single_file)
# genome_df = load_background(genome_file)
#
# shadow_enh = unique_enhancers(shadow_df)
# single_enh = unique_enhancers(single_df)
#
# shadow_obs_mat = observed_presence_matrix(shadow_df, shadow_enh)
# single_obs_mat = observed_presence_matrix(single_df, single_enh)
#
# shadow_obs_counts = observed_overlap_insert_counts(shadow_df, shadow_enh)
# single_obs_counts = observed_overlap_insert_counts(single_df, single_enh)
#
# shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
# single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)
# genome_counts = genome_insert_counts(genome_df)
#
# shadow_obs = observed_overlap_insert_composition(shadow_df, shadow_enh)
# single_obs = observed_overlap_insert_composition(single_df, single_enh)
# shadow_local = mean_local_insert_composition(shadow_enh, genome_df)
# single_local = mean_local_insert_composition(single_enh, genome_df)
# genome_comp = genome_insert_composition(genome_df)
#
# summary = pd.DataFrame({
#     "shadow_observed_insert_based": shadow_obs,
#     "shadow_local_insert_based": shadow_local,
#     "single_observed_insert_based": single_obs,
#     "single_local_insert_based": single_local,
#     "genome_insert_based": genome_comp
# }).reindex(cat_order)
#
# print("Unique shadow enhancers:", len(shadow_enh))
# print("Unique single enhancers:", len(single_enh))
# print("\nInsert-based composition table (%):")
# print(summary.round(2))
#
# summary.to_csv("stackedplot_local_insert_based_summary.tsv", sep="\t")
#
#
# # ----------------------------
# # overall composition tests for stacked plot
# # ----------------------------
# stacked_stats = pd.concat([
#     run_gtest_and_chisq(shadow_obs_counts, shadow_local_counts, "Shadow observed vs Shadow neighborhood"),
#     run_gtest_and_chisq(shadow_obs_counts, genome_counts, "Shadow observed vs Genome"),
#     run_gtest_and_chisq(single_obs_counts, single_local_counts, "Single observed vs Single neighborhood"),
#     run_gtest_and_chisq(single_obs_counts, genome_counts, "Single observed vs Genome")
# ], ignore_index=True)
#
# stacked_stats.to_csv("stacked_plot_gtest_chisq_stats.tsv", sep="\t", index=False)
#
# print("\n=== Stacked plot overall composition tests ===")
# print(stacked_stats.round(6))
#
#
# # ----------------------------
# # Fisher tests for low/high log2 bars
# # compare observed enhancer presence to neighborhood reference
# # ----------------------------
# shadow_fisher = fisher_te_enrichment(shadow_obs_mat, shadow_local_counts, "Shadow")
# single_fisher = fisher_te_enrichment(single_obs_mat, single_local_counts, "Single")
#
# fisher_stats = pd.concat([shadow_fisher, single_fisher], ignore_index=True)
# fisher_stats["q_BH"] = bh_fdr(fisher_stats["p_value"].values)
# fisher_stats["significant"] = fisher_stats["q_BH"] < 0.05
# fisher_stats["depleted_significant"] = (fisher_stats["direction"] == "depleted") & (fisher_stats["q_BH"] < 0.05)
# fisher_stats["enriched_significant"] = (fisher_stats["direction"] == "enriched") & (fisher_stats["q_BH"] < 0.05)
# fisher_stats["label"] = fisher_stats["q_BH"].map(p_to_star)
#
# fisher_stats.to_csv("fisher_TEclass_observed_vs_neighborhood_BH.tsv", sep="\t", index=False)
#
# print("\n=== Fisher TE-class tests ===")
# print(fisher_stats.round(6))
#
#
# # ----------------------------
# # plot 1
# # ----------------------------
# fig, ax = plt.subplots(figsize=(8, 3))
# ax.grid(False)
#
# all_values = [
#     shadow_obs.values,
#     shadow_local.values,
#     single_obs.values,
#     single_local.values,
#     genome_comp.values
# ]
#
# all_labels = [
#     "shadows",
#     "shadow\nneighborhood",
#     "singles",
#     "single\nneighborhood",
#     "genome"
# ]
#
# x_positions = [0, 1, 2.4, 3.4, 4.8]
#
# for xpos, vals in zip(x_positions, all_values):
#     plot_stacked_bar(ax, vals, xpos, width=0.62)
#
# ax.set_ylabel("TE insert percentage")
# ax.set_xticks(x_positions)
# ax.set_xticklabels(all_labels, fontsize=12)
# ax.set_ylim(0, 150)
#
# handles = [
#     Patch(facecolor=te_colors[0], edgecolor='black', label='LTR'),
#     Patch(facecolor=te_colors[1], edgecolor='black', label='LINE'),
#     Patch(facecolor=te_colors[2], edgecolor='black', label='SINE'),
#     Patch(facecolor=te_colors[3], edgecolor='black', label='DNA/RC'),
#     Patch(facecolor=te_colors[4], edgecolor='black', label='DNA/other')
# ]
# ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1), frameon=True)
#
# # annotate overall p-values
# sh_local_p = stacked_stats.loc[
#     stacked_stats["comparison"] == "Shadow observed vs Shadow neighborhood", "g_p"
# ].iloc[0]
# sh_genome_p = stacked_stats.loc[
#     stacked_stats["comparison"] == "Shadow observed vs Genome", "g_p"
# ].iloc[0]
# si_local_p = stacked_stats.loc[
#     stacked_stats["comparison"] == "Single observed vs Single neighborhood", "g_p"
# ].iloc[0]
# si_genome_p = stacked_stats.loc[
#     stacked_stats["comparison"] == "Single observed vs Genome", "g_p"
# ].iloc[0]
#
# # ax.text(0.5, 112, f"Shadow\nvs neigh G p={sh_local_p:.3g}\nvs genome G p={sh_genome_p:.3g}",
# #         ha="center", va="bottom", fontsize=10)
# # ax.text(2.9, 112, f"Single\nvs neigh G p={si_local_p:.3g}\nvs genome G p={si_genome_p:.3g}",
# #         ha="center", va="bottom", fontsize=10)
#
# for spine in ax.spines.values():
#     spine.set_visible(True)
#
# plt.tight_layout()
# legacy external save examples removed
# plt.show()
#
# print("\n=== Stacked plot overall composition tests ===")
# for _, row in stacked_stats.iterrows():
#     print(
#         f"{row['comparison']}: "
#         f"chi2 = {row['chi2_stat']:.4f}, p = {row['chi2_p']:.4g}; "
#         f"G = {row['g_stat']:.4f}, p = {row['g_p']:.4g}"
#     )
# # =========================================================
# # LOG2 ENRICHMENT PLOT: observed / neighborhood by TE type
# # =========================================================
# log2_df = pd.DataFrame({
#     "Shadow": np.log2((summary["shadow_observed_insert_based"] + EPS) / (summary["shadow_local_insert_based"] + EPS)),
#     "Single": np.log2((summary["single_observed_insert_based"] + EPS) / (summary["single_local_insert_based"] + EPS)),
# }, index=summary.index).reindex(cat_order)
#
# print("\nlog2(observed / neighborhood) by TE type:")
# print(log2_df.round(3))
#
# log2_df.to_csv("log2_observed_over_neighborhood_by_TEtype.tsv", sep="\t")
#
# labels = ['LTR', 'LINE', 'SINE', 'DNA/RC', 'DNA/other']
# te_color_map = {
#     'LTR': '#c6dbef',
#     'LINE': '#6baed6',
#     'SINE': '#1f78b4',
#     'DNA/RC': '#fdd49e',
#     'DNA/other': '#f16913'
# }
#
# shadow_vals = log2_df["Shadow"].to_numpy()
# single_vals = log2_df["Single"].to_numpy()
#
# x = np.arange(len(cat_order))
# bar_width = 0.46
#
# fig, ax = plt.subplots(figsize=(5, 3))
#
# bars1 = ax.bar(
#     x - bar_width/2,
#     shadow_vals,
#     width=bar_width,
#     color=[te_color_map[lbl] for lbl in labels],
#     edgecolor='black',
#     hatch='//'
# )
#
# bars2 = ax.bar(
#     x + bar_width/2,
#     single_vals,
#     width=bar_width,
#     color=[te_color_map[lbl] for lbl in labels],
#     edgecolor='black'
# )
#
# ax.set_xticks(x)
# ax.set_xticklabels(labels, fontsize=12)
# ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=12)
# ax.tick_params(axis='y', labelsize=11)
# ax.tick_params(axis='x', labelsize=11)
#
# ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
# ax.grid(axis='x', which='both', linestyle='', linewidth=0)
# ax.axhline(0, color='black', linewidth=1)
#
# finite_vals = np.concatenate([
#     shadow_vals[np.isfinite(shadow_vals)],
#     single_vals[np.isfinite(single_vals)]
# ])
# if len(finite_vals) > 0:
#     ymin = np.floor(finite_vals.min()) - 0.3
#     ymax = np.ceil(finite_vals.max()) + 0.8
#     ax.set_ylim(ymin, ymax)
# ax.set_ylim(-2.5, 3)
# # numeric labels
# for bars in [bars1, bars2]:
#     for bar in bars:
#         height = bar.get_height()
#         if np.isfinite(height):
#             offset = 0.1 if height >= 0 else -0.2
#             ax.text(
#                 bar.get_x() + bar.get_width() / 2,
#                 height + offset,
#                 f'{height:.2f}',
#                 ha='center',
#                 va='bottom' if height >= 0 else 'top',
#                 fontsize=11
#             )
#
# # stars from Fisher + BH
# for i, cat in enumerate(cat_order):
#     sh_row = fisher_stats[(fisher_stats["group"] == "Shadow") & (fisher_stats["te_class"] == cat)].iloc[0]
#     si_row = fisher_stats[(fisher_stats["group"] == "Single") & (fisher_stats["te_class"] == cat)].iloc[0]
#
#     if sh_row["label"] != "ns":
#         y = shadow_vals[i] + (0.32 if shadow_vals[i] >= 0 else -0.6)
#         ax.text(x[i] - bar_width/2, y, sh_row["label"],
#                 ha='center',
#                 va='bottom' if shadow_vals[i] >= 0 else 'top',
#                 fontsize=11, fontweight='bold')
#
#     if si_row["label"] != "ns":
#         y = single_vals[i] + (0.32 if single_vals[i] >= 0 else -0.6)
#         ax.text(x[i] + bar_width/2, y, si_row["label"],
#                 ha='center',
#                 va='bottom' if single_vals[i] >= 0 else 'top',
#                 fontsize=11, fontweight='bold')
#
# shadow_patch = mpatches.Patch(facecolor='white', edgecolor='black', hatch='//', label='Shadow enhancers')
# single_patch = mpatches.Patch(facecolor='white', edgecolor='black', label='Single enhancers')
# ax.legend(handles=[shadow_patch, single_patch], fontsize=11, loc='upper right')
#
# plt.tight_layout()
# legacy external save examples removed
# plt.show()
#
# =========================================================
# DOT PLOT: pooled fly shadows vs singles, no tissue split
# =========================================================

from matplotlib.lines import Line2D

# Use the ztest_stats table you already made above.
# It has:
#   group = Shadow / Single
#   te_class
#   log2_obs_over_ref
#   label
#
# This follows the no-tissue dotplot style from the mouse script,
# but only plots pooled Shadow and Single points.

dotplot_df = ztest_stats.copy()
dotplot_df["group"] = dotplot_df["group"].astype(str)

print("\nDot plot values:")
print(dotplot_df[["group", "te_class", "log2_obs_over_ref", "label", "direction"]].round(4))

dotplot_df.to_csv(
    OUTPUT / "fly_pooled_shadow_single_dotplot_values.tsv",
    sep="\t",
    index=False
)

# ----------------------------
# plot settings
# ----------------------------
DOTPLOT_FIGSIZE = (6.5, 3.2)
DOTPLOT_YLIM = (-2.5, 3.0)

x_base = np.arange(len(cat_order)) * 1.45

group_xoffset = {
    "Shadow": -0.13,
    "Single": 0.13
}

star_offset = 0.16

fig, ax = plt.subplots(figsize=DOTPLOT_FIGSIZE)

# zero line
ax.axhline(0, color="black", linewidth=1.1, zorder=1)

# light separators between TE classes
for x_sep in (x_base[:-1] + x_base[1:]) / 2:
    ax.axvline(x_sep, color="0.90", linewidth=0.8, zorder=0)

# ----------------------------
# points
# ----------------------------
for te_i, te in enumerate(cat_order):
    sh = dotplot_df[
        (dotplot_df["group"] == "Shadow") &
        (dotplot_df["te_class"] == te)
    ]

    si = dotplot_df[
        (dotplot_df["group"] == "Single") &
        (dotplot_df["te_class"] == te)
    ]

    if sh.empty or si.empty:
        continue

    x_sh = x_base[te_i] + group_xoffset["Shadow"]
    x_si = x_base[te_i] + group_xoffset["Single"]

    y_sh = float(sh.iloc[0]["log2_obs_over_ref"])
    y_si = float(si.iloc[0]["log2_obs_over_ref"])

    # Shadow = filled black circle
    ax.scatter(
        x_sh,
        y_sh,
        s=72,
        facecolor="black",
        edgecolor="black",
        linewidth=1.0,
        zorder=4
    )

    # Single = open circle
    ax.scatter(
        x_si,
        y_si,
        s=72,
        facecolor="white",
        edgecolor="black",
        linewidth=1.8,
        zorder=4
    )

    # Significance stars
    for x_pt, y_pt, row in [
        (x_sh, y_sh, sh.iloc[0]),
        (x_si, y_si, si.iloc[0])
    ]:
        label = str(row["label"])

        if label != "ns":
            ax.text(
                x_pt,
                y_pt + star_offset,
                label,
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
                zorder=6,
                clip_on=False
            )

# ----------------------------
# axes
# ----------------------------
ax.set_xticks(x_base)
ax.set_xticklabels(cat_order, fontsize=13)
ax.set_xlim(x_base[0] - 0.55, x_base[-1] + 0.55)

ax.set_ylabel(
    r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$',
    fontsize=13
)

ax.set_ylim(*DOTPLOT_YLIM)
ax.grid(False)

# clean axes like mouse plot
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
    labelcolor="black"
)

# ----------------------------
# legend
# ----------------------------
legend_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="black",
        markeredgecolor="black",
        markersize=7,
        label="Shadows"
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="white",
        markeredgecolor="black",
        markeredgewidth=1.8,
        markersize=7,
        label="Singles"
    )
]

ax.legend(
    handles=legend_handles,
    title="Enhancer class",
    frameon=False,
    fontsize=10,
    title_fontsize=11,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.00)
)

plt.tight_layout()

plt.savefig(
    PLOT_FIGURES / "fly_pooled_shadow_single_dotplot_log2_observed_over_neighborhood.png",
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.1
)

plt.savefig(
    PLOT_FIGURES / "fly_pooled_shadow_single_dotplot_log2_observed_over_neighborhood.pdf",
    bbox_inches="tight",
    pad_inches=0.1
)

plt.close(fig)
