# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from matplotlib.lines import Line2D
# from matplotlib.offsetbox import AnchoredOffsetbox, VPacker, HPacker, TextArea, DrawingArea
# import matplotlib.patches as patches
# from scipy.stats import chisquare, power_divergence, chi2_contingency, chi2 as chi2_dist
#
#
# # ----------------------------
# # file inputs
# # ----------------------------
# genome_file = "/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/filtered_mm10.fa.bed"
# shadow_file = "/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/new_20percent/Shadow_Teoverlap(1)_lastcol_ge50_overlapAtLeast20pctOfEnhancer.bed"
# single_file = "/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/new_20percent/singlesmouse_TEcooption.bed"
# output_dir = "/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/figures"
#
#
# # ----------------------------
# # parameters
# # ----------------------------
# MIN_OVERLAP_BP = 50
# FLANK_BP = 30000
# EPS = 1e-9
# PLOT_METRIC = "log2_obs_over_ref"   # or "log2_odds_ratio"
#
# cat_order = ["LTR", "LINE", "SINE", "TIR"]
# tissue_order = ["forebrain", "heart", "limb"]
# te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#f16913']
# te_color_map = {
#     "LTR": "#c6dbef",
#     "LINE": "#6baed6",
#     "SINE": "#1f78b4",
#     "TIR": "#f16913"
# }
# tissue_colors = {
#     "forebrain": "#4C78A8",
#     "heart": "#E45756",
#     "limb": "#54A24B",
#     "pooled": "black"
# }
#
# # Shared physical height so stacked plots and tissue enrichment plots match.
# PANEL_HEIGHT = 3.0
# STACKED_FIGSIZE = (8.0, PANEL_HEIGHT)
# ENRICHMENT_BAR_FIGSIZE = (5.0, PANEL_HEIGHT)
# TISSUE_SPLIT_FIGSIZE = (8.3, PANEL_HEIGHT)
# COMBINED_DOTPLOT_FIGSIZE = (7.8, PANEL_HEIGHT)
#
# # Shared enrichment y-scale so all enrichment-style plots use the same vertical range.
# ENRICHMENT_YLIM = (-3, 3.2)
#
# plt.style.use("seaborn-v0_8-whitegrid")
# plt.rcParams.update({
#     "font.size": 12,
#     "xtick.labelsize": 12,
#     "ytick.labelsize": 12,
#     "font.family": "Arial"
# })
#
#
# # ----------------------------
# # helpers
# # ----------------------------
# def map_te_bin(te_type):
#     te_type = str(te_type)
#
#     if te_type.startswith("LTR"):
#         return "LTR"
#     if te_type.startswith("LINE"):
#         return "LINE"
#     if te_type.startswith("SINE"):
#         return "SINE"
#     if te_type.startswith("DNA"):
#         return "TIR"
#     return None
#
#
# def assign_tissue_from_enh_id(enh_id):
#     s = str(enh_id).lower()
#
#     if "forebrain" in s or "fb" in s:
#         return "forebrain"
#     if "heart" in s:
#         return "heart"
#     if "limb" in s:
#         return "limb"
#     return "other"
#
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
#
# def p_to_star(p):
#     if pd.isna(p):
#         return "ns"
#     if p < 0.001:
#         return "***"
#     if p < 0.01:
#         return "**"
#     if p < 0.05:
#         return "*"
#     return "ns"
#
#
# def style_axes(ax):
#     ax.grid(False)
#     ax.spines["top"].set_visible(False)
#     ax.spines["right"].set_visible(False)
#     ax.spines["left"].set_visible(True)
#     ax.spines["bottom"].set_visible(True)
#     ax.spines["left"].set_color("black")
#     ax.spines["bottom"].set_color("black")
#     ax.spines["left"].set_linewidth(1.1)
#     ax.spines["bottom"].set_linewidth(1.1)
#     ax.tick_params(
#         axis="both",
#         which="both",
#         left=True,
#         right=False,
#         bottom=True,
#         top=False,
#         length=4,
#         width=1.0,
#         color="black",
#         labelcolor="black"
#     )
#
#
# def save_figure(fig, png_path=None, pdf_path=None, left=0.10, right=0.82, bottom=0.22, top=0.92):
#     fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
#     if png_path:
#         fig.savefig(png_path, dpi=600)
#     if pdf_path:
#         fig.savefig(pdf_path)
#     plt.close(fig)
#
#
# def make_grouped_te_legend(ax):
#     box_w = 34
#     box_text_sep = 6
#
#     def legend_row(color, label):
#         da = DrawingArea(box_w, 14, 0, 0)
#         rect = patches.Rectangle(
#             (0, 2), 30, 10,
#             facecolor=color,
#             edgecolor="black",
#             linewidth=1.2
#         )
#         da.add_artist(rect)
#         txt = TextArea(label, textprops=dict(size=12, family="Arial"))
#         return HPacker(children=[da, txt], align="center", pad=0, sep=box_text_sep)
#
#     def header_row(label):
#         spacer = DrawingArea(box_w + box_text_sep, 1, 0, 0)
#         txt = TextArea(label, textprops=dict(size=14, weight="bold", family="Arial"))
#         return HPacker(children=[spacer, txt], align="center", pad=0, sep=0)
#
#     legend_box = VPacker(
#         children=[
#             header_row("RNA"),
#             legend_row(te_colors[0], "LTR"),
#             legend_row(te_colors[1], "LINE"),
#             legend_row(te_colors[2], "SINE"),
#             header_row("DNA"),
#             legend_row(te_colors[3], "TIR"),
#         ],
#         align="left",
#         pad=0,
#         sep=4
#     )
#
#     anchored_box = AnchoredOffsetbox(
#         loc="upper left",
#         child=legend_box,
#         pad=0.3,
#         frameon=False,
#         bbox_to_anchor=(1.02, 1.0),
#         bbox_transform=ax.transAxes,
#         borderpad=0
#     )
#     ax.add_artist(anchored_box)
#
#
# def plot_stacked_bar(ax, values, xpos, width=0.62):
#     bottom = 0
#     for val, color in zip(values, te_colors):
#         ax.bar(xpos, val, width=width, bottom=bottom, color=color, edgecolor="black")
#         if val > 0:
#             if val < 6:
#                 ax.text(xpos + 0.33, bottom + val / 2, f"{val:.1f}%", ha="left", va="center", fontsize=11)
#             else:
#                 ax.text(xpos, bottom + val / 2, f"{val:.1f}%", ha="center", va="center", fontsize=11)
#         bottom += val
#
#
# def apply_enrichment_axis(ax):
#     ax.set_ylim(*ENRICHMENT_YLIM)
#     ax.axhline(0, color="black", linewidth=1)
#
# def composition_effect_size(obs_counts, ref_counts, label):
#     obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
#     ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)
#
#     obs_total = obs_counts.sum()
#     ref_total = ref_counts.sum()
#
#     obs_props = obs_counts / obs_total
#     ref_props = ref_counts / ref_total
#
#     keep = ref_props > 0
#     obs_props = obs_props[keep]
#     ref_props = ref_props[keep]
#
#     prop_chi2 = (((obs_props - ref_props) ** 2) / ref_props).sum()
#     cohens_w = np.sqrt(prop_chi2)
#
#     return pd.DataFrame([{
#         "comparison": label,
#         "obs_total": obs_total,
#         "ref_total": ref_total,
#         "proportional_chi2": prop_chi2,
#         "cohens_w": cohens_w,
#         "note": "Effect size from proportions only; not a chi-square p-value"
#     }])
# # ----------------------------
# # loading
# # ----------------------------
# def load_background(path):
#     raw = pd.read_csv(path, sep="\t", header=None)
#
#     if raw.shape[1] >= 5:
#         bg = pd.DataFrame({
#             "chrom": raw.iloc[:, 0],
#             "start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
#             "end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
#             "te_name": raw.iloc[:, 3].astype(str),
#             "te_type": raw.iloc[:, 4],
#         })
#     elif raw.shape[1] == 4:
#         bg = pd.DataFrame({
#             "chrom": raw.iloc[:, 0],
#             "start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
#             "end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
#             "te_name": ".",
#             "te_type": raw.iloc[:, 3],
#         })
#     else:
#         raise ValueError("Background file must have at least 4 columns.")
#
#     bg = bg.dropna(subset=["chrom", "start", "end", "te_type"]).copy()
#     bg = bg[bg["end"] > bg["start"]].copy()
#     bg["te_bin"] = bg["te_type"].map(map_te_bin)
#     bg = bg[bg["te_bin"].notna()].copy()
#     bg["te_id"] = (
#         bg["chrom"].astype(str) + ":" +
#         bg["start"].astype(int).astype(str) + "-" +
#         bg["end"].astype(int).astype(str) + ":" +
#         bg["te_name"].astype(str) + ":" +
#         bg["te_type"].astype(str)
#     )
#     return bg
#
#
# def load_mouse_overlap(path):
#     raw = pd.read_csv(path, sep="\t", header=None)
#     if raw.shape[1] < 10:
#         raise ValueError(
#             "Overlap file has fewer columns than expected. "
#             "Expected enhancer columns + TE_chrom TE_start TE_end TE_name TE_type overlap_bp."
#         )
#
#     df = pd.DataFrame({
#         "enh_chrom": raw.iloc[:, 0],
#         "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
#         "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
#         "enh_id": raw.iloc[:, 3].astype(str),
#         "te_chrom": raw.iloc[:, -6],
#         "te_start": pd.to_numeric(raw.iloc[:, -5], errors="coerce"),
#         "te_end": pd.to_numeric(raw.iloc[:, -4], errors="coerce"),
#         "te_name": raw.iloc[:, -3].astype(str),
#         "te_type": raw.iloc[:, -2],
#         "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
#     })
#
#     df = df.dropna(subset=[
#         "enh_chrom", "enh_start", "enh_end", "enh_id",
#         "te_chrom", "te_start", "te_end", "te_type", "overlap_bp"
#     ]).copy()
#
#     df = df[df["te_type"] != "."].copy()
#     df = df[df["overlap_bp"] >= MIN_OVERLAP_BP].copy()
#     df["te_bin"] = df["te_type"].map(map_te_bin)
#     df = df[df["te_bin"].notna()].copy()
#     df["tissue"] = df["enh_id"].map(assign_tissue_from_enh_id)
#
#     df["te_id"] = (
#         df["te_chrom"].astype(str) + ":" +
#         df["te_start"].astype(int).astype(str) + "-" +
#         df["te_end"].astype(int).astype(str) + ":" +
#         df["te_name"].astype(str) + ":" +
#         df["te_type"].astype(str)
#     )
#     return df
#
#
# def unique_enhancers(df):
#     enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
#     enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
#     enh["win_end"] = enh["enh_end"] + FLANK_BP
#     return enh
#
#
# # ----------------------------
# # observed and reference counts
# # ----------------------------
# def observed_enhancer_te_pair_counts(df):
#     pairs = df.drop_duplicates(["enh_id", "te_id"]).copy()
#     counts = pairs["te_bin"].value_counts().reindex(cat_order, fill_value=0).astype(int)
#     return counts
#
#
# def observed_enhancer_te_pair_composition(df):
#     counts = observed_enhancer_te_pair_counts(df).astype(float)
#     return counts / counts.sum() * 100
#
#
# def genome_insert_counts(bg):
#     unique_tes = bg.drop_duplicates("te_id").copy()
#     return unique_tes["te_bin"].value_counts().reindex(cat_order, fill_value=0).astype(int)
#
#
# def genome_insert_composition(bg):
#     counts = genome_insert_counts(bg).astype(float)
#     return counts / counts.sum() * 100
#
#
# def neighborhood_insert_counts_per_enhancer(enhancers, bg):
#     rows = []
#     bg_unique = bg.drop_duplicates("te_id").copy()
#
#     for chrom, enh_chr in enhancers.groupby("enh_chrom"):
#         te_chr = bg_unique[bg_unique["chrom"] == chrom].copy()
#         te_st = te_chr["start"].to_numpy()
#         te_en = te_chr["end"].to_numpy()
#         te_bin = te_chr["te_bin"].to_numpy()
#
#         for _, row in enh_chr.iterrows():
#             ws = row["win_start"]
#             we = row["win_end"]
#             mask = (te_en > ws) & (te_st < we)
#
#             if not np.any(mask):
#                 counts = pd.Series(0, index=cat_order, name=row["enh_id"])
#                 rows.append(counts)
#                 continue
#
#             sub_bin = te_bin[mask]
#             counts = pd.Series(sub_bin).value_counts().reindex(cat_order, fill_value=0).astype(int)
#             counts.name = row["enh_id"]
#             rows.append(counts)
#
#     out = pd.DataFrame(rows)
#     out.index = enhancers["enh_id"].tolist()
#     out = out.reindex(columns=cat_order, fill_value=0)
#     return out.astype(int)
#
#
# def neighborhood_total_insert_counts(enhancers, bg):
#     per_enh = neighborhood_insert_counts_per_enhancer(enhancers, bg)
#     return per_enh.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)
#
#
# def neighborhood_total_insert_composition(enhancers, bg):
#     counts = neighborhood_total_insert_counts(enhancers, bg).astype(float)
#     return counts / counts.sum() * 100
#
#
# # ----------------------------
# # statistics
# # ----------------------------
# def run_gtest_and_chisq(obs_counts, ref_counts, label):
#     obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
#     ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)
#
#     obs_total = obs_counts.sum()
#     ref_total = ref_counts.sum()
#     ref_props = ref_counts / ref_total
#     expected = ref_props * obs_total
#     expected = expected * (obs_total / expected.sum())
#
#     keep = ~((obs_counts == 0) & (expected == 0))
#     obs_use = obs_counts[keep]
#     exp_use = expected[keep]
#
#     if np.any((exp_use == 0) & (obs_use > 0)):
#         summary = pd.DataFrame([{
#             "comparison": label,
#             "obs_total": int(obs_total),
#             "ref_total": int(ref_total),
#             "n_categories_tested": int(keep.sum()),
#             "chi2_stat": np.nan,
#             "chi2_p": np.nan,
#             "chi2_log10_p": np.nan,
#             "g_stat": np.nan,
#             "g_p": np.nan,
#             "g_log10_p": np.nan,
#             "note": "Invalid because at least one category has expected=0 but observed>0"
#         }])
#
#         diagnostic = pd.DataFrame({
#             "comparison": label,
#             "te_class": cat_order,
#             "observed": obs_counts.values,
#             "expected": expected.values,
#             "obs_prop": obs_counts.values / obs_total,
#             "ref_prop": ref_props.values,
#             "chi_contribution": np.nan
#         })
#         return summary, diagnostic
#
#     chi_stat, chi_p = chisquare(f_obs=obs_use.values, f_exp=exp_use.values)
#     g_stat, g_p = power_divergence(f_obs=obs_use.values, f_exp=exp_use.values, lambda_="log-likelihood")
#
#     df = len(obs_use) - 1
#     chi_log10_p = chi2_dist.logsf(chi_stat, df) / np.log(10)
#     g_log10_p = chi2_dist.logsf(g_stat, df) / np.log(10)
#
#     diagnostic = pd.DataFrame({
#         "comparison": label,
#         "te_class": obs_use.index,
#         "observed": obs_use.values,
#         "expected": exp_use.values,
#         "obs_prop": obs_use.values / obs_total,
#         "ref_prop": ref_props[obs_use.index].values,
#         "chi_contribution": ((obs_use.values - exp_use.values) ** 2) / exp_use.values
#     })
#
#     summary = pd.DataFrame([{
#         "comparison": label,
#         "obs_total": int(obs_total),
#         "ref_total": int(ref_total),
#         "n_categories_tested": int(keep.sum()),
#         "chi2_stat": chi_stat,
#         "chi2_p": chi_p,
#         "chi2_log10_p": chi_log10_p,
#         "g_stat": g_stat,
#         "g_p": g_p,
#         "g_log10_p": g_log10_p,
#         "note": ""
#     }])
#     return summary, diagnostic
#
#
# def per_class_insert_enrichment(obs_counts, ref_counts, group_name):
#     obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
#     ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)
#
#     obs_total = obs_counts.sum()
#     ref_total = ref_counts.sum()
#     rows = []
#
#     for cat in cat_order:
#         a = obs_counts[cat]
#         b = obs_total - a
#         c = ref_counts[cat]
#         d = ref_total - c
#
#         obs_prop = a / obs_total
#         ref_prop = c / ref_total
#         log2_enrichment = np.log2((obs_prop + EPS) / (ref_prop + EPS))
#         log2_or = np.log2(((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5)))
#
#         chi2, p, _, _ = chi2_contingency([[a, b], [c, d]], correction=False)
#
#         rows.append({
#             "group": group_name,
#             "te_class": cat,
#             "obs_count": int(a),
#             "obs_total": int(obs_total),
#             "ref_count": int(c),
#             "ref_total": int(ref_total),
#             "obs_prop": obs_prop,
#             "ref_prop": ref_prop,
#             "log2_obs_over_ref": log2_enrichment,
#             "log2_odds_ratio": log2_or,
#             "chi2": chi2,
#             "p_value": p,
#             "direction": "enriched" if obs_prop > ref_prop else "depleted"
#         })
#
#     out = pd.DataFrame(rows)
#     out["q_value"] = bh_fdr(out["p_value"].values)
#     out["label"] = out["q_value"].map(p_to_star)
#     return out
#
#
# # ----------------------------
# # pooled plots
# # ----------------------------
# def run_pooled_analysis(shadow_df, single_df, genome_df):
#     shadow_enh = unique_enhancers(shadow_df)
#     single_enh = unique_enhancers(single_df)
#
#     shadow_obs_counts = observed_enhancer_te_pair_counts(shadow_df)
#     single_obs_counts = observed_enhancer_te_pair_counts(single_df)
#     shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
#     single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)
#     genome_counts = genome_insert_counts(genome_df)
#
#     shadow_obs = shadow_obs_counts / shadow_obs_counts.sum() * 100
#     single_obs = single_obs_counts / single_obs_counts.sum() * 100
#     shadow_local = shadow_local_counts / shadow_local_counts.sum() * 100
#     single_local = single_local_counts / single_local_counts.sum() * 100
#     genome_comp = genome_counts / genome_counts.sum() * 100
#
#     summary = pd.DataFrame({
#         "shadow_observed_enhancer_TE_pairs": shadow_obs,
#         "shadow_neighborhood_enhancer_window_TE_pairs": shadow_local,
#         "single_observed_enhancer_TE_pairs": single_obs,
#         "single_neighborhood_enhancer_window_TE_pairs": single_local,
#         "genome_TE_inserts": genome_comp
#     }).reindex(cat_order)
#     summary.to_csv(f"{output_dir}/pooled_optionB_enhancer_TE_pair_summary.tsv", sep="\t")
#
#     stats_list = []
#     diag_list = []
#     comparisons = [
#         (shadow_obs_counts, shadow_local_counts, "Pooled shadow observed enhancer-TE pairs vs shadow neighborhood TE opportunities"),
#         (shadow_obs_counts, genome_counts, "Pooled shadow observed enhancer-TE pairs vs genome TE inserts"),
#         (single_obs_counts, single_local_counts, "Pooled single observed enhancer-TE pairs vs single neighborhood TE opportunities"),
#         (single_obs_counts, genome_counts, "Pooled single observed enhancer-TE pairs vs genome TE inserts"),
#         (shadow_obs_counts, single_obs_counts, "Pooled shadow observed enhancer-TE pairs vs single observed enhancer-TE pairs")
#     ]
#
#     for obs, ref, label in comparisons:
#         s, d = run_gtest_and_chisq(obs, ref, label)
#         stats_list.append(s)
#         diag_list.append(d)
#
#     pd.concat(stats_list, ignore_index=True).to_csv(
#         f"{output_dir}/pooled_optionB_gtest_chisq_stats.tsv", sep="\t", index=False
#     )
#     pd.concat(diag_list, ignore_index=True).to_csv(
#         f"{output_dir}/pooled_optionB_gof_diagnostics_by_TEclass.tsv", sep="\t", index=False
#     )
#
#     shadow_estats = per_class_insert_enrichment(shadow_obs_counts, shadow_local_counts, "Shadow")
#     single_estats = per_class_insert_enrichment(single_obs_counts, single_local_counts, "Single")
#     enrichment_stats = pd.concat([shadow_estats, single_estats], ignore_index=True)
#     enrichment_stats.to_csv(
#         f"{output_dir}/pooled_optionB_TEclass_enhancer_TE_pairs_vs_neighborhood_BH.tsv",
#         sep="\t",
#         index=False
#     )
#
#     plot_pooled_stacked(summary)
#     plot_pooled_enrichment_bars(enrichment_stats)
#
#
# def plot_pooled_stacked(summary):
#     fig, ax = plt.subplots(figsize=STACKED_FIGSIZE)
#
#     all_values = [
#         summary["shadow_observed_enhancer_TE_pairs"].values,
#         summary["shadow_neighborhood_enhancer_window_TE_pairs"].values,
#         summary["single_observed_enhancer_TE_pairs"].values,
#         summary["single_neighborhood_enhancer_window_TE_pairs"].values,
#         summary["genome_TE_inserts"].values
#     ]
#
#     all_labels = [
#         "shadows",
#         "shadow\nneighborhood",
#         "single",
#         "single\nneighborhood",
#         "genome"
#     ]
#
#     x_positions = [0.0, 1.0, 2.6, 3.6, 5.2]
#
#     for xpos, vals in zip(x_positions, all_values):
#         plot_stacked_bar(ax, vals, xpos, width=0.62)
#
#     ax.set_ylabel("TE insert %")
#     ax.set_xticks(x_positions)
#     ax.set_xticklabels(all_labels, fontsize=12)
#     ax.set_ylim(0, 100)
#     make_grouped_te_legend(ax)
#     style_axes(ax)
#
#     save_figure(
#         fig,
#         png_path=f"{output_dir}/pooled_optionB_stacked_TE_observed_local_genome_insert_based.png",
#         pdf_path=f"{output_dir}/pooled_optionB_stacked_TE_observed_local_genome_insert_based.pdf",
#         left=0.10,
#         right=0.80,
#         bottom=0.24,
#         top=0.92
#     )
#
#
# def plot_pooled_enrichment_bars(enrichment_stats):
#     labels = cat_order
#     shadow_vals = enrichment_stats[enrichment_stats["group"] == "Shadow"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()
#     single_vals = enrichment_stats[enrichment_stats["group"] == "Single"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()
#
#     x = np.arange(len(cat_order))
#     bar_width = 0.46
#
#     fig, ax = plt.subplots(figsize=ENRICHMENT_BAR_FIGSIZE)
#
#     bars1 = ax.bar(
#         x - bar_width / 2,
#         shadow_vals,
#         width=bar_width,
#         color=[te_color_map[lbl] for lbl in labels],
#         edgecolor="black",
#         hatch="//"
#     )
#
#     bars2 = ax.bar(
#         x + bar_width / 2,
#         single_vals,
#         width=bar_width,
#         color=[te_color_map[lbl] for lbl in labels],
#         edgecolor="black"
#     )
#
#     ax.set_xticks(x)
#     ax.set_xticklabels(labels, fontsize=12)
#     if PLOT_METRIC == "log2_obs_over_ref":
#         ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=12)
#     else:
#         ax.set_ylabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=12)
#     ax.tick_params(axis="y", labelsize=11)
#     ax.tick_params(axis="x", labelsize=11)
#     apply_enrichment_axis(ax)
#
#     for bars in [bars1, bars2]:
#         for bar in bars:
#             height = bar.get_height()
#             if np.isfinite(height):
#                 offset = 0.10 if height >= 0 else -0.20
#                 ax.text(
#                     bar.get_x() + bar.get_width() / 2,
#                     height + offset,
#                     f"{height:.2f}",
#                     ha="center",
#                     va="bottom" if height >= 0 else "top",
#                     fontsize=11
#                 )
#
#     for i, cat in enumerate(cat_order):
#         sh_row = enrichment_stats[(enrichment_stats["group"] == "Shadow") & (enrichment_stats["te_class"] == cat)].iloc[0]
#         si_row = enrichment_stats[(enrichment_stats["group"] == "Single") & (enrichment_stats["te_class"] == cat)].iloc[0]
#
#         if sh_row["label"] != "ns":
#             ax.text(
#                 x[i] - bar_width / 2,
#                 shadow_vals[i] + 0.28,
#                 sh_row["label"],
#                 ha="center",
#                 va="bottom",
#                 fontsize=11,
#                 fontweight="bold"
#             )
#
#         if si_row["label"] != "ns":
#             ax.text(
#                 x[i] + bar_width / 2,
#                 single_vals[i] + 0.28,
#                 si_row["label"],
#                 ha="center",
#                 va="bottom",
#                 fontsize=11,
#                 fontweight="bold"
#             )
#
#     shadow_patch = mpatches.Patch(facecolor="white", edgecolor="black", hatch="//", label="shadows")
#     single_patch = mpatches.Patch(facecolor="white", edgecolor="black", label="singles")
#     ax.legend(handles=[shadow_patch, single_patch], fontsize=11, loc="upper right", frameon=False)
#     style_axes(ax)
#
#     save_figure(
#         fig,
#         png_path=f"{output_dir}/pooled_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.png",
#         pdf_path=f"{output_dir}/pooled_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.pdf",
#         left=0.16,
#         right=0.98,
#         bottom=0.24,
#         top=0.92
#     )
#
#
# # ----------------------------
# # tissue-specific analysis
# # ----------------------------
# def run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df):
#     shadow_df = shadow_df_all[shadow_df_all["tissue"] == tissue].copy()
#     single_df = single_df_all[single_df_all["tissue"] == tissue].copy()
#
#     if shadow_df.empty or single_df.empty:
#         print(f"Skipping {tissue}: missing data in one group")
#         return None
#
#     shadow_enh = unique_enhancers(shadow_df)
#     single_enh = unique_enhancers(single_df)
#
#     shadow_obs_counts = observed_enhancer_te_pair_counts(shadow_df)
#     single_obs_counts = observed_enhancer_te_pair_counts(single_df)
#     shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
#     single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)
#     genome_counts = genome_insert_counts(genome_df)
#
#     shadow_obs = shadow_obs_counts / shadow_obs_counts.sum() * 100
#     single_obs = single_obs_counts / single_obs_counts.sum() * 100
#     shadow_local = shadow_local_counts / shadow_local_counts.sum() * 100
#     single_local = single_local_counts / single_local_counts.sum() * 100
#     genome_comp = genome_counts / genome_counts.sum() * 100
#
#     summary = pd.DataFrame({
#         "shadow_observed_enhancer_TE_pairs": shadow_obs,
#         "shadow_neighborhood_enhancer_window_TE_pairs": shadow_local,
#         "single_observed_enhancer_TE_pairs": single_obs,
#         "single_neighborhood_enhancer_window_TE_pairs": single_local,
#         "genome_TE_inserts": genome_comp
#     }).reindex(cat_order)
#
#     summary.to_csv(f"{output_dir}/{tissue}_optionB_enhancer_TE_pair_summary.tsv", sep="\t")
#
#     stats_list = []
#     diag_list = []
#     comparisons = [
#         (shadow_obs_counts, shadow_local_counts, f"{tissue} shadow observed enhancer-TE pairs vs shadow neighborhood TE opportunities"),
#         (shadow_obs_counts, genome_counts, f"{tissue} shadow observed enhancer-TE pairs vs genome TE inserts"),
#         (single_obs_counts, single_local_counts, f"{tissue} single observed enhancer-TE pairs vs single neighborhood TE opportunities"),
#         (single_obs_counts, genome_counts, f"{tissue} single observed enhancer-TE pairs vs genome TE inserts"),
#     ]
#
#     for obs, ref, label in comparisons:
#         s, d = run_gtest_and_chisq(obs, ref, label)
#         stats_list.append(s)
#         diag_list.append(d)
#
#     pd.concat(stats_list, ignore_index=True).to_csv(
#         f"{output_dir}/{tissue}_optionB_gtest_chisq_stats.tsv", sep="\t", index=False
#     )
#     pd.concat(diag_list, ignore_index=True).to_csv(
#         f"{output_dir}/{tissue}_optionB_gof_diagnostics_by_TEclass.tsv", sep="\t", index=False
#     )
#
#     shadow_estats = per_class_insert_enrichment(shadow_obs_counts, shadow_local_counts, "Shadow")
#     single_estats = per_class_insert_enrichment(single_obs_counts, single_local_counts, "Single")
#     enrichment_stats = pd.concat([shadow_estats, single_estats], ignore_index=True)
#
#     enrichment_stats.to_csv(
#         f"{output_dir}/{tissue}_optionB_TEclass_enhancer_TE_pairs_vs_neighborhood_BH.tsv",
#         sep="\t",
#         index=False
#     )
#
#     plot_tissue_stacked(tissue, summary)
#     plot_tissue_enrichment_bars(tissue, enrichment_stats)
#     return enrichment_stats
#
#
# def plot_tissue_stacked(tissue, summary):
#     fig, ax = plt.subplots(figsize=STACKED_FIGSIZE)
#
#     all_values = [
#         summary["shadow_observed_enhancer_TE_pairs"].values,
#         summary["shadow_neighborhood_enhancer_window_TE_pairs"].values,
#         summary["single_observed_enhancer_TE_pairs"].values,
#         summary["single_neighborhood_enhancer_window_TE_pairs"].values,
#         summary["genome_TE_inserts"].values
#     ]
#     all_labels = [
#         "shadows",
#         "shadow\nneighborhood",
#         "single",
#         "single\nneighborhood",
#         "genome"
#     ]
#     x_positions = [0.0, 1.0, 2.6, 3.6, 5.2]
#
#     for xpos, vals in zip(x_positions, all_values):
#         plot_stacked_bar(ax, vals, xpos, width=0.62)
#
#     ax.set_title(tissue.capitalize(), fontsize=13)
#     ax.set_ylabel("TE insert %")
#     ax.set_xticks(x_positions)
#     ax.set_xticklabels(all_labels, fontsize=12)
#     ax.set_ylim(0, 100)
#     make_grouped_te_legend(ax)
#     style_axes(ax)
#
#     save_figure(
#         fig,
#         png_path=f"{output_dir}/{tissue}_optionB_stacked_TE_observed_local_genome_insert_based.png",
#         pdf_path=f"{output_dir}/{tissue}_optionB_stacked_TE_observed_local_genome_insert_based.pdf",
#         left=0.10,
#         right=0.80,
#         bottom=0.24,
#         top=0.90
#     )
#
#
# def plot_tissue_enrichment_bars(tissue, enrichment_stats):
#     shadow_vals = enrichment_stats[enrichment_stats["group"] == "Shadow"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()
#     single_vals = enrichment_stats[enrichment_stats["group"] == "Single"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()
#
#     x = np.arange(len(cat_order))
#     bar_width = 0.46
#     fig, ax = plt.subplots(figsize=ENRICHMENT_BAR_FIGSIZE)
#
#     bars1 = ax.bar(
#         x - bar_width / 2,
#         shadow_vals,
#         width=bar_width,
#         color=[te_color_map[lbl] for lbl in cat_order],
#         edgecolor="black",
#         hatch="//"
#     )
#     bars2 = ax.bar(
#         x + bar_width / 2,
#         single_vals,
#         width=bar_width,
#         color=[te_color_map[lbl] for lbl in cat_order],
#         edgecolor="black"
#     )
#
#     ax.set_title(tissue.capitalize(), fontsize=13)
#     ax.set_xticks(x)
#     ax.set_xticklabels(cat_order, fontsize=12)
#     if PLOT_METRIC == "log2_obs_over_ref":
#         ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=11)
#     else:
#         ax.set_ylabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=12)
#     apply_enrichment_axis(ax)
#
#     for bars in [bars1, bars2]:
#         for bar in bars:
#             height = bar.get_height()
#             if np.isfinite(height):
#                 offset = 0.10 if height >= 0 else -0.20
#                 ax.text(
#                     bar.get_x() + bar.get_width() / 2,
#                     height + offset,
#                     f"{height:.2f}",
#                     ha="center",
#                     va="bottom" if height >= 0 else "top",
#                     fontsize=11
#                 )
#
#     for i, cat in enumerate(cat_order):
#         sh_row = enrichment_stats[(enrichment_stats["group"] == "Shadow") & (enrichment_stats["te_class"] == cat)].iloc[0]
#         si_row = enrichment_stats[(enrichment_stats["group"] == "Single") & (enrichment_stats["te_class"] == cat)].iloc[0]
#
#         if sh_row["label"] != "ns":
#             ax.text(x[i] - bar_width / 2, shadow_vals[i] + 0.28, sh_row["label"], ha="center", va="bottom", fontsize=11, fontweight="bold")
#         if si_row["label"] != "ns":
#             ax.text(x[i] + bar_width / 2, single_vals[i] + 0.28, si_row["label"], ha="center", va="bottom", fontsize=11, fontweight="bold")
#
#     shadow_patch = mpatches.Patch(facecolor="white", edgecolor="black", hatch="//", label="shadows")
#     single_patch = mpatches.Patch(facecolor="white", edgecolor="black", label="singles")
#     ax.legend(handles=[shadow_patch, single_patch], fontsize=11, loc="upper right", frameon=False)
#     style_axes(ax)
#
#     save_figure(
#         fig,
#         png_path=f"{output_dir}/{tissue}_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.png",
#         pdf_path=f"{output_dir}/{tissue}_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.pdf",
#         left=0.16,
#         right=0.98,
#         bottom=0.24,
#         top=0.90
#     )
#
#
# # ----------------------------
# # tissue enrichment summary tables
# # ----------------------------
# def compute_tissue_enrichment_tables(shadow_df_all, single_df_all, genome_df):
#     rows = []
#
#     for tissue in tissue_order:
#         shadow_df = shadow_df_all[shadow_df_all["tissue"] == tissue].copy()
#         single_df = single_df_all[single_df_all["tissue"] == tissue].copy()
#
#         if shadow_df.empty or single_df.empty:
#             print(f"Skipping {tissue}: no data in one of the groups")
#             continue
#
#         shadow_enh = unique_enhancers(shadow_df)
#         single_enh = unique_enhancers(single_df)
#
#         shadow_obs_counts = observed_enhancer_te_pair_counts(shadow_df)
#         single_obs_counts = observed_enhancer_te_pair_counts(single_df)
#         shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
#         single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)
#
#         shadow_stats = per_class_insert_enrichment(shadow_obs_counts, shadow_local_counts, group_name="shadows").copy()
#         shadow_stats["tissue"] = tissue
#
#         single_stats = per_class_insert_enrichment(single_obs_counts, single_local_counts, group_name="singles").copy()
#         single_stats["tissue"] = tissue
#
#         rows.append(shadow_stats)
#         rows.append(single_stats)
#
#     shadow_enh_all = unique_enhancers(shadow_df_all)
#     single_enh_all = unique_enhancers(single_df_all)
#     shadow_obs_counts_all = observed_enhancer_te_pair_counts(shadow_df_all)
#     single_obs_counts_all = observed_enhancer_te_pair_counts(single_df_all)
#     shadow_local_counts_all = neighborhood_total_insert_counts(shadow_enh_all, genome_df)
#     single_local_counts_all = neighborhood_total_insert_counts(single_enh_all, genome_df)
#
#     shadow_stats_all = per_class_insert_enrichment(shadow_obs_counts_all, shadow_local_counts_all, group_name="shadows").copy()
#     shadow_stats_all["tissue"] = "pooled"
#
#     single_stats_all = per_class_insert_enrichment(single_obs_counts_all, single_local_counts_all, group_name="singles").copy()
#     single_stats_all["tissue"] = "pooled"
#
#     rows.append(shadow_stats_all)
#     rows.append(single_stats_all)
#
#     out = pd.concat(rows, ignore_index=True)
#     out["te_class"] = pd.Categorical(out["te_class"], categories=cat_order, ordered=True)
#     out["tissue"] = pd.Categorical(out["tissue"], categories=tissue_order + ["pooled"], ordered=True)
#     return out.sort_values(["group", "te_class", "tissue"])
#
#
# # ----------------------------
# # tissue split dot plot
# # ----------------------------
# def plot_tissue_dotplot(enrichment_df, metric="log2_obs_over_ref", outfile_png=None, outfile_pdf=None):
#     fig, axes = plt.subplots(1, 2, figsize=TISSUE_SPLIT_FIGSIZE, sharey=True)
#
#     group_order = ["shadows", "singles"]
#     group_titles = {"shadows": "Shadows", "singles": "Singles"}
#     yvals = np.arange(len(cat_order))
#     tissue_yoffset = {"forebrain": -0.18, "heart": 0.00, "limb": 0.18, "pooled": 0.00}
#
#     for ax, group in zip(axes, group_order):
#         sub = enrichment_df[enrichment_df["group"] == group].copy()
#         ax.axvline(0, color="black", linewidth=1)
#
#         for tissue in tissue_order:
#             tissue_sub = sub[sub["tissue"] == tissue].copy()
#             for i, te in enumerate(cat_order):
#                 row = tissue_sub[tissue_sub["te_class"] == te]
#                 if row.empty:
#                     continue
#
#                 x = float(row.iloc[0][metric])
#                 y = i + tissue_yoffset[tissue]
#                 ax.scatter(
#                     x,
#                     y,
#                     s=55,
#                     color=tissue_colors[tissue],
#                     edgecolor="black",
#                     linewidth=0.8,
#                     zorder=3
#                 )
#
#                 label = str(row.iloc[0]["label"])
#                 if label != "ns":
#                     ax.text(
#                         x + 0.08,
#                         y,
#                         label,
#                         ha="left",
#                         va="center",
#                         fontsize=9,
#                         fontweight="bold"
#                     )
#
#         pooled_sub = sub[sub["tissue"] == "pooled"].copy()
#         for i, te in enumerate(cat_order):
#             row = pooled_sub[pooled_sub["te_class"] == te]
#             if row.empty:
#                 continue
#
#             x = float(row.iloc[0][metric])
#             y = i + tissue_yoffset["pooled"]
#             ax.scatter(
#                 x,
#                 y,
#                 s=75,
#                 color=tissue_colors["pooled"],
#                 edgecolor="black",
#                 linewidth=0.8,
#                 marker="D",
#                 zorder=4
#             )
#
#         ax.set_title(group_titles[group], fontsize=14)
#         ax.set_xlim(*ENRICHMENT_YLIM)
#         ax.set_yticks(yvals)
#         ax.set_yticklabels(cat_order, fontsize=12)
#         ax.invert_yaxis()
#         style_axes(ax)
#
#     if metric == "log2_obs_over_ref":
#         axes[0].set_ylabel("TE class", fontsize=12)
#         fig.supxlabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=13, y=0.08)
#     else:
#         axes[0].set_ylabel("TE class", fontsize=12)
#         fig.supxlabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=13, y=0.08)
#
#     legend_handles = [
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["forebrain"], markeredgecolor="black", markersize=7, label="Forebrain"),
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["heart"], markeredgecolor="black", markersize=7, label="Heart"),
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["limb"], markeredgecolor="black", markersize=7, label="Limb"),
#         Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Pooled")
#     ]
#
#     fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=4, frameon=False, fontsize=11)
#
#     save_figure(fig, png_path=outfile_png, pdf_path=outfile_pdf, left=0.10, right=0.98, bottom=0.24, top=0.82)
#
#
# # ----------------------------
# # final combined vertical dot plot
# # ----------------------------
# def plot_combined_vertical_dotplot_clean(
#     enrichment_df,
#     metric="log2_obs_over_ref",
#     outfile_png=None,
#     outfile_pdf=None,
#     show_pooled=True,
#     show_stars=True
# ):
#     fig, ax = plt.subplots(figsize=COMBINED_DOTPLOT_FIGSIZE)
#
#     # Make category spacing larger so points do not bleed into neighboring TE classes.
#     x_base = np.arange(len(cat_order)) * 1.45
#     tissue_xoffset = {
#         "forebrain": -0.30,
#         "heart": 0.00,
#         "limb": 0.30,
#         "pooled": 0.55
#     }
#     group_xoffset = {
#         "shadows": -0.05,
#         "singles": 0.05
#     }
#     star_offset = 0.16
#
#     ax.axhline(0, color="black", linewidth=1.1, zorder=1)
#     for x_sep in (x_base[:-1] + x_base[1:]) / 2:
#         ax.axvline(x_sep, color="0.90", linewidth=0.8, zorder=0)
#
#     for tissue in tissue_order:
#         color = tissue_colors[tissue]
#
#         for te_i, te in enumerate(cat_order):
#             sh = enrichment_df[
#                 (enrichment_df["tissue"] == tissue) &
#                 (enrichment_df["group"] == "shadows") &
#                 (enrichment_df["te_class"] == te)
#             ]
#             si = enrichment_df[
#                 (enrichment_df["tissue"] == tissue) &
#                 (enrichment_df["group"] == "singles") &
#                 (enrichment_df["te_class"] == te)
#             ]
#
#             if sh.empty or si.empty:
#                 continue
#
#             x_sh = x_base[te_i] + tissue_xoffset[tissue] + group_xoffset["shadows"]
#             x_si = x_base[te_i] + tissue_xoffset[tissue] + group_xoffset["singles"]
#             y_sh = float(sh.iloc[0][metric])
#             y_si = float(si.iloc[0][metric])
#
#             # No connector line between shadow and single points.
#             ax.scatter(
#                 x_sh,
#                 y_sh,
#                 s=58,
#                 facecolor=color,
#                 edgecolor="black",
#                 linewidth=0.9,
#                 zorder=4
#             )
#             ax.scatter(
#                 x_si,
#                 y_si,
#                 s=58,
#                 facecolor="white",
#                 edgecolor=color,
#                 linewidth=1.8,
#                 zorder=4
#             )
#
#             if show_stars:
#                 for x_pt, y_pt, row in [
#                     (x_sh, y_sh, sh.iloc[0]),
#                     (x_si, y_si, si.iloc[0])
#                 ]:
#                     label = str(row["label"])
#                     if label != "ns":
#                         ax.text(
#                             x_pt,
#                             y_pt + star_offset,
#                             label,
#                             ha="center",
#                             va="bottom",
#                             fontsize=8,
#                             fontweight="bold",
#                             zorder=6,
#                             clip_on=False
#                         )
#
#     if show_pooled:
#         for te_i, te in enumerate(cat_order):
#             for group in ["shadows", "singles"]:
#                 row = enrichment_df[
#                     (enrichment_df["tissue"] == "pooled") &
#                     (enrichment_df["group"] == group) &
#                     (enrichment_df["te_class"] == te)
#                 ]
#                 if row.empty:
#                     continue
#
#                 x = x_base[te_i] + tissue_xoffset["pooled"] + group_xoffset[group]
#                 y = float(row.iloc[0][metric])
#
#                 if group == "shadows":
#                     ax.scatter(x, y, s=84, marker="D", facecolor="black", edgecolor="black", linewidth=1.0, zorder=5)
#                 else:
#                     ax.scatter(x, y, s=84, marker="D", facecolor="white", edgecolor="black", linewidth=1.8, zorder=5)
#
#                 if show_stars:
#                     label = str(row.iloc[0]["label"])
#                     if label != "ns":
#                         ax.text(
#                             x,
#                             y + star_offset,
#                             label,
#                             ha="center",
#                             va="bottom",
#                             fontsize=8,
#                             fontweight="bold",
#                             zorder=6,
#                             clip_on=False
#                         )
#
#     ax.set_xticks(x_base)
#     ax.set_xticklabels(cat_order, fontsize=13)
#     ax.set_xlim(x_base[0] - 0.65, x_base[-1] + 0.90)
#
#     if metric == "log2_obs_over_ref":
#         ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=13)
#     else:
#         ax.set_ylabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=13)
#
#     apply_enrichment_axis(ax)
#     style_axes(ax)
#
#     tissue_handles = [
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["forebrain"], markeredgecolor="black", markersize=7, label="Forebrain"),
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["heart"], markeredgecolor="black", markersize=7, label="Heart"),
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["limb"], markeredgecolor="black", markersize=7, label="Limb")
#     ]
#     if show_pooled:
#         tissue_handles.append(
#             Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Pooled")
#         )
#
#     class_handles = [
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Shadows"),
#         Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="white", markeredgecolor="black", markersize=7, label="Singles")
#     ]
#
#     leg1 = ax.legend(
#         handles=tissue_handles,
#         title="Tissue",
#         frameon=False,
#         fontsize=10,
#         title_fontsize=11,
#         loc="upper left",
#         bbox_to_anchor=(1.02, 1.00)
#     )
#     ax.add_artist(leg1)
#
#     ax.legend(
#         handles=class_handles,
#         title="Enhancer class",
#         frameon=False,
#         fontsize=10,
#         title_fontsize=11,
#         loc="upper left",
#         bbox_to_anchor=(1.02, 0.50)
#     )
#
#     save_figure(fig, png_path=outfile_png, pdf_path=outfile_pdf, left=0.10, right=0.77, bottom=0.24, top=0.92)
#
#
# # ----------------------------
# # run
# # ----------------------------
# def main():
#     genome_df = load_background(genome_file)
#     shadow_df_all = load_mouse_overlap(shadow_file)
#     single_df_all = load_mouse_overlap(single_file)
#
#     shadow_df_all = shadow_df_all[shadow_df_all["tissue"].isin(tissue_order)].copy()
#     single_df_all = single_df_all[single_df_all["tissue"].isin(tissue_order)].copy()
#
#     print("Shadow tissues:")
#     print(shadow_df_all["tissue"].value_counts())
#     print("\nSingle tissues:")
#     print(single_df_all["tissue"].value_counts())
#
#     run_pooled_analysis(shadow_df_all, single_df_all, genome_df)
#
#     for tissue in tissue_order:
#         run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df)
#
#     enrichment_df = compute_tissue_enrichment_tables(
#         shadow_df_all=shadow_df_all,
#         single_df_all=single_df_all,
#         genome_df=genome_df
#     )
#
#     print("\nCombined tissue enrichment table:")
#     print(enrichment_df.round(4))
#
#     enrichment_df.to_csv(
#         f"{output_dir}/tissue_dotplot_enrichment_table.tsv",
#         sep="\t",
#         index=False
#     )
#
#
#
#     plot_combined_vertical_dotplot_clean(
#         enrichment_df,
#         metric=PLOT_METRIC,
#         outfile_png=f"{output_dir}/combined_vertical_dotplot_TEclass_enrichment_clean.png",
#         outfile_pdf=f"{output_dir}/combined_vertical_dotplot_TEclass_enrichment_clean.pdf",
#         show_pooled=True,
#         show_stars=True
#     )
#
#
# if __name__ == "__main__":
#     main()
##break
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, VPacker, HPacker, TextArea, DrawingArea
import matplotlib.patches as patches
from scipy.stats import chisquare, power_divergence, chi2_contingency, chi2 as chi2_dist

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "input"
OUTPUT_FILES = ROOT / "outputs" / "files"
OUTPUT_FIGURES = ROOT / "outputs" / "figures"
OUTPUT_FILES.mkdir(parents=True, exist_ok=True)
OUTPUT_FIGURES.mkdir(parents=True, exist_ok=True)

# ----------------------------
# file inputs
# ----------------------------
genome_file = str(INPUT_DIR / "FINAL_TE_genomefile_merged_dedup.bed")
shadow_file = str(INPUT_DIR / "Shadows_mouse_TEcooption_lastcol_Final.bed")
single_file = str(INPUT_DIR / "singlesmouse_TEcooption.bed")
output_dir = str(OUTPUT_FIGURES)


# ----------------------------
# parameters
# ----------------------------
MIN_OVERLAP_BP = 50
FLANK_BP = 30000
EPS = 1e-9
PLOT_METRIC = "log2_obs_over_ref"   # or "log2_odds_ratio"

cat_order = ["LTR", "LINE", "SINE", "TIR"]
tissue_order = ["forebrain", "heart", "limb"]
te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#f16913']
te_color_map = {
    "LTR": "#c6dbef",
    "LINE": "#6baed6",
    "SINE": "#1f78b4",
    "TIR": "#f16913"
}
tissue_colors = {
    "forebrain": "#4C78A8",
    "heart": "#E45756",
    "limb": "#54A24B",
    "pooled": "black"
}

# Shared physical height so stacked plots and tissue enrichment plots match.
PANEL_HEIGHT = 3.0
STACKED_FIGSIZE = (8.0, PANEL_HEIGHT)
ENRICHMENT_BAR_FIGSIZE = (5.0, PANEL_HEIGHT)
TISSUE_SPLIT_FIGSIZE = (8.3, PANEL_HEIGHT)
COMBINED_DOTPLOT_FIGSIZE = (7.8, PANEL_HEIGHT)

# Shared enrichment y-scale so all enrichment-style plots use the same vertical range.
ENRICHMENT_YLIM = (-3, 3.2)

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial"
})


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
    if te_type.startswith("DNA"):
        return "TIR"
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
        labelcolor="black"
    )


def save_figure(fig, png_path=None, pdf_path=None, left=0.10, right=0.82, bottom=0.22, top=0.92):
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    if png_path:
        fig.savefig(png_path, dpi=600)
    if pdf_path:
        fig.savefig(pdf_path)
    plt.close(fig)


def make_grouped_te_legend(ax):
    box_w = 34
    box_text_sep = 6

    def legend_row(color, label):
        da = DrawingArea(box_w, 14, 0, 0)
        rect = patches.Rectangle(
            (0, 2), 30, 10,
            facecolor=color,
            edgecolor="black",
            linewidth=1.2
        )
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
        bbox_to_anchor=(1.02, 1.0),
        bbox_transform=ax.transAxes,
        borderpad=0
    )
    ax.add_artist(anchored_box)


def plot_stacked_bar(ax, values, xpos, width=0.62):
    bottom = 0
    for val, color in zip(values, te_colors):
        ax.bar(xpos, val, width=width, bottom=bottom, color=color, edgecolor="black")
        if val > 0:
            if val < 6:
                ax.text(xpos + 0.33, bottom + val / 2, f"{val:.1f}%", ha="left", va="center", fontsize=11)
            else:
                ax.text(xpos, bottom + val / 2, f"{val:.1f}%", ha="center", va="center", fontsize=11)
        bottom += val


def apply_enrichment_axis(ax):
    ax.set_ylim(*ENRICHMENT_YLIM)
    ax.axhline(0, color="black", linewidth=1)

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
# ----------------------------
# loading
# ----------------------------
def load_background(path):
    raw = pd.read_csv(path, sep="\t", header=None)

    if raw.shape[1] >= 5:
        bg = pd.DataFrame({
            "chrom": raw.iloc[:, 0],
            "start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
            "end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
            "te_name": raw.iloc[:, 3].astype(str),
            "te_type": raw.iloc[:, 4],
        })
    elif raw.shape[1] == 4:
        bg = pd.DataFrame({
            "chrom": raw.iloc[:, 0],
            "start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
            "end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
            "te_name": ".",
            "te_type": raw.iloc[:, 3],
        })
    else:
        raise ValueError("Background file must have at least 4 columns.")

    bg = bg.dropna(subset=["chrom", "start", "end", "te_type"]).copy()
    bg = bg[bg["end"] > bg["start"]].copy()
    bg["te_bin"] = bg["te_type"].map(map_te_bin)
    bg = bg[bg["te_bin"].notna()].copy()
    bg["te_id"] = (
        bg["chrom"].astype(str) + ":" +
        bg["start"].astype(int).astype(str) + "-" +
        bg["end"].astype(int).astype(str) + ":" +
        bg["te_name"].astype(str) + ":" +
        bg["te_type"].astype(str)
    )
    return bg


def load_mouse_overlap(path):
    raw = pd.read_csv(path, sep="\t", header=None)
    if raw.shape[1] < 10:
        raise ValueError(
            "Overlap file has fewer columns than expected. "
            "Expected enhancer columns + TE_chrom TE_start TE_end TE_name TE_type overlap_bp."
        )

    df = pd.DataFrame({
        "enh_chrom": raw.iloc[:, 0],
        "enh_start": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        "enh_end": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
        "enh_id": raw.iloc[:, 3].astype(str),
        "te_chrom": raw.iloc[:, -6],
        "te_start": pd.to_numeric(raw.iloc[:, -5], errors="coerce"),
        "te_end": pd.to_numeric(raw.iloc[:, -4], errors="coerce"),
        "te_name": raw.iloc[:, -3].astype(str),
        "te_type": raw.iloc[:, -2],
        "overlap_bp": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
    })

    df = df.dropna(subset=[
        "enh_chrom", "enh_start", "enh_end", "enh_id",
        "te_chrom", "te_start", "te_end", "te_type", "overlap_bp"
    ]).copy()

    df = df[df["te_type"] != "."].copy()
    df = df[df["overlap_bp"] >= MIN_OVERLAP_BP].copy()
    df["te_bin"] = df["te_type"].map(map_te_bin)
    df = df[df["te_bin"].notna()].copy()
    df["tissue"] = df["enh_id"].map(assign_tissue_from_enh_id)

    df["te_id"] = (
        df["te_chrom"].astype(str) + ":" +
        df["te_start"].astype(int).astype(str) + "-" +
        df["te_end"].astype(int).astype(str) + ":" +
        df["te_name"].astype(str) + ":" +
        df["te_type"].astype(str)
    )
    return df


def unique_enhancers(df):
    enh = df[["enh_chrom", "enh_start", "enh_end", "enh_id"]].drop_duplicates().copy()
    enh["win_start"] = (enh["enh_start"] - FLANK_BP).clip(lower=0)
    enh["win_end"] = enh["enh_end"] + FLANK_BP
    return enh


# ----------------------------
# observed and reference counts
# ----------------------------
def observed_enhancer_te_pair_counts(df):
    pairs = df.drop_duplicates(["enh_id", "te_id"]).copy()
    counts = pairs["te_bin"].value_counts().reindex(cat_order, fill_value=0).astype(int)
    return counts


def observed_enhancer_te_pair_composition(df):
    counts = observed_enhancer_te_pair_counts(df).astype(float)
    return counts / counts.sum() * 100


def genome_insert_counts(bg):
    unique_tes = bg.drop_duplicates("te_id").copy()
    return unique_tes["te_bin"].value_counts().reindex(cat_order, fill_value=0).astype(int)


def genome_insert_composition(bg):
    counts = genome_insert_counts(bg).astype(float)
    return counts / counts.sum() * 100


def neighborhood_insert_counts_per_enhancer(enhancers, bg):
    rows = []
    bg_unique = bg.drop_duplicates("te_id").copy()

    for chrom, enh_chr in enhancers.groupby("enh_chrom"):
        te_chr = bg_unique[bg_unique["chrom"] == chrom].copy()
        te_st = te_chr["start"].to_numpy()
        te_en = te_chr["end"].to_numpy()
        te_bin = te_chr["te_bin"].to_numpy()

        for _, row in enh_chr.iterrows():
            ws = row["win_start"]
            we = row["win_end"]
            mask = (te_en > ws) & (te_st < we)

            if not np.any(mask):
                counts = pd.Series(0, index=cat_order, name=row["enh_id"])
                rows.append(counts)
                continue

            sub_bin = te_bin[mask]
            counts = pd.Series(sub_bin).value_counts().reindex(cat_order, fill_value=0).astype(int)
            counts.name = row["enh_id"]
            rows.append(counts)

    out = pd.DataFrame(rows)
    out.index = enhancers["enh_id"].tolist()
    out = out.reindex(columns=cat_order, fill_value=0)
    return out.astype(int)


def neighborhood_total_insert_counts(enhancers, bg):
    per_enh = neighborhood_insert_counts_per_enhancer(enhancers, bg)
    return per_enh.sum(axis=0).reindex(cat_order, fill_value=0).astype(int)


def neighborhood_total_insert_composition(enhancers, bg):
    counts = neighborhood_total_insert_counts(enhancers, bg).astype(float)
    return counts / counts.sum() * 100


# ----------------------------
# statistics
# ----------------------------
def run_gtest_and_chisq(obs_counts, ref_counts, label):
    obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
    ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)

    obs_total = obs_counts.sum()
    ref_total = ref_counts.sum()
    ref_props = ref_counts / ref_total
    expected = ref_props * obs_total
    expected = expected * (obs_total / expected.sum())

    keep = ~((obs_counts == 0) & (expected == 0))
    obs_use = obs_counts[keep]
    exp_use = expected[keep]

    if np.any((exp_use == 0) & (obs_use > 0)):
        summary = pd.DataFrame([{
            "comparison": label,
            "obs_total": int(obs_total),
            "ref_total": int(ref_total),
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
            "observed": obs_counts.values,
            "expected": expected.values,
            "obs_prop": obs_counts.values / obs_total,
            "ref_prop": ref_props.values,
            "chi_contribution": np.nan
        })
        return summary, diagnostic

    chi_stat, chi_p = chisquare(f_obs=obs_use.values, f_exp=exp_use.values)
    g_stat, g_p = power_divergence(f_obs=obs_use.values, f_exp=exp_use.values, lambda_="log-likelihood")

    df = len(obs_use) - 1
    chi_log10_p = chi2_dist.logsf(chi_stat, df) / np.log(10)
    g_log10_p = chi2_dist.logsf(g_stat, df) / np.log(10)

    diagnostic = pd.DataFrame({
        "comparison": label,
        "te_class": obs_use.index,
        "observed": obs_use.values,
        "expected": exp_use.values,
        "obs_prop": obs_use.values / obs_total,
        "ref_prop": ref_props[obs_use.index].values,
        "chi_contribution": ((obs_use.values - exp_use.values) ** 2) / exp_use.values
    })

    summary = pd.DataFrame([{
        "comparison": label,
        "obs_total": int(obs_total),
        "ref_total": int(ref_total),
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


def per_class_insert_enrichment(obs_counts, ref_counts, group_name):
    obs_counts = pd.Series(obs_counts, index=cat_order, dtype=float)
    ref_counts = pd.Series(ref_counts, index=cat_order, dtype=float)

    obs_total = obs_counts.sum()
    ref_total = ref_counts.sum()
    rows = []

    for cat in cat_order:
        a = obs_counts[cat]
        b = obs_total - a
        c = ref_counts[cat]
        d = ref_total - c

        obs_prop = a / obs_total
        ref_prop = c / ref_total
        log2_enrichment = np.log2((obs_prop + EPS) / (ref_prop + EPS))
        log2_or = np.log2(((a + 0.5) / (b + 0.5)) / ((c + 0.5) / (d + 0.5)))

        chi2, p, _, _ = chi2_contingency([[a, b], [c, d]], correction=False)

        rows.append({
            "group": group_name,
            "te_class": cat,
            "obs_count": int(a),
            "obs_total": int(obs_total),
            "ref_count": int(c),
            "ref_total": int(ref_total),
            "obs_prop": obs_prop,
            "ref_prop": ref_prop,
            "log2_obs_over_ref": log2_enrichment,
            "log2_odds_ratio": log2_or,
            "chi2": chi2,
            "p_value": p,
            "direction": "enriched" if obs_prop > ref_prop else "depleted"
        })

    out = pd.DataFrame(rows)
    out["q_value"] = bh_fdr(out["p_value"].values)
    out["label"] = out["q_value"].map(p_to_star)
    return out


# ----------------------------
# pooled plots
# ----------------------------
def run_pooled_analysis(shadow_df, single_df, genome_df):
    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    shadow_obs_counts = observed_enhancer_te_pair_counts(shadow_df)
    single_obs_counts = observed_enhancer_te_pair_counts(single_df)
    shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
    single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)
    genome_counts = genome_insert_counts(genome_df)

    shadow_obs = shadow_obs_counts / shadow_obs_counts.sum() * 100
    single_obs = single_obs_counts / single_obs_counts.sum() * 100
    shadow_local = shadow_local_counts / shadow_local_counts.sum() * 100
    single_local = single_local_counts / single_local_counts.sum() * 100
    genome_comp = genome_counts / genome_counts.sum() * 100

    summary = pd.DataFrame({
        "shadow_observed_enhancer_TE_pairs": shadow_obs,
        "shadow_neighborhood_enhancer_window_TE_pairs": shadow_local,
        "single_observed_enhancer_TE_pairs": single_obs,
        "single_neighborhood_enhancer_window_TE_pairs": single_local,
        "genome_TE_inserts": genome_comp
    }).reindex(cat_order)
    summary.to_csv(OUTPUT_FILES / "pooled_optionB_enhancer_TE_pair_summary.tsv", sep="\t")

    stats_list = []
    diag_list = []
    comparisons = [
        (shadow_obs_counts, shadow_local_counts, "Pooled shadow observed enhancer-TE pairs vs shadow neighborhood TE opportunities"),
        (shadow_obs_counts, genome_counts, "Pooled shadow observed enhancer-TE pairs vs genome TE inserts"),
        (single_obs_counts, single_local_counts, "Pooled single observed enhancer-TE pairs vs single neighborhood TE opportunities"),
        (single_obs_counts, genome_counts, "Pooled single observed enhancer-TE pairs vs genome TE inserts"),
        (shadow_obs_counts, single_obs_counts, "Pooled shadow observed enhancer-TE pairs vs single observed enhancer-TE pairs")
    ]

    for obs, ref, label in comparisons:
        s, d = run_gtest_and_chisq(obs, ref, label)
        stats_list.append(s)
        diag_list.append(d)

    pd.concat(stats_list, ignore_index=True).to_csv(
        OUTPUT_FILES / "pooled_optionB_gtest_chisq_stats.tsv", sep="\t", index=False
    )
    pd.concat(diag_list, ignore_index=True).to_csv(
        OUTPUT_FILES / "pooled_optionB_gof_diagnostics_by_TEclass.tsv", sep="\t", index=False
    )

    shadow_estats = per_class_insert_enrichment(shadow_obs_counts, shadow_local_counts, "Shadow")
    single_estats = per_class_insert_enrichment(single_obs_counts, single_local_counts, "Single")
    enrichment_stats = pd.concat([shadow_estats, single_estats], ignore_index=True)
    enrichment_stats.to_csv(
        OUTPUT_FILES / "pooled_optionB_TEclass_enhancer_TE_pairs_vs_neighborhood_BH.tsv",
        sep="\t",
        index=False
    )

    plot_pooled_stacked(summary)
    plot_pooled_enrichment_bars(enrichment_stats)


def plot_pooled_stacked(summary):
    fig, ax = plt.subplots(figsize=STACKED_FIGSIZE)

    all_values = [
        summary["shadow_observed_enhancer_TE_pairs"].values,
        summary["shadow_neighborhood_enhancer_window_TE_pairs"].values,
        summary["single_observed_enhancer_TE_pairs"].values,
        summary["single_neighborhood_enhancer_window_TE_pairs"].values,
        summary["genome_TE_inserts"].values
    ]

    all_labels = [
        "shadows",
        "shadow\nneighborhood",
        "single",
        "single\nneighborhood",
        "genome"
    ]

    x_positions = [0.0, 1.0, 2.6, 3.6, 5.2]

    for xpos, vals in zip(x_positions, all_values):
        plot_stacked_bar(ax, vals, xpos, width=0.62)

    ax.set_ylabel("TE insert %")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(all_labels, fontsize=12)
    ax.set_ylim(0, 100)
    make_grouped_te_legend(ax)
    style_axes(ax)

    save_figure(
        fig,
        png_path=f"{output_dir}/pooled_optionB_stacked_TE_observed_local_genome_insert_based.png",
        pdf_path=f"{output_dir}/pooled_optionB_stacked_TE_observed_local_genome_insert_based.pdf",
        left=0.10,
        right=0.80,
        bottom=0.24,
        top=0.92
    )


def plot_pooled_enrichment_bars(enrichment_stats):
    labels = cat_order
    shadow_vals = enrichment_stats[enrichment_stats["group"] == "Shadow"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()
    single_vals = enrichment_stats[enrichment_stats["group"] == "Single"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()

    x = np.arange(len(cat_order))
    bar_width = 0.46

    fig, ax = plt.subplots(figsize=ENRICHMENT_BAR_FIGSIZE)

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
    ax.set_xticklabels(labels, fontsize=12)
    if PLOT_METRIC == "log2_obs_over_ref":
        ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=12)
    else:
        ax.set_ylabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=12)
    ax.tick_params(axis="y", labelsize=11)
    ax.tick_params(axis="x", labelsize=11)
    apply_enrichment_axis(ax)

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if np.isfinite(height):
                offset = 0.10 if height >= 0 else -0.20
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    height + offset,
                    f"{height:.2f}",
                    ha="center",
                    va="bottom" if height >= 0 else "top",
                    fontsize=11
                )

    for i, cat in enumerate(cat_order):
        sh_row = enrichment_stats[(enrichment_stats["group"] == "Shadow") & (enrichment_stats["te_class"] == cat)].iloc[0]
        si_row = enrichment_stats[(enrichment_stats["group"] == "Single") & (enrichment_stats["te_class"] == cat)].iloc[0]

        if sh_row["label"] != "ns":
            ax.text(
                x[i] - bar_width / 2,
                shadow_vals[i] + 0.28,
                sh_row["label"],
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight="bold"
            )

        if si_row["label"] != "ns":
            ax.text(
                x[i] + bar_width / 2,
                single_vals[i] + 0.28,
                si_row["label"],
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight="bold"
            )

    shadow_patch = mpatches.Patch(facecolor="white", edgecolor="black", hatch="//", label="shadows")
    single_patch = mpatches.Patch(facecolor="white", edgecolor="black", label="singles")
    ax.legend(handles=[shadow_patch, single_patch], fontsize=11, loc="upper right", frameon=False)
    style_axes(ax)

    save_figure(
        fig,
        png_path=f"{output_dir}/pooled_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.png",
        pdf_path=f"{output_dir}/pooled_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.pdf",
        left=0.16,
        right=0.98,
        bottom=0.24,
        top=0.92
    )


# ----------------------------
# tissue-specific analysis
# ----------------------------
def run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df):
    shadow_df = shadow_df_all[shadow_df_all["tissue"] == tissue].copy()
    single_df = single_df_all[single_df_all["tissue"] == tissue].copy()

    if shadow_df.empty or single_df.empty:
        print(f"Skipping {tissue}: missing data in one group")
        return None

    shadow_enh = unique_enhancers(shadow_df)
    single_enh = unique_enhancers(single_df)

    shadow_obs_counts = observed_enhancer_te_pair_counts(shadow_df)
    single_obs_counts = observed_enhancer_te_pair_counts(single_df)
    shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
    single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)
    genome_counts = genome_insert_counts(genome_df)

    shadow_obs = shadow_obs_counts / shadow_obs_counts.sum() * 100
    single_obs = single_obs_counts / single_obs_counts.sum() * 100
    shadow_local = shadow_local_counts / shadow_local_counts.sum() * 100
    single_local = single_local_counts / single_local_counts.sum() * 100
    genome_comp = genome_counts / genome_counts.sum() * 100

    summary = pd.DataFrame({
        "shadow_observed_enhancer_TE_pairs": shadow_obs,
        "shadow_neighborhood_enhancer_window_TE_pairs": shadow_local,
        "single_observed_enhancer_TE_pairs": single_obs,
        "single_neighborhood_enhancer_window_TE_pairs": single_local,
        "genome_TE_inserts": genome_comp
    }).reindex(cat_order)

    summary.to_csv(OUTPUT_FILES / f"{tissue}_optionB_enhancer_TE_pair_summary.tsv", sep="\t")

    stats_list = []
    diag_list = []
    comparisons = [
        (shadow_obs_counts, shadow_local_counts, f"{tissue} shadow observed enhancer-TE pairs vs shadow neighborhood TE opportunities"),
        (shadow_obs_counts, genome_counts, f"{tissue} shadow observed enhancer-TE pairs vs genome TE inserts"),
        (single_obs_counts, single_local_counts, f"{tissue} single observed enhancer-TE pairs vs single neighborhood TE opportunities"),
        (single_obs_counts, genome_counts, f"{tissue} single observed enhancer-TE pairs vs genome TE inserts"),
    ]

    for obs, ref, label in comparisons:
        s, d = run_gtest_and_chisq(obs, ref, label)
        stats_list.append(s)
        diag_list.append(d)

    pd.concat(stats_list, ignore_index=True).to_csv(
        OUTPUT_FILES / f"{tissue}_optionB_gtest_chisq_stats.tsv", sep="\t", index=False
    )
    pd.concat(diag_list, ignore_index=True).to_csv(
        OUTPUT_FILES / f"{tissue}_optionB_gof_diagnostics_by_TEclass.tsv", sep="\t", index=False
    )

    shadow_estats = per_class_insert_enrichment(shadow_obs_counts, shadow_local_counts, "Shadow")
    single_estats = per_class_insert_enrichment(single_obs_counts, single_local_counts, "Single")
    enrichment_stats = pd.concat([shadow_estats, single_estats], ignore_index=True)

    enrichment_stats.to_csv(
        OUTPUT_FILES / f"{tissue}_optionB_TEclass_enhancer_TE_pairs_vs_neighborhood_BH.tsv",
        sep="\t",
        index=False
    )

    plot_tissue_stacked(tissue, summary)
    plot_tissue_enrichment_bars(tissue, enrichment_stats)
    return enrichment_stats


def plot_tissue_stacked(tissue, summary):
    fig, ax = plt.subplots(figsize=STACKED_FIGSIZE)

    all_values = [
        summary["shadow_observed_enhancer_TE_pairs"].values,
        summary["shadow_neighborhood_enhancer_window_TE_pairs"].values,
        summary["single_observed_enhancer_TE_pairs"].values,
        summary["single_neighborhood_enhancer_window_TE_pairs"].values,
        summary["genome_TE_inserts"].values
    ]
    all_labels = [
        "shadows",
        "shadow\nneighborhood",
        "single",
        "single\nneighborhood",
        "genome"
    ]
    x_positions = [0.0, 1.0, 2.6, 3.6, 5.2]

    for xpos, vals in zip(x_positions, all_values):
        plot_stacked_bar(ax, vals, xpos, width=0.62)

    ax.set_title(tissue.capitalize(), fontsize=13)
    ax.set_ylabel("TE insert %")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(all_labels, fontsize=12)
    ax.set_ylim(0, 100)
    make_grouped_te_legend(ax)
    style_axes(ax)

    save_figure(
        fig,
        png_path=f"{output_dir}/{tissue}_optionB_stacked_TE_observed_local_genome_insert_based.png",
        pdf_path=f"{output_dir}/{tissue}_optionB_stacked_TE_observed_local_genome_insert_based.pdf",
        left=0.10,
        right=0.80,
        bottom=0.24,
        top=0.90
    )


def plot_tissue_enrichment_bars(tissue, enrichment_stats):
    shadow_vals = enrichment_stats[enrichment_stats["group"] == "Shadow"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()
    single_vals = enrichment_stats[enrichment_stats["group"] == "Single"].set_index("te_class").reindex(cat_order)[PLOT_METRIC].to_numpy()

    x = np.arange(len(cat_order))
    bar_width = 0.46
    fig, ax = plt.subplots(figsize=ENRICHMENT_BAR_FIGSIZE)

    bars1 = ax.bar(
        x - bar_width / 2,
        shadow_vals,
        width=bar_width,
        color=[te_color_map[lbl] for lbl in cat_order],
        edgecolor="black",
        hatch="//"
    )
    bars2 = ax.bar(
        x + bar_width / 2,
        single_vals,
        width=bar_width,
        color=[te_color_map[lbl] for lbl in cat_order],
        edgecolor="black"
    )

    ax.set_title(tissue.capitalize(), fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(cat_order, fontsize=12)
    if PLOT_METRIC == "log2_obs_over_ref":
        ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=11)
    else:
        ax.set_ylabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=12)
    apply_enrichment_axis(ax)

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if np.isfinite(height):
                offset = 0.10 if height >= 0 else -0.20
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    height + offset,
                    f"{height:.2f}",
                    ha="center",
                    va="bottom" if height >= 0 else "top",
                    fontsize=11
                )

    for i, cat in enumerate(cat_order):
        sh_row = enrichment_stats[(enrichment_stats["group"] == "Shadow") & (enrichment_stats["te_class"] == cat)].iloc[0]
        si_row = enrichment_stats[(enrichment_stats["group"] == "Single") & (enrichment_stats["te_class"] == cat)].iloc[0]

        if sh_row["label"] != "ns":
            ax.text(x[i] - bar_width / 2, shadow_vals[i] + 0.28, sh_row["label"], ha="center", va="bottom", fontsize=11, fontweight="bold")
        if si_row["label"] != "ns":
            ax.text(x[i] + bar_width / 2, single_vals[i] + 0.28, si_row["label"], ha="center", va="bottom", fontsize=11, fontweight="bold")

    shadow_patch = mpatches.Patch(facecolor="white", edgecolor="black", hatch="//", label="shadows")
    single_patch = mpatches.Patch(facecolor="white", edgecolor="black", label="singles")
    ax.legend(handles=[shadow_patch, single_patch], fontsize=11, loc="upper right", frameon=False)
    style_axes(ax)

    save_figure(
        fig,
        png_path=f"{output_dir}/{tissue}_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.png",
        pdf_path=f"{output_dir}/{tissue}_optionB_{PLOT_METRIC}_observed_over_neighborhood_by_TEtype.pdf",
        left=0.16,
        right=0.98,
        bottom=0.24,
        top=0.90
    )


# ----------------------------
# tissue enrichment summary tables
# ----------------------------
def compute_tissue_enrichment_tables(shadow_df_all, single_df_all, genome_df):
    rows = []

    for tissue in tissue_order:
        shadow_df = shadow_df_all[shadow_df_all["tissue"] == tissue].copy()
        single_df = single_df_all[single_df_all["tissue"] == tissue].copy()

        if shadow_df.empty or single_df.empty:
            print(f"Skipping {tissue}: no data in one of the groups")
            continue

        shadow_enh = unique_enhancers(shadow_df)
        single_enh = unique_enhancers(single_df)

        shadow_obs_counts = observed_enhancer_te_pair_counts(shadow_df)
        single_obs_counts = observed_enhancer_te_pair_counts(single_df)
        shadow_local_counts = neighborhood_total_insert_counts(shadow_enh, genome_df)
        single_local_counts = neighborhood_total_insert_counts(single_enh, genome_df)

        shadow_stats = per_class_insert_enrichment(shadow_obs_counts, shadow_local_counts, group_name="shadows").copy()
        shadow_stats["tissue"] = tissue

        single_stats = per_class_insert_enrichment(single_obs_counts, single_local_counts, group_name="singles").copy()
        single_stats["tissue"] = tissue

        rows.append(shadow_stats)
        rows.append(single_stats)

    shadow_enh_all = unique_enhancers(shadow_df_all)
    single_enh_all = unique_enhancers(single_df_all)
    shadow_obs_counts_all = observed_enhancer_te_pair_counts(shadow_df_all)
    single_obs_counts_all = observed_enhancer_te_pair_counts(single_df_all)
    shadow_local_counts_all = neighborhood_total_insert_counts(shadow_enh_all, genome_df)
    single_local_counts_all = neighborhood_total_insert_counts(single_enh_all, genome_df)

    shadow_stats_all = per_class_insert_enrichment(shadow_obs_counts_all, shadow_local_counts_all, group_name="shadows").copy()
    shadow_stats_all["tissue"] = "pooled"

    single_stats_all = per_class_insert_enrichment(single_obs_counts_all, single_local_counts_all, group_name="singles").copy()
    single_stats_all["tissue"] = "pooled"

    rows.append(shadow_stats_all)
    rows.append(single_stats_all)

    out = pd.concat(rows, ignore_index=True)
    out["te_class"] = pd.Categorical(out["te_class"], categories=cat_order, ordered=True)
    out["tissue"] = pd.Categorical(out["tissue"], categories=tissue_order + ["pooled"], ordered=True)
    return out.sort_values(["group", "te_class", "tissue"])


# ----------------------------
# tissue split dot plot
# ----------------------------
def plot_tissue_dotplot(enrichment_df, metric="log2_obs_over_ref", outfile_png=None, outfile_pdf=None):
    fig, axes = plt.subplots(1, 2, figsize=TISSUE_SPLIT_FIGSIZE, sharey=True)

    group_order = ["shadows", "singles"]
    group_titles = {"shadows": "Shadows", "singles": "Singles"}
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
                    x,
                    y,
                    s=55,
                    color=tissue_colors[tissue],
                    edgecolor="black",
                    linewidth=0.8,
                    zorder=3
                )

                label = str(row.iloc[0]["label"])
                if label != "ns":
                    ax.text(
                        x + 0.08,
                        y,
                        label,
                        ha="left",
                        va="center",
                        fontsize=9,
                        fontweight="bold"
                    )

        pooled_sub = sub[sub["tissue"] == "pooled"].copy()
        for i, te in enumerate(cat_order):
            row = pooled_sub[pooled_sub["te_class"] == te]
            if row.empty:
                continue

            x = float(row.iloc[0][metric])
            y = i + tissue_yoffset["pooled"]
            ax.scatter(
                x,
                y,
                s=75,
                color=tissue_colors["pooled"],
                edgecolor="black",
                linewidth=0.8,
                marker="D",
                zorder=4
            )

        ax.set_title(group_titles[group], fontsize=14)
        ax.set_xlim(*ENRICHMENT_YLIM)
        ax.set_yticks(yvals)
        ax.set_yticklabels(cat_order, fontsize=12)
        ax.invert_yaxis()
        style_axes(ax)

    if metric == "log2_obs_over_ref":
        axes[0].set_ylabel("TE class", fontsize=12)
        fig.supxlabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=13, y=0.08)
    else:
        axes[0].set_ylabel("TE class", fontsize=12)
        fig.supxlabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=13, y=0.08)

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["forebrain"], markeredgecolor="black", markersize=7, label="Forebrain"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["heart"], markeredgecolor="black", markersize=7, label="Heart"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["limb"], markeredgecolor="black", markersize=7, label="Limb"),
        Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Pooled")
    ]

    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=4, frameon=False, fontsize=11)

    save_figure(fig, png_path=outfile_png, pdf_path=outfile_pdf, left=0.10, right=0.98, bottom=0.24, top=0.82)


# ----------------------------
# final combined vertical dot plot
# ----------------------------
def plot_combined_vertical_dotplot_clean(
    enrichment_df,
    metric="log2_obs_over_ref",
    outfile_png=None,
    outfile_pdf=None,
    show_pooled=True,
    show_stars=True
):
    fig, ax = plt.subplots(figsize=COMBINED_DOTPLOT_FIGSIZE)

    # Horizontal layout controls.
    # Tissue offsets are spread out to prevent overlap within each TE class.
    # Shadow/single offsets are intentionally small so each pair stays close.
    category_spacing = 1.45
    x_base = np.arange(len(cat_order)) * category_spacing

    tissue_xoffset = {
        "forebrain": -0.40,
        "heart": 0.00,
        "limb": 0.40,
        "pooled": 0.66
    }
    group_xoffset = {
        "shadows": -0.035,
        "singles": 0.035
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
                (enrichment_df["group"] == "shadows") &
                (enrichment_df["te_class"] == te)
            ]
            si = enrichment_df[
                (enrichment_df["tissue"] == tissue) &
                (enrichment_df["group"] == "singles") &
                (enrichment_df["te_class"] == te)
            ]

            if sh.empty or si.empty:
                continue

            x_sh = x_base[te_i] + tissue_xoffset[tissue] + group_xoffset["shadows"]
            x_si = x_base[te_i] + tissue_xoffset[tissue] + group_xoffset["singles"]
            y_sh = float(sh.iloc[0][metric])
            y_si = float(si.iloc[0][metric])

            # No connector line between shadow and single points.
            ax.scatter(
                x_sh,
                y_sh,
                s=48,
                facecolor=color,
                edgecolor="black",
                linewidth=0.9,
                zorder=4
            )
            ax.scatter(
                x_si,
                y_si,
                s=48,
                facecolor="white",
                edgecolor=color,
                linewidth=1.8,
                zorder=4
            )

            if show_stars:
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
                            fontsize=8,
                            fontweight="bold",
                            zorder=6,
                            clip_on=False
                        )

    if show_pooled:
        for te_i, te in enumerate(cat_order):
            for group in ["shadows", "singles"]:
                row = enrichment_df[
                    (enrichment_df["tissue"] == "pooled") &
                    (enrichment_df["group"] == group) &
                    (enrichment_df["te_class"] == te)
                ]
                if row.empty:
                    continue

                x = x_base[te_i] + tissue_xoffset["pooled"] + group_xoffset[group]
                y = float(row.iloc[0][metric])

                if group == "shadows":
                    ax.scatter(x, y, s=70, marker="D", facecolor="black", edgecolor="black", linewidth=1.0, zorder=5)
                else:
                    ax.scatter(x, y, s=70, marker="D", facecolor="white", edgecolor="black", linewidth=1.8, zorder=5)

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
                            clip_on=False
                        )

    ax.set_xticks(x_base)
    ax.set_xticklabels(cat_order, fontsize=13)

    # Trim the unused left-side whitespace so the first TE-class points start closer
    # to the y-axis without clipping the leftmost marker or star.
    leftmost_point = x_base[0] + min(tissue_xoffset.values()) + min(group_xoffset.values())
    rightmost_point = x_base[-1] + max(tissue_xoffset.values()) + max(group_xoffset.values())
    ax.set_xlim(leftmost_point - 0.14, rightmost_point + 0.18)

    if metric == "log2_obs_over_ref":
        ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Neighborhood}}\right)$', fontsize=13)
    else:
        ax.set_ylabel(r'$\log_2(\mathrm{odds\ ratio})$', fontsize=13)

    apply_enrichment_axis(ax)
    style_axes(ax)

    tissue_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["forebrain"], markeredgecolor="black", markersize=7, label="Forebrain"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["heart"], markeredgecolor="black", markersize=7, label="Heart"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=tissue_colors["limb"], markeredgecolor="black", markersize=7, label="Limb")
    ]
    if show_pooled:
        tissue_handles.append(
            Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Pooled")
        )

    class_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="black", markeredgecolor="black", markersize=7, label="Shadows"),
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="white", markeredgecolor="black", markersize=7, label="Singles")
    ]

    leg1 = ax.legend(
        handles=tissue_handles,
        title="Tissue",
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00)
    )
    ax.add_artist(leg1)

    ax.legend(
        handles=class_handles,
        title="Enhancer class",
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.50)
    )

    save_figure(fig, png_path=outfile_png, pdf_path=outfile_pdf, left=0.10, right=0.77, bottom=0.24, top=0.92)


# ----------------------------
# run
# ----------------------------
def main():
    genome_df = load_background(genome_file)
    shadow_df_all = load_mouse_overlap(shadow_file)
    single_df_all = load_mouse_overlap(single_file)

    shadow_df_all = shadow_df_all[shadow_df_all["tissue"].isin(tissue_order)].copy()
    single_df_all = single_df_all[single_df_all["tissue"].isin(tissue_order)].copy()

    print("Shadow tissues:")
    print(shadow_df_all["tissue"].value_counts())
    print("\nSingle tissues:")
    print(single_df_all["tissue"].value_counts())

    run_pooled_analysis(shadow_df_all, single_df_all, genome_df)

    for tissue in tissue_order:
        run_tissue_analysis(tissue, shadow_df_all, single_df_all, genome_df)

    enrichment_df = compute_tissue_enrichment_tables(
        shadow_df_all=shadow_df_all,
        single_df_all=single_df_all,
        genome_df=genome_df
    )

    print("\nCombined tissue enrichment table:")
    print(enrichment_df.round(4))

    enrichment_df.to_csv(
        OUTPUT_FILES / "tissue_dotplot_enrichment_table.tsv",
        sep="\t",
        index=False
    )



    plot_combined_vertical_dotplot_clean(
        enrichment_df,
        metric=PLOT_METRIC,
        outfile_png=f"{output_dir}/combined_vertical_dotplot_TEclass_enrichment_clean.png",
        outfile_pdf=f"{output_dir}/combined_vertical_dotplot_TEclass_enrichment_clean.pdf",
        show_pooled=True,
        show_stars=True
    )


if __name__ == "__main__":
    main()
