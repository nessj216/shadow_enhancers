import collections
import matplotlib.pyplot as plt
import numpy as np

##############################################################################
# 1) Load & bucket SHADOW file
##############################################################################
shadow_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/011925_all_shadowsets_DM6.bed"

shadow_counts = collections.Counter()
with open(shadow_file, 'r') as sf:
    for line in sf:
        if line.strip():
            cols = line.strip().split('\t')
            set_name = cols[3]  # the "column 4" ID
            shadow_counts[set_name] += 1

shadow_buckets = {}
for set_name, count in shadow_counts.items():
    if count == 2:
        shadow_buckets[set_name] = "2 shadows"
    elif 3 <= count <= 4:
        shadow_buckets[set_name] = "3-4 shadows"
    elif count >= 5:
        shadow_buckets[set_name] = ">=5 shadows"
    # else: handle 1 or 0 shadows if you wish

##############################################################################
# 2) Load & bucket TE file
##############################################################################
te_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/overlap_enhancer_TE.bed"

te_counts = collections.Counter()
with open(te_file, 'r') as tf:
    for line in tf:
        if line.strip():
            cols = line.strip().split('\t')
            set_name = cols[3]
            te_counts[set_name] += 1

te_buckets = {}
for set_name, count in te_counts.items():
    if count == 1:
        te_buckets[set_name] = "1 TE"
    elif count == 2:
        te_buckets[set_name] = "2 TE"
    elif 3 <= count <= 4:
        te_buckets[set_name] = "3-4 TE"
    else:  # count >= 5
        te_buckets[set_name] = "5+ TE"

##############################################################################
# 3) Cross-tabulate shadow vs. TE (raw counts + store set names)
##############################################################################
shadow_categories = ["2 shadows", "3-4 shadows", ">=5 shadows"]
te_categories     = ["1 TE", "2 TE", "3-4 TE", "5+ TE"]

# combo_counts holds the raw integer counts
combo_counts = {}
# sets_in_bin holds which set_names fall in each bin
from collections import defaultdict
sets_in_bin = defaultdict(list)

# Initialize counts to 0
for s_cat in shadow_categories:
    for t_cat in te_categories:
        combo_counts[(s_cat, t_cat)] = 0

# Fill combo_counts and sets_in_bin
for set_name, s_bucket in shadow_buckets.items():
    if set_name in te_buckets:
        t_bucket = te_buckets[set_name]
        combo_counts[(s_bucket, t_bucket)] += 1
        sets_in_bin[(s_bucket, t_bucket)].append(set_name)

##############################################################################
# 4) Convert to row-wise fractions where denominator = sum of TE bins in row
##############################################################################
heatmap_data = []
for s_cat in shadow_categories:
    row_te_sum = sum(combo_counts[(s_cat, t_cat)] for t_cat in te_categories)
    row_fractions = []
    for t_cat in te_categories:
        numerator = combo_counts[(s_cat, t_cat)]
        fraction = 0.0
        if row_te_sum > 0:
            fraction = numerator / row_te_sum
        row_fractions.append(fraction)
    heatmap_data.append(row_fractions)

##############################################################################
# 5) Print numerator, denominator, fraction, and the set names in each bin
##############################################################################
print("Detailed bin info (Shadow x TE):\n")
for row_idx, s_cat in enumerate(shadow_categories):
    row_te_sum = sum(combo_counts[(s_cat, t_cat)] for t_cat in te_categories)
    print(f"Shadow bucket: {s_cat}  (row_te_sum = {row_te_sum})")
    for col_idx, t_cat in enumerate(te_categories):
        numerator = combo_counts[(s_cat, t_cat)]
        fraction = heatmap_data[row_idx][col_idx]
        # The list of sets in this bin
        set_names_here = sets_in_bin[(s_cat, t_cat)]
        set_names_str = ", ".join(set_names_here) if set_names_here else "None"
        print(f"  {t_cat}: {numerator}/{row_te_sum} = {fraction:.3f}")
        print(f"    Sets in this bin: {set_names_str}")
    print()

##############################################################################
# 6) Plot the heatmap
##############################################################################
##############################################################################
# 6) Plot the heatmap and add row sums to each Y-axis label
##############################################################################
fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(heatmap_data, cmap="Blues", aspect='auto')

# X-axis: TE categories
ax.set_xticks(np.arange(len(te_categories)))
ax.set_xticklabels(te_categories)

# Y-axis: Shadow categories with row sums
y_labels = []
for s_cat in shadow_categories:
    # Compute how many sets total in this shadow category (sum across all TE bins)
    row_te_sum = sum(combo_counts[(s_cat, t_cat)] for t_cat in te_categories)
    # Add a line break to put the sum underneath
    y_labels.append(f"{s_cat}\nn={row_te_sum}")

ax.set_yticks(np.arange(len(shadow_categories)))
ax.set_yticklabels(y_labels)

ax.set_xlabel("TE-Derived Categories")
ax.set_ylabel("Total Shadows in Set")
ax.set_title("Row-Wise Proportions (Denominator = Sum of TE bins in each shadow row)")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Proportion of TE-Labeled Sets", rotation=270, labelpad=15)

plt.tight_layout()
plt.show()
