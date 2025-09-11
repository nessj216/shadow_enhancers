'''The denominator across the row is the total shadows whether or not they had a +TE enhancer'''
import collections
import matplotlib.pyplot as plt
import numpy as np

##############################################################################
# 1) Read & bucket SHADOW file
#    Bucket sets based on how many times they appear in the "shadow" file:
#       - exactly 2          -> "2 shadows"
#       - between 3 and 4    -> "3-4 shadows"
#       - >= 5               -> ">=5 shadows"
##############################################################################
shadow_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/all_shadowsets_DM6.bed"

shadow_counts = collections.Counter()

with open(shadow_file, 'r') as sf:
    for line in sf:
        if line.strip():
            cols = line.strip().split('\t')
            # Assume set name is in column 4 (index 3)
            set_name = cols[3]
            shadow_counts[set_name] += 1

# Assign each set to a shadow bucket
shadow_buckets = {}
for set_name, count in shadow_counts.items():
    if count == 2:
        shadow_buckets[set_name] = "2 shadows"
    elif 3 <= count <= 4:
        shadow_buckets[set_name] = "3-4 shadows"
    elif count >= 5:
        shadow_buckets[set_name] = ">=5 shadows"
    # If you have sets with 1 or 0 shadows, decide how to handle or skip them

# For convenience, count how many sets fall into each shadow category
shadow_cat_counts = collections.Counter(shadow_buckets.values())
# e.g. shadow_cat_counts["2 shadows"] = number of sets that have exactly 2 shadows

##############################################################################
# 2) Read & bucket TE file
#    Bucket sets based on how many TEs they have:
#       - exactly 1          -> "1 TE"
#       - exactly 2          -> "2 TE"
#       - between 3 and 4    -> "3-4 TE"
#       - >= 5               -> "5+ TE"
##############################################################################
te_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/cooption_enhancer_TE.bed"

te_counts = collections.Counter()

with open(te_file, 'r') as tf:
    for line in tf:
        if line.strip():
            cols = line.strip().split('\t')
            # Assume set name is in column 4 (index 3)
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
    else:
        # count >= 5
        te_buckets[set_name] = "5+ TE"

##############################################################################
# 3) Cross‐tabulate into a 3×4 table: (shadow bucket) × (TE bucket)
##############################################################################
shadow_categories = ["2 shadows", "3-4 shadows", ">=5 shadows"]
te_categories     = ["1 TE", "2 TE", "3-4 TE", "5+ TE"]

# combo_counts[(s_cat, t_cat)] = how many sets are in shadow bucket s_cat + TE bucket t_cat
combo_counts = {}
for s_cat in shadow_categories:
    for t_cat in te_categories:
        combo_counts[(s_cat, t_cat)] = 0

for set_name, s_bucket in shadow_buckets.items():
    if set_name in te_buckets:
        t_bucket = te_buckets[set_name]
        combo_counts[(s_bucket, t_bucket)] += 1
    # If a set doesn't appear in te_buckets, it won't increment any TE bucket
    # That effectively counts as "no TE" or "not TE+" and won't show in the matrix.

##############################################################################
# 4) Convert raw counts to PROPORTIONS:
#    Each cell is: combo_counts[(s_cat, t_cat)] / total # sets in s_cat
#    where "total # sets in s_cat" = shadow_cat_counts[s_cat] from shadow_file
##############################################################################
heatmap_data = []
for s_cat in shadow_categories:
    row = []
    row_total = shadow_cat_counts[s_cat]  # total sets in this shadow bucket
    for t_cat in te_categories:
        num_in_both = combo_counts[(s_cat, t_cat)]
        prop = 0.0
        if row_total > 0:
            prop = num_in_both / row_total
        row.append(prop)
    heatmap_data.append(row)

# At this point, each row's sum will be <= 1.0 (often less), because
# not all sets in that shadow category necessarily appear in the TE file.

##############################################################################
# 5) Plot the PROPORTION heatmap
#    - X‐axis = TE categories (1 TE, 2 TE, 3-4 TE, 5+ TE)
#    - Y‐axis = Shadow categories (2 shadows, 3-4 shadows, >=5 shadows)
##############################################################################

print("\nProportion matrix (rows = shadow categories, cols = TE categories):")
for s_cat, row_vals in zip(shadow_categories, heatmap_data):
    print(f"{s_cat}\t{row_vals}")


# Suppose you have already defined:
# shadow_categories = ["2 shadows", "3-4 shadows", ">=5 shadows"]
# te_categories = ["1 TE", "2 TE", "3-4 TE", "5+ TE"]
# heatmap_data = ... (a list of lists: 3 rows × 4 columns)

print("Heatmap fractions (row = Shadow Category, col = TE Category):")
for row_idx, s_cat in enumerate(shadow_categories):
    for col_idx, t_cat in enumerate(te_categories):
        fraction = heatmap_data[row_idx][col_idx]
        print(f"  {s_cat} | {t_cat} -> {fraction:.3f}")
    print()  # blank line after each shadow category row

print("Detailed counts per box:")
for s_cat in shadow_categories:
    denom = shadow_cat_counts[s_cat]  # total sets in this shadow bucket
    for t_cat in te_categories:
        numerator = combo_counts[(s_cat, t_cat)]
        fraction = 0.0
        if denom > 0:
            fraction = numerator / denom
        print(f"  {s_cat} & {t_cat}:  {numerator}/{denom} = {fraction:.3f}")
    print()  # blank line after each shadow category





fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(heatmap_data, cmap="Blues", aspect='auto')

# Set up axis ticks
ax.set_xticks(np.arange(len(te_categories)))
ax.set_yticks(np.arange(len(shadow_categories)))
ax.set_xticklabels(te_categories)
ax.set_yticklabels(shadow_categories)

ax.set_xlabel("TE‐Derived Categories")
ax.set_ylabel("Total Shadows in Set")
ax.set_title("Proportion of Sets in Each Shadow Category by TE Category")

# Colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Proportion of Sets", rotation=270, labelpad=15)

plt.tight_layout()
plt.show()

##############################################################################
# 6) (Optional) Print the proportion matrix for debugging
##############################################################################
print("\nProportion matrix (rows = shadow categories, cols = TE categories):")
for s_cat, row_vals in zip(shadow_categories, heatmap_data):
    print(f"{s_cat}\t{row_vals}")


# Suppose you have already defined:
# shadow_categories = ["2 shadows", "3-4 shadows", ">=5 shadows"]
# te_categories = ["1 TE", "2 TE", "3-4 TE", "5+ TE"]
# heatmap_data = ... (a list of lists: 3 rows × 4 columns)

print("Heatmap fractions (row = Shadow Category, col = TE Category):")
for row_idx, s_cat in enumerate(shadow_categories):
    for col_idx, t_cat in enumerate(te_categories):
        fraction = heatmap_data[row_idx][col_idx]
        print(f"  {s_cat} | {t_cat} -> {fraction:.3f}")
    print()  # blank line after each shadow category row
