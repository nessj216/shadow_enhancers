import collections
import matplotlib.pyplot as plt
import numpy as np

##############################################################################
# 1) Read & compute the number of PAIRS for each set in your "shadow" file
#    Pairs = n*(n-1)/2 if a set has n shadows.
#    Then bin sets by exactly that number of pairs.
##############################################################################
shadow_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/011925_all_shadowsets_DM6.bed"

shadow_counts = collections.Counter()

# First, count how many shadows each set has
with open(shadow_file, 'r') as sf:
    for line in sf:
        if line.strip():
            cols = line.strip().split('\t')
            # Suppose set name is in column 4 (index 3)
            set_name = cols[3]
            shadow_counts[set_name] += 1

# Now map each set to the bin = total number of pairs
# If a set has n shadows, the number of pairs is n*(n-1)/2.
shadow_pairs = {}
for set_name, n_shadows in shadow_counts.items():
    # Option 1: If you want to *exclude* sets that have fewer than 2 shadows, do:
    if n_shadows < 2:
        continue

    # Option 2: If you want to *include* them in a "0 pairs" bin, you can do:
    # pairs = (n_shadows * (n_shadows - 1)) // 2
    # shadow_pairs[set_name] = pairs
    # (and remove the "if n_shadows < 2: continue" above)

    # For now, let's exclude sets that have < 2 shadows:
    pairs = (n_shadows * (n_shadows - 1)) // 2
    shadow_pairs[set_name] = pairs

# Count how many sets fall into each *pairs* bin
pairs_counts = collections.Counter(shadow_pairs.values())

# We'll make a sorted list of all distinct pairs values
all_pairs = sorted(pairs_counts.keys())  # e.g. [1, 3, 6, 10, ...]

##############################################################################
# 2) Read & bucket TE file
#    Bucket sets based on how many BLAST hits they have:
#       - exactly 1     -> "1 BLAST"
#       - exactly 2     -> "2 BLAST"
#       - exactly 3     -> "3 BLAST"
#       - >= 4          -> "4+ BLAST"
##############################################################################
te_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/HEATMAP_stuff/BLAStHITS_deduplicated_lines.csv"

te_counts = collections.Counter()

with open(te_file, 'r') as tf:
    for line in tf:
        if line.strip():
            cols = line.strip().split(',')
            # Suppose set name is in column 0
            set_name = cols[0]
            te_counts[set_name] += 1


def get_te_bucket(count):
    if count == 1:
        return "1 BLAST"
    elif count == 2:
        return "2 BLAST"
    elif count == 3:
        return "3 BLAST"
    else:
        return "4+ BLAST"


te_buckets = {}
for set_name, count in te_counts.items():
    te_buckets[set_name] = get_te_bucket(count)

##############################################################################
# 3) Create a cross‐tab for (#pairs bin) × (TE bucket), i.e. how many sets
#    simultaneously fall into a given pairs bin and a given TE category.
##############################################################################
te_categories = ["1 BLAST", "2 BLAST", "3 BLAST", "4+ BLAST"]

# Initialize a dictionary to hold raw counts for each (pairs_bin, TE category)
combo_counts = {}
for p_bin in all_pairs:
    for t_cat in te_categories:
        combo_counts[(p_bin, t_cat)] = 0

# Populate combo_counts
for set_name, p_bin in shadow_pairs.items():
    # Check if we have a TE bucket for this set. If not, it means 0 BLAST (?)
    # Decide if you want to skip or treat them as "No BLAST" category, etc.
    if set_name in te_buckets:
        t_bucket = te_buckets[set_name]
        combo_counts[(p_bin, t_bucket)] += 1
    # else:
    #    # Possibly count them in a "0 BLAST" bin or skip them
    #    pass

##############################################################################
# 4) Convert raw counts to PROPORTIONS:
#    For each pair bin = p_bin, total # sets in that bin = pairs_counts[p_bin].
#    Then each cell in the row is:
#         combo_counts[(p_bin, t_cat)] / pairs_counts[p_bin]
##############################################################################
heatmap_data = []
for p_bin in all_pairs:
    row = []
    row_total = pairs_counts[p_bin]  # how many sets have this # of pairs
    for t_cat in te_categories:
        num_in_both = combo_counts[(p_bin, t_cat)]
        prop = 0.0
        if row_total > 0:
            prop = num_in_both / row_total
        row.append(prop)
    heatmap_data.append(row)

##############################################################################
# 5) Plot the PROPORTION heatmap
#    - X-axis = TE categories
#    - Y-axis = the distinct # of pairs bins
##############################################################################
print("\nProportion matrix (rows = # of pairs bin, cols = TE categories):")
for p_bin, row_vals in zip(all_pairs, heatmap_data):
    print(f"{p_bin} pairs\t{row_vals}")

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(heatmap_data, cmap="Blues", aspect='auto')

# X-axis: TE categories
ax.set_xticks(np.arange(len(te_categories)))
ax.set_xticklabels(te_categories)

# Y-axis: # of pairs (with row sum under the label, if desired)
row_labels = []
for p_bin in all_pairs:
    total_in_bin = pairs_counts[p_bin]
    row_labels.append(f"{p_bin} pairs\nn={total_in_bin}")

ax.set_yticks(np.arange(len(all_pairs)))
ax.set_yticklabels(row_labels)

ax.set_xlabel("BLAST (TE) Categories", fontsize=14)
ax.set_ylabel("Number of enhancer pairs in set", fontsize=14)

# Colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Proportion of Sets", rotation=270, labelpad=15)

plt.tight_layout()
plt.show()

##############################################################################
# 6) (Optional) Detailed breakdown
##############################################################################
print("\nDetailed counts per (#pairs bin) × (TE category) (numerator/denominator):")
for p_bin in all_pairs:
    denom = pairs_counts[p_bin]
    for t_cat in te_categories:
        numerator = combo_counts[(p_bin, t_cat)]
        fraction = numerator / denom if denom > 0 else 0.0
        print(f"  {p_bin} pairs & {t_cat}: {numerator}/{denom} = {fraction:.3f}")
    print()  # blank line after each row

