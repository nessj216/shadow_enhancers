'''This code reads in the total shadow file and organizes the gene names into buckets 2 shadows/set, 3-5 or >5
then it finds the # enhancers with hits/ total enhancers in that bucket

this was an extra test after doing #pairs with hits/shadow bucket'''



import collections
import csv

##############################################################################
# 1) Read & bucket SHADOW file
##############################################################################
shadow_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/011925_all_shadowsets_DM6.bed"

# Step A: Count how many times each set appears
shadow_counts = collections.Counter()
with open(shadow_file, 'r') as sf:
    for line in sf:
        if line.strip():
            cols = line.strip().split('\t')
            set_name = cols[3]  # Assume set name is in column 4 (index 3)
            shadow_counts[set_name] += 1

# Step B: Assign each set to a shadow bucket
shadow_buckets = {}
for set_name, count in shadow_counts.items():
    if count == 2:
        shadow_buckets[set_name] = "2 shadows"
    elif 3 <= count <= 4:
        shadow_buckets[set_name] = "3-4 shadows"
    elif count >= 5:
        shadow_buckets[set_name] = ">=5 shadows"

# Step C: Count how many sets fall into each shadow category
shadow_cat_counts = collections.Counter(shadow_buckets.values())

# Step D: Count total *enhancers/lines* for each bucket
bucket_enhancer_counts = collections.Counter()
for set_name, bucket_name in shadow_buckets.items():
    bucket_enhancer_counts[bucket_name] += shadow_counts[set_name]

# Print shadow set results
print("Bucket name -> (Number of sets, Total enhancers in bucket)")
for bucket_name in ["2 shadows", "3-4 shadows", ">=5 shadows"]:
    print(f"{bucket_name} -> ({shadow_cat_counts[bucket_name]} sets, {bucket_enhancer_counts[bucket_name]} enhancers)")

##############################################################################
# 2) Read & Process HIT file, Extract Unique IDs, and Count per Bucket
##############################################################################
hit_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/HEATMAP_stuff/BLAStHITS_deduplicated_lines.csv"

# Initialize bucket storage for unique IDs
shadow_ids = {"2 shadows": set(), "3-4 shadows": set(), ">=5 shadows": set()}

# Read hit file and classify each row into the correct bucket
with open(hit_file, 'r') as hf:
    reader = csv.reader(hf)
    next(reader)  # Skip header if present

    for row in reader:
        if len(row) >= 2:
            set_name, col2_value = row[0], row[1]  # First column is set_name, second is filename

            # Extract ID1 and ID2 from the filename (remove ".txt" and split by "_")
            if col2_value.endswith(".txt"):
                col2_value = col2_value[:-4]  # Remove ".txt"
            id_parts = col2_value.split("_")  # Split into ID1 and ID2

            # Assign IDs to the correct shadow bucket
            if set_name in shadow_buckets:
                bucket = shadow_buckets[set_name]
                shadow_ids[bucket].update(id_parts)  # Add both IDs to the bucket's set (avoiding duplicates)

# Print total unique IDs per bucket
print("\nTotal Unique IDs in Each Shadow Bucket:")
for bucket, unique_ids in shadow_ids.items():
    print(f"{bucket}: {len(unique_ids)} unique IDs")

##############################################################################
# 3) Compute Portion of Hit Enhancers in Each Bucket
##############################################################################
portion_hit_enhancers = {}

for bucket in ["2 shadows", "3-4 shadows", ">=5 shadows"]:
    total_ids = len(shadow_ids[bucket])
    total_enhancers = bucket_enhancer_counts[bucket]

    # Avoid division by zero
    if total_enhancers > 0:
        portion_hit_enhancers[bucket] = total_ids / total_enhancers
    else:
        portion_hit_enhancers[bucket] = 0  # Assign 0 if no enhancers

# Print portion results
print("\nPortion of Hit Enhancers in Each Bucket:")
for bucket, portion in portion_hit_enhancers.items():
    print(f"{bucket}: {portion:.3f}")  # Format to 3 decimal places

import itertools
import scipy.stats as stats
import statsmodels.stats.multitest as smm
import collections
import csv

# Define observed proportions and total counts
buckets = ["2 shadows", "3-4 shadows", ">=5 shadows"]

# Given data
bucket_enhancer_counts = {
    "2 shadows": 434,
    "3-4 shadows": 520,
    ">=5 shadows": 296
}

shadow_ids = {
    "2 shadows": set(range(24)),  # Simulating 24 unique IDs
    "3-4 shadows": set(range(53)),  # Simulating 53 unique IDs
    ">=5 shadows": set(range(71))  # Simulating 71 unique IDs
}



import matplotlib.pyplot as plt

# Data for plotting
buckets = ["2 shadows/set", "3-4 shadows/set", "≥5 shadows/set"]
proportions = [24 / 434, 53 / 520, 71 / 296]



import matplotlib.pyplot as plt

# Data for plotting
buckets = ["2 shadows/set", "3-4 shadows/set", "≥5 shadows/set"]
proportions = [24 / 434, 53 / 520, 71 / 296]
total_enhancers = [434, 520, 296]  # Total enhancers in each bucket

# Define custom colors
colors = ["lightblue", "royalblue", "navy"]

# Creating the bar plot
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(buckets, proportions, color=colors, alpha=0.75)

# Adding labels and title
#ax.set_xlabel("Shadow Bucket")
ax.set_ylabel("Proportion of Hit Enhancers")

# Significance lines and p-values
comparisons = [(0, 1), (0, 2), (1, 2)]
y_max = max(proportions)  # Get the highest bar for placement
y_spacing = 0.02  # Space between lines

# Bonferroni-corrected p-values (replace with actual values)
p_values_corrected = ['.02', 'p<1e-6', 'p<1e-6']

for idx, (i, j) in enumerate(comparisons):
    x1, x2 = i, j
    y = y_max + (idx + 1) * y_spacing - 0.02  # Place the line slightly above the highest bar

    # Draw horizontal line
    ax.plot([x1, x2], [y, y], color="black", lw=1.5)

    # Add p-value
    ax.text((x1 + x2) / 2, y + 0.005, f"p = {p_values_corrected[idx]}", ha="center", fontsize=10)

# Add total enhancer count below each x label
for i, (bucket, enh_count) in enumerate(zip(buckets, total_enhancers)):
    ax.text(i, -0.10, f"(n={enh_count} enhancers)", ha="center", va="top", fontsize=13, transform=ax.get_xaxis_transform())

# Adjust y-limits to accommodate annotations
ax.set_ylim(0, y_max + 3 * y_spacing)
ax.set_ylabel("Proportion of Hit Enhancers", fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=15)

# Display the plot
plt.show()
