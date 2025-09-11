
import collections
import csv

# Define file paths
shadow_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/011925_all_shadowsets_DM6.bed"
hit_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/HEATMAP_stuff/BLAStHITS_deduplicated_lines.csv"

# Count occurrences of each set name
shadow_counts = collections.Counter()

# Read the shadow file and populate the counter
with open(shadow_file, 'r') as sf:
    for line in sf:
        if line.strip():
            cols = line.strip().split('\t')
            set_name = cols[3]  # Assuming set name is in column 4 (index 3)
            shadow_counts[set_name] += 1

# Assign each set to a shadow bucket
shadow_buckets = {}
shadow_pairs = {}

for set_name, count in shadow_counts.items():
    if count == 2:
        bucket = "2 shadows"
    elif 3 <= count <= 4:
        bucket = "3-4 shadows"
    elif count >= 5:
        bucket = ">=5 shadows"
    else:
        continue  # Skip sets with 1 or 0 shadows

    shadow_buckets[set_name] = bucket

    # Calculate total number of possible pairs
    total_pairs = count * (count - 1) // 2
    shadow_pairs[set_name] = total_pairs

# Aggregate total pairs by shadow bucket
bucket_pair_totals = collections.defaultdict(int)
for set_name, bucket in shadow_buckets.items():
    bucket_pair_totals[bucket] += shadow_pairs[set_name]

# Print the total pairs for each shadow bucket
for bucket, total_pairs in bucket_pair_totals.items():
    print(f"{bucket}: {total_pairs} total pairs")

# Process the hit file and sort col2 into shadow buckets
sorted_hits = collections.defaultdict(list)

with open(hit_file, 'r') as hf:
    reader = csv.reader(hf)
    next(reader)  # Skip header if present
    for row in reader:
        if len(row) < 2:
            continue  # Skip malformed lines
        set_name, col2_value = row[0], row[1]  # col1 = set_name, col2 = value

        # Find the corresponding shadow bucket
        bucket = shadow_buckets.get(set_name, "Uncategorized")
        sorted_hits[bucket].append(col2_value)

# Print sorted results
for bucket, values in sorted_hits.items():
    print(f"{bucket}: {len(values)} entries")
