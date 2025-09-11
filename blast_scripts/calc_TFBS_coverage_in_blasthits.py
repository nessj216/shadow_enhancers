
##fixed - strand TFBS. should be good
##plots the % coverage of TFBS over each BLAST hit

import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

# Load the BLAST and FIMO data
blast_df = pd.read_csv('/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons_f.csv')
fimo_df = pd.read_csv('/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_SCAN_withinblasthits/fimo_blasthits.tsv', sep='\t')

# Rename columns dynamically based on the structure of the BLAST file
blast_df.columns = [
    "gene_name", "enhancer_id", "qseqid", "sseqid", "pident", "length",
    "qstart", "qend", "sstart", "send", "evalue", "sseq", "Unnamed_12", "Unnamed_13"
][:len(blast_df.columns)]

# Dictionary to store each BLAST hit with a unique key by adding a counter to the `enhancer_id`
blast_lengths = {}
for idx, row in blast_df.iterrows():
    unique_key = f"{row['enhancer_id']}_{idx}"  # Create a unique key with `enhancer_id` and row index
    blast_lengths[unique_key] = row['length']
enhancer_coverage = defaultdict(int)
tfbs_counts = defaultdict(int)

# Iterate through each unique BLAST hit in the modified BLAST file
for unique_key, enhancer_length in blast_lengths.items():
    # Filter FIMO rows where full unique BLAST ID appears as a substring in FIMO's `sequence_name` column
    matching_fimo_rows = fimo_df[fimo_df['sequence_name'].str.contains(unique_key, na=False)]

    # Check if matches were found in FIMO

    '''commented out is for - strand correction'''
   # if not matching_fimo_rows.empty:
        # Adjust coordinates for negative strand hits in FIMO
        # matching_fimo_rows['adjusted_start'] = matching_fimo_rows.apply(
        #     lambda row: row['start'] if row['strand'] == '+' else (enhancer_length - row['stop'] + 1), axis=1
        # )
        # matching_fimo_rows['adjusted_stop'] = matching_fimo_rows.apply(
        #     lambda row: row['stop'] if row['strand'] == '+' else (enhancer_length - row['start'] + 1), axis=1
        # )
        #
        # # Sort TFBS hits by adjusted start position for easier overlap handling
        # tfbs_hits = matching_fimo_rows[['adjusted_start', 'adjusted_stop']].sort_values(by='adjusted_start').to_numpy()
        # Check if matches were found in FIMO
    if not matching_fimo_rows.empty:
            # Use the start and stop coordinates as-is without strand adjustments
        tfbs_hits = matching_fimo_rows[['start', 'stop']].sort_values(by='start').to_numpy()

        # Merge overlapping intervals and calculate total covered length
        merged_intervals = []
        current_start, current_end = tfbs_hits[0]

        for start, end in tfbs_hits[1:]:
            if start <= current_end:  # Overlapping interval
                current_end = max(current_end, end)
            else:
                merged_intervals.append((current_start, current_end))
                current_start, current_end = start, end
        merged_intervals.append((current_start, current_end))

        # Calculate the total covered length from merged intervals
        total_covered_length = sum(end - start + 1 for start, end in merged_intervals)

        # Calculate percentage coverage
        percent_coverage = (total_covered_length / enhancer_length) * 100
        enhancer_coverage[unique_key] = percent_coverage

        # Count the number of TFBS for each BLAST hit and store in tfbs_counts dictionary
        tfbs_counts[unique_key] = len(matching_fimo_rows)
    else:
        print(f"[WARNING] No match found for BLAST hit '{unique_key}' in FIMO file.")

# Convert results to a DataFrame
coverage_df = pd.DataFrame(list(enhancer_coverage.items()), columns=['unique_blast_hit_id', 'percent_coverage'])

# Add columns for BLAST hit lengths and TFBS counts
coverage_df['length'] = coverage_df['unique_blast_hit_id'].apply(lambda x: blast_lengths[x])
coverage_df['tfbs_count'] = coverage_df['unique_blast_hit_id'].apply(lambda x: tfbs_counts[x])

# Calculate the median number of TFBS per BLAST hit
median_tfbs_count = np.median(list(tfbs_counts.values()))
print(f"Median number of TFBS per BLAST hit: {median_tfbs_count}")

# Save results to CSV for further analysis
coverage_df.to_csv('enhancer_tfbs_percent_coverage.csv', index=False)
print("\nCoverage calculation completed. Results saved to 'enhancer_tfbs_percent_coverage.csv'")

# Plotting BLAST length vs. TFBS percent coverage
plt.figure(figsize=(10, 6))
plt.scatter(coverage_df['length'], coverage_df['percent_coverage'], alpha=0.7)
plt.title("BLAST Hit Length vs. TFBS Percent Coverage", fontsize=18)
plt.xlabel("BLAST Hit Length", fontsize=16)
plt.ylabel("TFBS Percent Coverage (%)", fontsize=16)
plt.grid(True)
plt.savefig('blast_length_vs_tfbs_coverage.png')
plt.show()

# Plotting histogram of TFBS counts per BLAST hit
plt.figure(figsize=(10, 6))
plt.hist(coverage_df['tfbs_count'], bins=20, alpha=0.7)
plt.title("Histogram of TFBS Counts per BLAST Hit", fontsize=18)
plt.xlabel("Number of TFBS per BLAST Hit", fontsize=16)
plt.ylabel("Frequency", fontsize=16)
