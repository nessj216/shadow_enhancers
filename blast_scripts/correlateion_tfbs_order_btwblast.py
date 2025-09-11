
#only plots the similarity between blast tfbs , not entire sequences

import pandas as pd
import numpy as np

# Parameters
window_size = 50  # window size in bp for qseq
step_size = 20    # step size in bp for qseq sliding
import pandas as pd
import matplotlib.pyplot as plt

# File paths
blast_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_v2/CollatedComparisons_f.csv' # Replace with your BLAST file path
fimo_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/TFBS_scan_allshadows/fimo.tsv' # Replace with your FIMO file path

# Read the BLAST CSV file
blast_columns = ["Gene Name", "Comparisons", "qseqid", "sseqid", "pident", "length",
                 "qstart", "qend", "sstart", "send", "evalue", "sseq"]
blast_data = pd.read_csv(blast_file_path, names=blast_columns, sep=",", header=0)

# Read the FIMO TSV file
fimo_columns = ["motif_id", "motif_alt_id", "sequence_name", "start", "stop",
                "strand", "score", "p-value", "q-value", "matched_sequence"]
fimo_data = pd.read_csv(fimo_file_path, names=fimo_columns, sep="\t", header=0)


# Assume you have blast_data and fimo_data from before
# blast_data: contains qseqid, sseqid, and alignment coordinates (qstart, qend, sstart, send)
# fimo_data: contains motif_id and coordinates for each TFBS in each sequence

def parse_coordinates(seq):
    try:
        _, coord = seq.split(":")
        start, end = map(int, coord.split("-"))
        return start, end
    except Exception as e:
        print(f"Error parsing sequence: {seq} - {e}")
        return None, None

# Make sure blast_data has qseq_start, qseq_end, sseq_start, sseq_end
blast_data['qseq_start'], blast_data['qseq_end'] = zip(*blast_data['qseqid'].map(parse_coordinates))
blast_data['sseq_start'], blast_data['sseq_end'] = zip(*blast_data['sseqid'].map(parse_coordinates))

# Sort coordinates to ensure start < end
blast_data['qstart'], blast_data['qend'] = zip(*blast_data.apply(
    lambda row: (min(row['qstart'], row['qend']), max(row['qstart'], row['qend'])), axis=1))
blast_data['sstart'], blast_data['send'] = zip(*blast_data.apply(
    lambda row: (min(row['sstart'], row['send']), max(row['sstart'], row['send'])), axis=1))

# Group by pairs
grouped_blast = blast_data.groupby(['qseqid', 'sseqid'])

# Function to get TFBS in a given sequence and coordinate range
def get_tfbs_in_range(sequence_name, start_coord, end_coord):
    subset = fimo_data[(fimo_data['sequence_name'] == sequence_name) &
                       (fimo_data['start'] >= start_coord) &
                       (fimo_data['stop'] <= end_coord)]
    return subset['motif_id'].tolist(), subset[['motif_id', 'start', 'stop']]

# Function to map qseq window coordinates to sseq coordinates based on BLAST alignment
def map_coordinates(q_coord, qstart, qend, sstart, send):
    # linear mapping:
    # If q_coord is relative to qstart, normalized coordinate = (q_coord - qstart)/(qend - qstart)
    # sseq_coord = sstart + normalized * (send - sstart)
    if qend == qstart:
        # degenerate case, avoid division by zero
        return sstart
    fraction = (q_coord - qstart) / (qend - qstart)
    return int(round(sstart + fraction * (send - sstart)))

# Prepare a list to store window-level results
results = []

for (qseqid, sseqid), group in grouped_blast:
    # For simplicity, take the first row to get alignment info (assuming one alignment per pair)
    row = group.iloc[0]
    qstart, qend = row['qstart'], row['qend']
    sstart, send = row['sstart'], row['send']

    qseq_region_start = row['qseq_start']
    qseq_region_end = row['qseq_end']

    # We'll slide windows along the qseq coordinates
    # Assume we only analyze the aligned region within qstart/qend to be relevant
    # Adjust as needed.
    q_aligned_length = qend - qstart

    # Start from qstart and slide in steps until qend-window_size
    for q_window_start in range(qstart, qend - window_size + 1, step_size):
        q_window_end = q_window_start + window_size

        # Map q_window_start and q_window_end to sseq coordinates
        s_window_start = map_coordinates(q_window_start, qstart, qend, sstart, send)
        s_window_end = map_coordinates(q_window_end, qstart, qend, sstart, send)

        # Get TFBS in qseq window
        q_tfbs_list, q_tfbs_details = get_tfbs_in_range(qseqid, q_window_start + row['qseq_start'] - 1,
                                                        q_window_end + row['qseq_start'] - 1)
        # Note: Adjust indexing by qseq_start if necessary depending on how coordinates are defined.
        # The original code uses absolute coordinates from qseqid. If q_window_* are relative to qstart,
        # and qstart refers to qseqid's actual coordinates, you might need to carefully check indexing.
        # Here we assume qstart/qend are absolute to qseqid coordinates. If not, adjust accordingly.

        # Get TFBS in sseq window
        s_tfbs_list, s_tfbs_details = get_tfbs_in_range(sseqid, s_window_start + row['sseq_start'] - 1,
                                                        s_window_end + row['sseq_start'] - 1)

        # Measure intersection or order
        q_tfbs_set = set(q_tfbs_list)
        s_tfbs_set = set(s_tfbs_list)
        intersection = q_tfbs_set.intersection(s_tfbs_set)

        # One way to measure "grammar" would be just intersection count:
        # More sophisticated measures could consider the order:
        # For example, sorting TFBS by start position and measuring correlation of their ranks.
        # Here, we just record intersection count:
        intersection_count = len(intersection)

        # Order-based measure (example):
        # Sort q_tfbs by start, create a list of motif_ids in order
        q_tfbs_order = q_tfbs_details.sort_values('start')['motif_id'].tolist()
        s_tfbs_order = s_tfbs_details.sort_values('start')['motif_id'].tolist()

        # If you wanted to measure how similar the order is, you could, for example,
        # count how many motifs in intersection appear in the same relative order:
        # A simple metric: how many intersection TFBS appear in the same order in both?
        # This is a simplistic approach:
        common_order = [m for m in q_tfbs_order if m in intersection]
        s_common_order = [m for m in s_tfbs_order if m in intersection]

        # Measure order similarity (for example, fraction of intersection motifs that appear in the same sorted order)
        # A simple approach: check if the order of intersection motifs in s_common_order matches q_common_order
        # We'll count how many pairs of motifs are in the same relative order.
        # For a more advanced measure, consider metrics like Kendall’s tau or Spearman correlation of indices.
        if len(common_order) > 1:
            # Map motifs to their indices in s_common_order
            s_indices = [s_common_order.index(m) for m in common_order]
            # If the order is the same, s_indices should be increasing.
            # Compute a rank correlation:
            order_cor = pd.Series(range(len(s_indices))).corr(pd.Series(s_indices), method='spearman')
        else:
            order_cor = np.nan  # Not enough TFBS to compute order correlation

        results.append({
            'qseqid': qseqid,
            'sseqid': sseqid,
            'q_window_start': q_window_start,
            'q_window_end': q_window_end,
            's_window_start': s_window_start,
            's_window_end': s_window_end,
            'intersection_count': intersection_count,
            'order_correlation': order_cor,
            'q_tfbs_count': len(q_tfbs_list),
            's_tfbs_count': len(s_tfbs_list),
            'intersection_size': len(intersection)
        })

# Convert results into a DataFrame
window_results = pd.DataFrame(results)

# Now you can analyze `window_results`:
# For example, plot intersection_count vs q_window_start, or histogram of order_correlation
# Example plotting (pseudocode):
#
# import matplotlib.pyplot as plt
plt.figure()
plt.hist(window_results['intersection_count'], bins=20)
plt.title("Distribution of Intersection Counts in Sliding Windows")
plt.xlabel("Intersection Count")
plt.ylabel("Frequency")
plt.show()

#Similarly, for order_correlation:
plt.figure()
plt.hist(window_results['order_correlation'].dropna(), bins=20)
plt.title("Distribution of TFBS Order Correlations in Sliding Windows")
plt.xlabel("Spearman Correlation")
plt.ylabel("Frequency")
plt.show()

plt.figure(figsize=(6,4))
plt.scatter(window_results['intersection_size'], window_results['order_correlation'],
            color='teal', alpha=0.7)
plt.title("TFBS Order Correlation vs. Intersection Size")
plt.xlabel("Intersection Size (Number of Shared TFBS)")
plt.ylabel("Spearman Correlation of TFBS Order")
plt.tight_layout()
plt.show()
