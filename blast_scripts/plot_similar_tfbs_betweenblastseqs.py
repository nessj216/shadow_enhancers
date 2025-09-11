import pandas as pd
import matplotlib.pyplot as plt
#checked; it's between all tfbs actoss the seqqeucnes not just blast
# File paths
'''blast_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_v2/CollatedComparisons_f.csv' # Replace with your BLAST file path
fimo_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/TFBS_scan_allshadows/fimo.tsv' # Replace with your FIMO file path

# Read the BLAST CSV file
blast_columns = ["Gene Name", "Comparisons", "qseqid", "sseqid", "pident", "length",
                 "qstart", "qend", "sstart", "send", "evalue", "sseq"]
blast_data = pd.read_csv(blast_file_path, names=blast_columns, sep=",", header=0)

# Read the FIMO TSV file
fimo_columns = ["motif_id", "motif_alt_id", "sequence_name", "start", "stop",
                "strand", "score", "p-value", "q-value", "matched_sequence"]
fimo_data = pd.read_csv(fimo_file_path, names=fimo_columns, sep="\t", header=0)


# Function to parse sequence coordinates
def parse_coordinates(seq):
    try:
        _, coord = seq.split(":")
        start, end = map(int, coord.split("-"))
        return start, end
    except Exception as e:
        print(f"Error parsing sequence: {seq} - {e}")
        return None, None


# Parse start and end coordinates for BLAST sequences
blast_data['qseq_start'], blast_data['qseq_end'] = zip(*blast_data['qseqid'].map(parse_coordinates))
blast_data['sseq_start'], blast_data['sseq_end'] = zip(*blast_data['sseqid'].map(parse_coordinates))

# Calculate lengths
blast_data['qseq_length'] = blast_data['qseq_end'] - blast_data['qseq_start']
blast_data['sseq_length'] = blast_data['sseq_end'] - blast_data['sseq_start']

# Ensure start < end for BLAST hits
blast_data['qstart'], blast_data['qend'] = zip(*blast_data.apply(
    lambda row: (min(row['qstart'], row['qend']), max(row['qstart'], row['qend'])), axis=1))
blast_data['sstart'], blast_data['send'] = zip(*blast_data.apply(
    lambda row: (min(row['sstart'], row['send']), max(row['sstart'], row['send'])), axis=1))

# Ensure start < stop for FIMO TFBS
fimo_data['start'], fimo_data['stop'] = zip(*fimo_data.apply(
    lambda row: (min(row['start'], row['stop']), max(row['start'], row['stop'])), axis=1))

# Group BLAST data by qseqid and sseqid
grouped_blast = blast_data.groupby(['qseqid', 'sseqid'])

# List to store percentage shared TFBS for each pair
shared_percentages = []

for (qseqid, sseqid), group in grouped_blast:
    # Get TFBS sets for qseqid and sseqid
    tfbs_qseq = fimo_data[fimo_data['sequence_name'] == qseqid]['motif_id'].unique()
    tfbs_sseq = fimo_data[fimo_data['sequence_name'] == sseqid]['motif_id'].unique()

    # Convert to sets
    set_qseq = set(tfbs_qseq)
    set_sseq = set(tfbs_sseq)

    # Compute intersection and union
    intersection = set_qseq.intersection(set_sseq)
    union = set_qseq.union(set_sseq)

    if len(union) > 0:
        shared_percentage = (len(intersection) / len(union)) * 100.0
    else:
        # If no TFBS found in either sequence, define shared as 0%
        shared_percentage = 0.0

    shared_percentages.append(shared_percentage)


# Now plot the histogram of shared TFBS percentages
plt.figure(figsize=(10, 6))
plt.hist(shared_percentages, bins=20, color='blue', edgecolor='black', alpha=0.7)
plt.title("Histogram of % Shared TFBS Between qseq and sseq Pairs")
plt.xlabel("% Shared TFBS")
plt.ylabel("Number of Pairs")
plt.tight_layout()
plt.show()

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

# Function to parse sequence coordinates
def parse_coordinates(seq):
    try:
        _, coord = seq.split(":")
        start, end = map(int, coord.split("-"))
        return start, end
    except Exception as e:
        print(f"Error parsing sequence: {seq} - {e}")
        return None, None

# Parse start and end coordinates for BLAST sequences
blast_data['qseq_start'], blast_data['qseq_end'] = zip(*blast_data['qseqid'].map(parse_coordinates))
blast_data['sseq_start'], blast_data['sseq_end'] = zip(*blast_data['sseqid'].map(parse_coordinates))

# Calculate lengths
blast_data['qseq_length'] = blast_data['qseq_end'] - blast_data['qseq_start']
blast_data['sseq_length'] = blast_data['sseq_end'] - blast_data['sseq_start']

# Ensure start < end for BLAST hits
blast_data['qstart'], blast_data['qend'] = zip(*blast_data.apply(
    lambda row: (min(row['qstart'], row['qend']), max(row['qstart'], row['qend'])), axis=1))
blast_data['sstart'], blast_data['send'] = zip(*blast_data.apply(
    lambda row: (min(row['sstart'], row['send']), max(row['sstart'], row['send'])), axis=1))

# Ensure start < stop for FIMO TFBS
fimo_data['start'], fimo_data['stop'] = zip(*fimo_data.apply(
    lambda row: (min(row['start'], row['stop']), max(row['start'], row['stop'])), axis=1))

# Group BLAST data by qseqid and sseqid
grouped_blast = blast_data.groupby(['qseqid', 'sseqid'])

# Lists to store calculated values
shared_percentages_union = []
shared_fraction_qseq = []
shared_fraction_sseq = []

for (qseqid, sseqid), group in grouped_blast:
    # Get TFBS sets for qseqid and sseqid
    tfbs_qseq = fimo_data[fimo_data['sequence_name'] == qseqid]['motif_id'].unique()
    tfbs_sseq = fimo_data[fimo_data['sequence_name'] == sseqid]['motif_id'].unique()

    set_qseq = set(tfbs_qseq)
    set_sseq = set(tfbs_sseq)

    # Compute intersection and union
    intersection = set_qseq.intersection(set_sseq)
    union = set_qseq.union(set_sseq)

    # Intersection/Union percentage
    if len(union) > 0:
        shared_percentage = (len(intersection) / len(union)) * 100.0
    else:
        shared_percentage = 0.0

    # Intersection as fraction of qseq TFBS
    if len(set_qseq) > 0:
        fraction_qseq = len(intersection) / len(set_qseq)
    else:
        fraction_qseq = 0.0

    # Intersection as fraction of sseq TFBS
    if len(set_sseq) > 0:
        fraction_sseq = len(intersection) / len(set_sseq)
    else:
        fraction_sseq = 0.0

    shared_percentages_union.append(shared_percentage)
    shared_fraction_qseq.append(fraction_qseq)
    shared_fraction_sseq.append(fraction_sseq)

# Plot the histograms
plt.figure(figsize=(15, 4))

# 1. Histogram of % Shared TFBS (intersection/union)
plt.subplot(1, 3, 1)
plt.hist(shared_percentages_union, bins=20, color='blue', edgecolor='black', alpha=0.7)
plt.title("% Shared TFBS (Intersection/Union)")
plt.xlabel("% Shared")
plt.ylabel("Number of Pairs")

# 2. Histogram of Intersection/Total qseq TFBS
plt.subplot(1, 3, 2)
plt.hist(shared_fraction_qseq, bins=20, color='green', edgecolor='black', alpha=0.7)
plt.title("Shared Fraction of qseq TFBS (Intersection/qseq)")
plt.xlabel("Fraction Shared")
plt.ylabel("Number of Pairs")

# 3. Histogram of Intersection/Total sseq TFBS
plt.subplot(1, 3, 3)
plt.hist(shared_fraction_sseq, bins=20, color='red', edgecolor='black', alpha=0.7)
plt.title("Shared Fraction of sseq TFBS (Intersection/sseq)")
plt.xlabel("Fraction Shared")
plt.ylabel("Number of Pairs")

plt.tight_layout()
plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

# Function to parse sequence coordinates
def parse_coordinates(seq):
    try:
        _, coord = seq.split(":")
        start, end = map(int, coord.split("-"))
        return start, end
    except Exception as e:
        print(f"Error parsing sequence: {seq} - {e}")
        return None, None

# Parse start and end coordinates for BLAST sequences
blast_data['qseq_start'], blast_data['qseq_end'] = zip(*blast_data['qseqid'].map(parse_coordinates))
blast_data['sseq_start'], blast_data['sseq_end'] = zip(*blast_data['sseqid'].map(parse_coordinates))

# Calculate lengths
blast_data['qseq_length'] = blast_data['qseq_end'] - blast_data['qseq_start']
blast_data['sseq_length'] = blast_data['sseq_end'] - blast_data['sseq_start']

# Ensure start < end for BLAST hits
blast_data['qstart'], blast_data['qend'] = zip(*blast_data.apply(
    lambda row: (min(row['qstart'], row['qend']), max(row['qstart'], row['qend'])), axis=1))
blast_data['sstart'], blast_data['send'] = zip(*blast_data.apply(
    lambda row: (min(row['sstart'], row['send']), max(row['sstart'], row['send'])), axis=1))

# Ensure start < stop for FIMO TFBS
fimo_data['start'], fimo_data['stop'] = zip(*fimo_data.apply(
    lambda row: (min(row['start'], row['stop']), max(row['start'], row['stop'])), axis=1))

# Group BLAST data by qseqid and sseqid (actual pairs)
grouped_blast = blast_data.groupby(['qseqid', 'sseqid'])

# Function to calculate intersection metrics for a pair
def calculate_metrics(qseqid, sseqid):
    tfbs_qseq = fimo_data[fimo_data['sequence_name'] == qseqid]['motif_id'].unique()
    tfbs_sseq = fimo_data[fimo_data['sequence_name'] == sseqid]['motif_id'].unique()

    set_qseq = set(tfbs_qseq)
    set_sseq = set(tfbs_sseq)
    intersection = set_qseq.intersection(set_sseq)
    union = set_qseq.union(set_sseq)

    # Intersection/Union percentage
    if len(union) > 0:
        shared_percentage = (len(intersection) / len(union)) * 100.0
    else:
        shared_percentage = 0.0

    # Fraction of qseq TFBS
    if len(set_qseq) > 0:
        fraction_qseq = len(intersection) / len(set_qseq)
    else:
        fraction_qseq = 0.0

    # Fraction of sseq TFBS
    if len(set_sseq) > 0:
        fraction_sseq = len(intersection) / len(set_sseq)
    else:
        fraction_sseq = 0.0

    return shared_percentage, fraction_qseq, fraction_sseq

# Calculate metrics for actual pairs
shared_percentages_union = []
shared_fraction_qseq = []
shared_fraction_sseq = []

for (qseqid, sseqid), group in grouped_blast:
    sp, fq, fs = calculate_metrics(qseqid, sseqid)
    shared_percentages_union.append(sp)
    shared_fraction_qseq.append(fq)
    shared_fraction_sseq.append(fs)

# Now create random pairs
unique_qseqids = blast_data['qseqid'].unique()
unique_sseqids = blast_data['sseqid'].unique()

num_pairs = len(shared_percentages_union)  # same number of random pairs as actual pairs

# To form random pairs, we randomly choose qseqids and sseqids from the unique sets
# If one set is smaller, consider replacement or ensure both sets are large enough.
random_qseq_choices = np.random.choice(unique_qseqids, size=num_pairs, replace=True)
random_sseq_choices = np.random.choice(unique_sseqids, size=num_pairs, replace=True)

random_shared_percentages_union = []

for qid, sid in zip(random_qseq_choices, random_sseq_choices):
    sp_rand, _, _ = calculate_metrics(qid, sid)
    random_shared_percentages_union.append(sp_rand)

# Plot the random pairs histogram
plt.figure(figsize=(10, 6))
plt.hist(random_shared_percentages_union, bins=20, color='cyan', edgecolor='black', alpha=0.7)
plt.title("Histogram of % Shared TFBS (Intersection/Union) for Randomly Paired Sequences")
plt.xlabel("% Shared TFBS")
plt.ylabel("Number of Random Pairs")
plt.tight_layout()
plt.show()

# Replot the original three histograms for actual pairs
plt.figure(figsize=(15, 4))

# 1. Histogram of % Shared TFBS (intersection/union)
plt.subplot(1, 3, 1)
plt.hist(shared_percentages_union, bins=20, color='blue', edgecolor='black', alpha=0.7)
plt.title("% Shared TFBS (Intersection/Union) - Actual Pairs")
plt.xlabel("% Shared")
plt.ylabel("Number of Pairs")

# 2. Histogram of Intersection/Total qseq TFBS
plt.subplot(1, 3, 2)
plt.hist(shared_fraction_qseq, bins=20, color='green', edgecolor='black', alpha=0.7)
plt.title("Shared Fraction of qseq TFBS (Intersection/qseq) - Actual Pairs")
plt.xlabel("Fraction Shared")
plt.ylabel("Number of Pairs")

# 3. Histogram of Intersection/Total sseq TFBS
plt.subplot(1, 3, 3)
plt.hist(shared_fraction_sseq, bins=20, color='red', edgecolor='black', alpha=0.7)
plt.title("Shared Fraction of sseq TFBS (Intersection/sseq) - Actual Pairs")
plt.xlabel("Fraction Shared")
plt.ylabel("Number of Pairs")

plt.tight_layout()
plt.show()
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# File paths
blast_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_v2/CollatedComparisons_f.csv'  # Replace with your BLAST file path
fimo_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/TFBS_scan_allshadows/fimo.tsv'  # Replace with your FIMO file path

# Read the BLAST CSV file
blast_columns = ["Gene Name", "Comparisons", "qseqid", "sseqid", "pident", "length",
                 "qstart", "qend", "sstart", "send", "evalue", "sseq"]
blast_data = pd.read_csv(blast_file_path, names=blast_columns, sep=",", header=0)

# Read the FIMO TSV file
fimo_columns = ["motif_id", "motif_alt_id", "sequence_name", "start", "stop",
                "strand", "score", "p-value", "q-value", "matched_sequence"]
fimo_data = pd.read_csv(fimo_file_path, names=fimo_columns, sep="\t", header=0)


# Function to parse sequence coordinates
def parse_coordinates(seq):
    try:
        _, coord = seq.split(":")
        start, end = map(int, coord.split("-"))
        return start, end
    except Exception as e:
        print(f"Error parsing sequence: {seq} - {e}")
        return None, None


# Parse start and end coordinates for BLAST sequences
blast_data['qseq_start'], blast_data['qseq_end'] = zip(*blast_data['qseqid'].map(parse_coordinates))
blast_data['sseq_start'], blast_data['sseq_end'] = zip(*blast_data['sseqid'].map(parse_coordinates))

# Calculate lengths
blast_data['qseq_length'] = blast_data['qseq_end'] - blast_data['qseq_start']
blast_data['sseq_length'] = blast_data['sseq_end'] - blast_data['sseq_start']

# Ensure start < end for BLAST hits
blast_data['qstart'], blast_data['qend'] = zip(*blast_data.apply(
    lambda row: (min(row['qstart'], row['qend']), max(row['qstart'], row['qend'])), axis=1))
blast_data['sstart'], blast_data['send'] = zip(*blast_data.apply(
    lambda row: (min(row['sstart'], row['send']), max(row['sstart'], row['send'])), axis=1))

# Ensure start < stop for FIMO TFBS
fimo_data['start'], fimo_data['stop'] = zip(*fimo_data.apply(
    lambda row: (min(row['start'], row['stop']), max(row['start'], row['stop'])), axis=1))

# Group BLAST data by qseqid and sseqid (actual pairs)
grouped_blast = blast_data.groupby(['qseqid', 'sseqid'])


def calculate_metrics(qseqid, sseqid):
    tfbs_qseq = fimo_data[fimo_data['sequence_name'] == qseqid]['motif_id'].unique()
    tfbs_sseq = fimo_data[fimo_data['sequence_name'] == sseqid]['motif_id'].unique()

    set_qseq = set(tfbs_qseq)
    set_sseq = set(tfbs_sseq)
    intersection = set_qseq.intersection(set_sseq)
    union = set_qseq.union(set_sseq)


    # # Debugging: Stop if intersection > 0
    # if len(intersection) > 0:
    #     print(set(tfbs_qseq))
    #     print(set(tfbs_sseq))
    #     print(f"Intersection found for qseqid: {qseqid}, sseqid: {sseqid}")
    #     print(f"Intersection: {intersection}")
    #     print(f"Union: {union}")
    #     input("Press Enter to continue...")  # Pause execution

    # Intersection/Union percentage
    if len(union) > 0:
        shared_percentage = (len(intersection) / len(union)) * 100.0
    else:
        shared_percentage = 0.0

    # Fraction of qseq TFBS
    if len(set_qseq) > 0:
        fraction_qseq = len(intersection) / len(set_qseq)
    else:
        fraction_qseq = 0.0

    # Fraction of sseq TFBS
    if len(set_sseq) > 0:
        fraction_sseq = len(intersection) / len(set_sseq)
    else:
        fraction_sseq = 0.0

    return shared_percentage, fraction_qseq, fraction_sseq


# Calculate metrics for actual pairs
shared_percentages_union_actual = []
shared_fraction_qseq_actual = []
shared_fraction_sseq_actual = []

for (qseqid, sseqid), group in grouped_blast:
    sp, fq, fs = calculate_metrics(qseqid, sseqid)
    shared_percentages_union_actual.append(sp)
    shared_fraction_qseq_actual.append(fq)
    shared_fraction_sseq_actual.append(fs)

# Prepare for randomization
unique_qseqids = blast_data['qseqid'].unique()
unique_sseqids = blast_data['sseqid'].unique()

# num_pairs = len(shared_percentages_union_actual)  # number of actual pairs
# num_randomizations = 100
#
# # Arrays to store random results from all randomizations
# random_shared_union_all = []
# random_fraction_qseq_all = []
# random_fraction_sseq_all = []
#
# for _ in range(num_randomizations):
#     # Random pairing
#     random_qseq_choices = np.random.choice(unique_qseqids, size=num_pairs, replace=True)
#     random_sseq_choices = np.random.choice(unique_sseqids, size=num_pairs, replace=True)
#
#     for qid, sid in zip(random_qseq_choices, random_sseq_choices):
#         sp_rand, fq_rand, fs_rand = calculate_metrics(qid, sid)
#         random_shared_union_all.append(sp_rand)
#         random_fraction_qseq_all.append(fq_rand)
#         random_fraction_sseq_all.append(fs_rand)
#
# # Now plot the histograms for random pairs (aggregated over 100 randomizations)
# plt.figure(figsize=(15, 4))
#
# # Random: % Shared TFBS (Intersection/Union)
# plt.subplot(1, 3, 1)
# plt.hist(random_shared_union_all, bins=20, color='cyan', edgecolor='black', alpha=0.7)
# plt.title("% Shared TFBS (Intersection/Union) - Random (100x)")
# plt.xlabel("% Shared")
# plt.ylabel("Count")
#
# # Random: Intersection/Total qseq TFBS
# plt.subplot(1, 3, 2)
# plt.hist(random_fraction_qseq_all, bins=20, color='purple', edgecolor='black', alpha=0.7)
# plt.title("Shared Fraction qseq (Intersection/qseq) - Random (100x)")
# plt.xlabel("Fraction Shared")
# plt.ylabel("Count")
#
# # Random: Intersection/Total sseq TFBS
# plt.subplot(1, 3, 3)
# plt.hist(random_fraction_sseq_all, bins=20, color='orange', edgecolor='black', alpha=0.7)
# plt.title("Shared Fraction sseq (Intersection/sseq) - Random (100x)")
# plt.xlabel("Fraction Shared")
# plt.ylabel("Count")
#
# plt.tight_layout()
# plt.show()

# Replot the original three histograms for actual pairs, for comparison
plt.figure(figsize=(15, 4))

# Actual: % Shared TFBS (Intersection/Union)
plt.subplot(1, 3, 1)
plt.hist(shared_percentages_union_actual, bins=20, color='blue', edgecolor='black', alpha=0.7)
plt.title("% Shared TFBS (Intersection/Union)")
plt.xlabel("% Shared")
plt.ylabel("Count")

# Actual: Intersection/Total qseq TFBS
plt.subplot(1, 3, 2)
plt.hist(shared_fraction_qseq_actual, bins=20, color='green', edgecolor='black', alpha=0.7)
plt.title("Shared Fraction qseq (Intersection/qseq)")
plt.xlabel("Fraction Shared")
plt.ylabel("Count")

# Actual: Intersection/Total sseq TFBS
plt.subplot(1, 3, 3)
plt.hist(shared_fraction_sseq_actual, bins=20, color='red', edgecolor='black', alpha=0.7)
plt.title("Shared Fraction sseq (Intersection/sseq)")
plt.xlabel("Fraction Shared")
plt.ylabel("Count")

plt.tight_layout()
plt.show()
