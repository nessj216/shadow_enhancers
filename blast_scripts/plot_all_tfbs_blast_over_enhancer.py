import matplotlib.pyplot as plt
import pandas as pd

# Load the uploaded file to examine the headers
#file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_v2/CollatedComparisons_f.csv'

import pandas as pd
import matplotlib.pyplot as plt

# File paths
blast_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_v2/CollatedComparisons_f.csv' # Replace with the BLAST file path
fimo_file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/TFBS_scan_allshadows/fimo.tsv' # Replace with the FIMO file path

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

# Plot grouped data with TFBS
for (qseqid, sseqid), group in grouped_blast:
    plt.figure(figsize=(10, 6))

    # Get unique gene names
    gene_names = " | ".join(group['Gene Name'].unique())

    # Plot qseqid
    qseq_length = group['qseq_length'].iloc[0]
    plt.plot([0, qseq_length], [1, 1], color='black', label='qseqid Sequence', alpha=0.5)
    for _, row in group.iterrows():
        plt.plot([row['qstart'], row['qend']], [1, 1], color='red', linewidth=2, label='qseqid BLAST Hit')

    # Add TFBS for qseqid
    tfbs_qseq = fimo_data[fimo_data['sequence_name'] == qseqid]
    for _, row in tfbs_qseq.iterrows():
        plt.plot([row['start'], row['stop']], [1.1, 1.1], color='purple', linewidth=2, label='qseqid TFBS')

    # Plot sseqid
    sseq_length = group['sseq_length'].iloc[0]
    plt.plot([0, sseq_length], [0, 0], color='black', label='sseqid Sequence', alpha=0.5)
    for _, row in group.iterrows():
        plt.plot([row['sstart'], row['send']], [0, 0], color='orange', linewidth=2, label='sseqid BLAST Hit')

    # Add TFBS for sseqid
    tfbs_sseq = fimo_data[fimo_data['sequence_name'] == sseqid]
    for _, row in tfbs_sseq.iterrows():
        plt.plot([row['start'], row['stop']], [-0.1, -0.1], color='purple', linewidth=2, label='sseqid TFBS')

    # Labels and title
    plt.title(f"BLAST Hits and TFBS for qseqid: {qseqid} and sseqid: {sseqid}\nGene Names: {gene_names}")
    plt.xlabel("Position")
    plt.yticks([-0.1, 0, 1, 1.1], ['TFBS (sseqid)', 'sseqid', 'qseqid', 'TFBS (qseqid)'])
    #plt.legend()
    # Adjust spacing between the plots
    plt.subplots_adjust(top=0.3, bottom=0.1, left=0.2, right=0.9, hspace=0.8)

    #plt.tight_layout()
    plt.show()
