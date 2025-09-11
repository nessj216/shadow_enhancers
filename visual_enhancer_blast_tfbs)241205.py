import pandas as pd
import matplotlib.pyplot as plt

# Load the BLAST and FIMO data
blast_df = pd.read_csv(
    '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_newID/Final_Updated_Collated_Comparisons_with_Indexed.csv')  # Update file paths as needed
fimo_df = pd.read_csv(
    '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_newID/fimo10-3.tsv',
    sep='\t')

# Rename columns in BLAST file
blast_df.columns = [
                       "gene_name", "enhancer_id", "qseqid", "sseqid", "pident", "length",
                       "qstart", "qend", "sstart", "send", "evalue", "sseq"
                   ][:len(blast_df.columns)]

# Extract enhancer_id prefix before ".txt" for grouping
blast_df['enhancer_prefix'] = blast_df['enhancer_id'].str.split('.txt').str[0]


# Helper function to calculate sequence length
def calculate_length(coord_string):
    try:
        _, positions = coord_string.split(':')
        start, end = map(int, positions.split('-'))
        return end - start
    except:
        return 0


# Add calculated lengths for qseqid and sseqid (for enhancers)
blast_df['qseq_length'] = blast_df['qseqid'].apply(calculate_length)
blast_df['sseq_length'] = blast_df['sseqid'].apply(calculate_length)

# Correct BLAST coordinates if start > stop
blast_df['qstart'], blast_df['qend'] = (
    blast_df[['qstart', 'qend']].min(axis=1),
    blast_df[['qstart', 'qend']].max(axis=1),
)
blast_df['sstart'], blast_df['send'] = (
    blast_df[['sstart', 'send']].min(axis=1),
    blast_df[['sstart', 'send']].max(axis=1),
)


# Function to merge overlapping intervals
def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    for start, end in sorted_intervals:
        if not merged or merged[-1][1] < start:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


# Plot each group of BLAST hits and TFBS
for enhancer_prefix, group in blast_df.groupby('enhancer_prefix'):
    plt.figure(figsize=(14, 6))

    for (enhancer_id, qseqid, sseqid), sub_group in group.groupby(['enhancer_id', 'qseqid', 'sseqid']):
        qseq_length = sub_group['qseq_length'].iloc[0]
        sseq_length = sub_group['sseq_length'].iloc[0]

        # Merge overlapping BLAST hits for qseqid and sseqid
        qseq_intervals = merge_intervals(sub_group[['qstart', 'qend']].values.tolist())
        sseq_intervals = merge_intervals(sub_group[['sstart', 'send']].values.tolist())

        plt.hlines(y=2, xmin=0, xmax=qseq_length, color='gray', linestyles='dotted', label='qseqid Region')
        plt.hlines(y=1, xmin=0, xmax=sseq_length, color='gray', linestyles='dotted', label='sseqid Region')

        # Plot merged BLAST hits for qseqid
        for qstart, qend in qseq_intervals:
            plt.plot(
                [qstart, qend],
                [2, 2],
                color='green',
                linewidth=2,
                label='Merged BLAST hit (qseqid)' if qseq_intervals.index([qstart, qend]) == 0 else ""
            )

        # Plot merged BLAST hits for sseqid
        for sstart, send in sseq_intervals:
            plt.plot(
                [sstart, send],
                [1, 1],
                color='orange',
                linewidth=2,
                label='Merged BLAST hit (sseqid)' if sseq_intervals.index([sstart, send]) == 0 else ""
            )

        # Add TFBS to plot, adjusting coordinates for both sequences
        for _, row in sub_group.iterrows():
            matching_fimo = fimo_df[fimo_df['sequence_name'].str.contains(row['enhancer_id'], na=False)]
            for _, fimo_row in matching_fimo.iterrows():
                # Adjust TFBS coordinates for qseqid
                tfbs_start_q = fimo_row['start'] + row['qstart'] - 1
                tfbs_stop_q = fimo_row['stop'] + row['qstart'] - 1
                # Adjust TFBS coordinates for sseqid
                tfbs_start_s = fimo_row['start'] + row['sstart'] - 1
                tfbs_stop_s = fimo_row['stop'] + row['sstart'] - 1

                # Plot TFBS on qseqid
                plt.plot(
                    [tfbs_start_q, tfbs_stop_q],
                    [2, 2],
                    color='blue' if fimo_row['strand'] == '+' else 'red',
                    linewidth=1.5,
                    linestyle='--',
                    label='TFBS (+)' if fimo_row['strand'] == '+' else 'TFBS (-)'
                )
                # Plot TFBS on sseqid
                plt.plot(
                    [tfbs_start_s, tfbs_stop_s],
                    [1, 1],
                    color='blue' if fimo_row['strand'] == '+' else 'red',
                    linewidth=1.5,
                    linestyle='--'
                )

    plt.title(f"Merged BLAST Hits and TFBS for {enhancer_prefix}", fontsize=14)
    plt.xlabel("Position", fontsize=12)
    plt.yticks([1, 2], labels=['sseqid', 'qseqid'])
    plt.legend(loc='upper right')
    plt.grid(axis='x')
    plt.tight_layout()
    plt.savefig(f"merged_blast_tfbs_{enhancer_prefix}.png")
    plt.show()