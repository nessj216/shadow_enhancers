import pandas as pd
import matplotlib.pyplot as plt

# Load the BLAST and FIMO data
blast_df = pd.read_csv(
    '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_newID/Final_Updated_Collated_Comparisons_with_Indexed.csv'
)  # Update file paths as needed
fimo_df = pd.read_csv(
    '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_newID/fimo.tsv',
    sep='\t'
)

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

# Track original orientation and correct BLAST coordinates if start > stop
blast_df['q_flipped'] = blast_df['qstart'] > blast_df['qend']
blast_df['s_flipped'] = blast_df['sstart'] > blast_df['send']

blast_df['qstart'], blast_df['qend'] = (
    blast_df[['qstart', 'qend']].min(axis=1),
    blast_df[['qstart', 'qend']].max(axis=1),
)
blast_df['sstart'], blast_df['send'] = (
    blast_df[['sstart', 'send']].min(axis=1),
    blast_df[['sstart', 'send']].max(axis=1),
)

# Group by gene_name, qseqid, and sseqid for plotting
for (gene_name, qseqid, sseqid), group in blast_df.groupby(['gene_name', 'qseqid', 'sseqid']):
    plt.figure(figsize=(14, 6))

    print(f"Plotting group for:")
    print(f"  Gene Name: {gene_name}")
    print(f"  qseqid: {qseqid}")
    print(f"  sseqid: {sseqid}")
    print(f"  Enhancer IDs in this group: {group['enhancer_id'].unique()}")

    qseq_length = group['qseq_length'].iloc[0]
    sseq_length = group['sseq_length'].iloc[0]

    plt.hlines(y=2, xmin=0, xmax=qseq_length, color='gray', linestyles='dotted', label='qseqid Region')
    plt.hlines(y=1, xmin=0, xmax=sseq_length, color='gray', linestyles='dotted', label='sseqid Region')

    # Plot raw BLAST hits for qseqid and sseqid
    for _, row in group.iterrows():
        plt.plot([row['qstart'], row['qend']], [2, 2], color='green', linewidth=2, label='BLAST hit (qseqid)')
        plt.plot([row['sstart'], row['send']], [1, 1], color='orange', linewidth=2, label='BLAST hit (sseqid)')

        # Add TFBS to plot, adjusting coordinates for both sequences
        matching_fimo = fimo_df[fimo_df['sequence_name'].str.contains(row['enhancer_id'], na=False)]
        for _, fimo_row in matching_fimo.iterrows():
            # Adjust TFBS coordinates for qseqid
            if row['q_flipped']:
                tfbs_start_q = row['qend'] - (fimo_row['start'] - 1)
                tfbs_stop_q = row['qend'] - (fimo_row['stop'] - 1)
            else:
                tfbs_start_q = fimo_row['start'] + row['qstart'] - 1
                tfbs_stop_q = fimo_row['stop'] + row['qstart'] - 1

            # Adjust TFBS coordinates for sseqid
            if row['s_flipped']:
                tfbs_start_s = row['send'] - (fimo_row['start'] - 1)
                tfbs_stop_s = row['send'] - (fimo_row['stop'] - 1)
            else:
                tfbs_start_s = fimo_row['start'] + row['sstart'] - 1
                tfbs_stop_s = fimo_row['stop'] + row['sstart'] - 1

            print(f"TFBS start/stop on qseqid: {tfbs_start_q}, {tfbs_stop_q}", row['enhancer_id'])
            print(f"TFBS start/stop on sseqid: {tfbs_start_s}, {tfbs_stop_s}", row['enhancer_id'])

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

    plt.title(f"BLAST Hits and TFBS for {gene_name} | {qseqid} | {sseqid}", fontsize=14)
    plt.xlabel("Position", fontsize=12)
    plt.yticks([1, 2], labels=['sseqid', 'qseqid'])
    plt.legend(loc='upper right')
    plt.grid(axis='x')
    plt.tight_layout()
    plt.savefig(f"blast_tfbs_{gene_name}_{qseqid}_{sseqid}.png")
    plt.show()

