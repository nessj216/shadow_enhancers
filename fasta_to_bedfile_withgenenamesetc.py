import os
import pandas as pd
from Bio import SeqIO

# Define your top-level directory where all subdirectories are
base_dir = '/Users/jillianness/Desktop/mouse_analysis_031925/Mouse_enhancers/Fasta_files/Shadow_fasta_files/FlilesOutput'  # <-- CHANGE THIS

# Storage for all rows
rows = []

# Walk through all subdirectories
for subdir, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith('.fa') or file.endswith('.fasta') or file.endswith('.txt'):  # Accept common fasta extensions
            file_path = os.path.join(subdir, file)
            subdir_name = os.path.basename(subdir)

            # Parse fasta
            for record in SeqIO.parse(file_path, 'fasta'):
                header = record.description  # Full header string
                if header.startswith('>'):
                    header = header[1:]

                # Parse the header: chr:start-end
                try:
                    chrom_part, pos_part = header.split(':')
                    start, end = pos_part.split('-')
                    rows.append([chrom_part, int(start), int(end), subdir_name, file])
                except ValueError:
                    print(f"Warning: couldn't parse header '{header}' in file {file_path}")

# Create DataFrame
bed_df = pd.DataFrame(rows, columns=['Chromosome', 'Start', 'End', 'Filename', 'Subdirectory'])

# Save to a tab-separated .bed-like file
output_path = '/Users/jillianness/Desktop/mouse_analysis_031925/Mouse_enhancers/Fasta_files/Shadow_fasta_files/enhancercombined.bed'  # <-- CHANGE THIS
bed_df.to_csv(output_path, sep='\t', header=False, index=False)

print(f"Done! BED file saved to: {output_path}")
