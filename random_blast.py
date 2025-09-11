import os
import random
from tqdm import tqdm

# Set the path to the directory containing all sequence files
sequences_dir = "/Users/jillianness/Desktop/Mouse_SE_birth_analysis/random_mouse/comparison"
sequences = [f for f in os.listdir(sequences_dir) if
             os.path.isfile(os.path.join(sequences_dir, f)) and f != '.DS_Store']

# Create Comparisons_0001 directory
comparisons_folder = os.path.join(sequences_dir, 'Comparisons_0001')
os.makedirs(comparisons_folder, exist_ok=True)

# Perform BLAST for approximately 10,000 iterations
for _ in tqdm(range(10000)):
    # Randomly select two different sequence files
    text1, text2 = random.sample(sequences, 2)
    file1 = os.path.join(sequences_dir, text1)
    file2 = os.path.join(sequences_dir, text2)

    # Remove file extensions for output file naming
    text1_base = os.path.splitext(text1)[0]
    text2_base = os.path.splitext(text2)[0]

    # Run BLASTn command
    comparison = os.popen(
        f"blastn -evalue '1' -word_size '7' -gapopen '5' -gapextend '2' -reward '2' -penalty '-3' -dust 'yes' "
        f"-query {file1} -subject {file2} -outfmt '6 qseqid sseqid pident length evalue sseq'"
    )

    output = comparison.read()

    # Write the BLAST output to a file
    output_file = os.path.join(comparisons_folder, f"{text1_base}_{text2_base}")
    with open(output_file, 'w') as f:
        f.write(output)

# Note: The blastn command assumes that BLAST+ is installed and accessible from the command line.
