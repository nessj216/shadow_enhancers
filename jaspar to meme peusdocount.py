def convert_jaspar_to_meme_with_pseudocount(jaspar_file_path, meme_file_path, pseudocount=0.1):
    """
    Convert a JASPAR PWM file to MEME format, including pseudo-counts to avoid zero probabilities.

    Parameters:
    - jaspar_file_path: Path to the input JASPAR file.
    - meme_file_path: Path to the output MEME file.
    - pseudocount: The pseudo-count to add to each count to avoid zero probabilities.
    """
    with open(jaspar_file_path, 'r') as jaspar_file:
        jaspar_lines = jaspar_file.readlines()

    meme_content = [
        "MEME version 4",
        "",
        "ALPHABET= ACGT",
        "",
        "strands: + -",
        "",
        "Background letter frequencies (from uniform background):",
        "A 0.25 C 0.25 G 0.25 T 0.25",
        ""
    ]

    motif_id, motif_name, matrix_lines = "", "", []

    for line in jaspar_lines:
        if line.startswith('>'):
            if motif_id:  # Save the previous motif before starting a new one
                meme_content.append(f"MOTIF {motif_id} {motif_name}")
                meme_content.append(
                    f"letter-probability matrix: alength= 4 w= {len(matrix_lines[0])} nsites= 20 E= 0e+0")
                for position in zip(*matrix_lines):  # Transpose the matrix
                    meme_content.append(" ".join(position) + " ")
                meme_content.append("")
                matrix_lines = []

            parts = line.strip().split('\t')
            motif_id = parts[0][1:]  # Remove '>' character
            motif_name = parts[1] if len(parts) > 1 else motif_id
        elif line.strip():
            counts = line.strip().split('[')[-1].split(']')[0].split()
            adjusted_counts = [int(count) + pseudocount for count in counts]
            total_counts = sum(adjusted_counts)
            probabilities = [f"{count / total_counts:.6f}" for count in adjusted_counts]
            matrix_lines.append(probabilities)

    # Save the last motif
    if motif_id:
        meme_content.append(f"MOTIF {motif_id} {motif_name}")
        meme_content.append(f"letter-probability matrix: alength= 4 w= {len(matrix_lines[0])} nsites= 20 E= 0e+0")
        for position in zip(*matrix_lines):
            meme_content.append(" ".join(position) + " ")
        meme_content.append("")

    with open(meme_file_path, 'w') as meme_file:
        meme_file.write("\n".join(meme_content))


# Example usage
jaspar_file_path = '/Users/jillianness/PycharmProjects/pythonProject1/data-tobias-2020/motifs.jaspar'  # Update this path to your JASPAR file
meme_file_path = 'output_meme_file_pseudo.meme'  # Update this path for your output MEME file
convert_jaspar_to_meme_with_pseudocount(jaspar_file_path, meme_file_path)
