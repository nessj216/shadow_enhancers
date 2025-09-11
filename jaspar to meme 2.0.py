def jaspar_to_meme(jaspar_file_path, meme_file_path, nsites_default=20):
    """
    Convert a JASPAR PWM file to MEME format.

    Parameters:
    - jaspar_file_path: Path to the input JASPAR file.
    - meme_file_path: Path to the output MEME file.
    - nsites_default: Default number of sites to use if not provided in JASPAR.
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

    motif_id, motif_name, matrix = "", "", []

    for line in jaspar_lines:
        if line.startswith('>'):
            if motif_id:  # Save the previous motif before starting a new one
                meme_content.extend(format_meme_motif(motif_id, motif_name, matrix, nsites_default))
                matrix = []

            parts = line.strip().split('\t')
            motif_id = parts[0][1:]  # Remove '>' character
            motif_name = parts[1] if len(parts) > 1 else motif_id
        elif line.strip():
            counts = line.strip().split('[')[-1].split(']')[0].split()
            matrix.append([int(count) for count in counts])

    # Save the last motif
    if motif_id:
        meme_content.extend(format_meme_motif(motif_id, motif_name, matrix, nsites_default))

    with open(meme_file_path, 'w') as meme_file:
        meme_file.write("\n".join(meme_content))


def format_meme_motif(motif_id, motif_name, matrix, nsites):
    """
    Format a single motif's matrix to MEME format.

    Parameters:
    - motif_id: ID of the motif.
    - motif_name: Name of the motif.
    - matrix: Matrix of counts for the motif.
    - nsites: Number of sites for normalization.

    Returns:
    - List of strings representing the motif in MEME format.
    """
    motif_content = [f"MOTIF {motif_id} {motif_name}",
                     f"letter-probability matrix: alength= 4 w= {len(matrix[0])} nsites= {nsites} E= 0"]

    transposed_matrix = list(zip(*matrix))  # Transpose to work with columns
    for column in transposed_matrix:
        total_counts = sum(column)
        probabilities = [count / total_counts for count in column]
        motif_content.append("  " + "\t  ".join(f"{prob:.6f}" for prob in probabilities))

    motif_content.append("")  # Blank line after each motif
    return motif_content


# Example usage
jaspar_file_path = '/Users/jillianness/PycharmProjects/pythonProject1/data-tobias-2020/motifs.jaspar'  # Update this with the path to your JASPAR file
meme_file_path = 'test_output_meme_file.meme'  # Update this with your desired output MEME file path
jaspar_to_meme(jaspar_file_path, meme_file_path)
