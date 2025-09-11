import random

def generate_random_dna_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def insert_tf_motifs(sequence, tf_sequences):
    sequence_with_motifs = sequence
    for tf_sequence, num_insertions in tf_sequences.items():
        for _ in range(num_insertions):
            position = random.randint(0, len(sequence_with_motifs) - len(tf_sequence))
            sequence_with_motifs = (
                sequence_with_motifs[:position] + tf_sequence + sequence_with_motifs[position:]
            )
    return sequence_with_motifs

# Define the DNA sequence length
sequence_length = 1000

# Generate a random neutral DNA sequence
neutral_sequence = generate_random_dna_sequence(sequence_length)

# Define TF sequences and the number of times to insert each one
tf_sequences = {
    "TFMOTIF1": 10,   # Insert "TFMOTIF1" 10 times
    "TFMOTIF2": 5,    # Insert "TFMOTIF2" 5 times
    # Add more TF motifs as needed
}

# Insert TF motifs into the sequence
sequence_with_tf = insert_tf_motifs(neutral_sequence, tf_sequences)

# Print or use the resulting sequence as needed
print(sequence_with_tf)
