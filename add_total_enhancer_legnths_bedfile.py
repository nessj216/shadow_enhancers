import os

# Path to the BED file
bed_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/single_enhancer_cooption/singleenhancer_50bpfilerTEhits.bed"

total_length = 0

with open(bed_file, 'r') as file:
    for line in file:
        # Skip empty lines or malformed lines
        if not line.strip():
            continue

        columns = line.strip().split()

        # Ensure at least 3 columns: chromosome, start, end
        if len(columns) >= 3:
            start = int(columns[1])
            end = int(columns[2])
            length = end - start
            total_length += length

print(f"Total sum of lengths: {total_length}")
