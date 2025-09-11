import matplotlib.pyplot as plt

# Input and output file paths
input_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_splitting/output_splitingmergedTE_enhancer.bed"
output_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/60_overlap_hits.bed"

# Function to calculate overlap length
def calculate_overlap_length(te_start, te_end, enhancer_start, enhancer_end):
    overlap_start = max(te_start, enhancer_start)
    overlap_end = min(te_end, enhancer_end)
    overlap_length = max(0, overlap_end - overlap_start)
    return overlap_length

# Store overlaps by enhancer
enhancer_overlaps = {}

# Read input file and calculate overlaps
with open(input_file, 'r') as infile:
    for line in infile:
        if line.strip():  # Skip empty lines
            columns = line.strip().split('\t')
            # Extract coordinates
            chr_name = columns[0]
            enhancer_start = int(columns[1])
            enhancer_end = int(columns[2])
            te_start = int(columns[6])
            te_end = int(columns[7])

            # Define enhancer key based on chr, start, end
            enhancer_key = (chr_name, enhancer_start, enhancer_end)
            enhancer_length = enhancer_end - enhancer_start

            # Calculate overlap length
            overlap_length = calculate_overlap_length(te_start, te_end, enhancer_start, enhancer_end)

            # Store overlap for the given enhancer
            if enhancer_key not in enhancer_overlaps:
                enhancer_overlaps[enhancer_key] = {'total_overlap': 0, 'enhancer_length': enhancer_length, 'line': line.strip()}

            enhancer_overlaps[enhancer_key]['total_overlap'] += overlap_length

# Calculate overlap percentages and write hits to output file
with open(output_file, 'w') as outfile:
    for data in enhancer_overlaps.values():
        enhancer_length = data['enhancer_length']
        total_overlap = data['total_overlap']
        overlap_percentage = (total_overlap / enhancer_length) * 100 if enhancer_length > 0 else 0

        # Write to output file if overlap percentage is greater than 60
        if overlap_percentage > 60:
            outfile.write(f"{data['line']}\n")

print(f"High overlap hits have been written to: {output_file}")

# Create the scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(
    [data['enhancer_length'] for data in enhancer_overlaps.values()],
    [(data['total_overlap'] / data['enhancer_length']) * 100 if data['enhancer_length'] > 0 else 0 for data in enhancer_overlaps.values()],
    alpha=0.6
)
plt.xlabel('Enhancer Size (bp)')
plt.ylabel('Overlap Percentage (%)')
plt.title('Enhancer Size vs. Overlap Percentage')
plt.grid(True)
plt.tight_layout()
plt.show()
