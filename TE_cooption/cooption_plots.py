import matplotlib.pyplot as plt
import statistics
# Input file
input_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/cooption_enhancer_TE.bed"
##input file is merged list
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
            te_start = int(columns[5])
            te_end = int(columns[6])

            # Define enhancer key based on chr, start, end
            enhancer_key = (chr_name, enhancer_start, enhancer_end)
            enhancer_length = enhancer_end - enhancer_start

            # Calculate overlap length
            overlap_length = calculate_overlap_length(te_start, te_end, enhancer_start, enhancer_end)

            # Store overlap for the given enhancer
            if enhancer_key not in enhancer_overlaps:
                enhancer_overlaps[enhancer_key] = {'total_overlap': 0, 'enhancer_length': enhancer_length}

            enhancer_overlaps[enhancer_key]['total_overlap'] += overlap_length

# Calculate overlap percentages
enhancer_sizes = []
overlap_percentages = []

for data in enhancer_overlaps.values():
    enhancer_length = data['enhancer_length']
    total_overlap = data['total_overlap']
    overlap_percentage = (total_overlap / enhancer_length) * 100 if enhancer_length > 0 else 0

    enhancer_sizes.append(enhancer_length)
    overlap_percentages.append(overlap_percentage)
    # Detect high overlaps
    if overlap_percentage > 60:
        print(f"Enhancer: {chr_name}, {enhancer_start}, {enhancer_end}")


# Calculate and print the median of overlap percentages
median_overlap = statistics.median(overlap_percentages)
print(f"Median overlap percentage: {median_overlap}")

# Create the scatter plot
plt.figure(figsize=(12, 6))

# Scatter plot of enhancer sizes vs. overlap percentages
plt.subplot(1, 2, 1)
plt.scatter(enhancer_sizes, overlap_percentages, alpha=0.6)
plt.xlabel('Enhancer Size (bp)', fontsize=14)  # Adjust font size here
plt.ylabel('Overlap Percentage (%)', fontsize=14)  # Adjust font size here
plt.title('Enhancer Size vs. Overlap Percentage', fontsize=16)  # Adjust font size here
plt.grid(True)


# Histogram of overlap percentages
plt.subplot(1, 2, 2)
plt.hist(overlap_percentages, bins=30, alpha=0.7, edgecolor='black')
plt.xlabel('Overlap Percentage (%)', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.title('Distribution of Overlap Percentages', fontsize=16)
plt.grid(True)

plt.tight_layout()
plt.show()
