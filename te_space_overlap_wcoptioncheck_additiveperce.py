import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

# File paths
input_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_splitting/output_splitingmergedTE_enhancer.bed"
additional_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/all_shadowsets_DM6.bed"  # Replace with the path to your additional file
output_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_splitting/output_splitting_60cov.bed"
cumulative_output_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_splitting/cumulative_overlap_by_col5.bed"

# Function to calculate overlap percentage
def calculate_overlap_percentage(start1, end1, start2, end2):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_length = max(0, overlap_end - overlap_start)
    reference_length = end2 - start2  # Length of the enhancer space in additional data
    return (overlap_length / reference_length) if reference_length > 0 else 0

# Load additional file into a list for overlap and name checking
additional_data = []
with open(additional_file, 'r') as addfile:
    for line in addfile:
        columns = line.strip().split('\t')
        if len(columns) < 4:
            continue
        enhancer_start = int(columns[1])
        enhancer_end = int(columns[2])
        enhancer_name = columns[3]
        additional_data.append((enhancer_start, enhancer_end, enhancer_name))

# Initialize lists and a dictionary to collect cumulative overlap per column 5 identifier
percentages = []
te_sizes = []
space_sizes = []
TEandpercet = []
cumulative_overlaps = defaultdict(float)  # Dictionary to accumulate overlap by column 5

# Read the input file and calculate overlap
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.strip():  # Skip empty lines
            columns = line.strip().split('\t')
            te_start = int(columns[6])
            te_end = int(columns[7])
            space_start = int(columns[1])
            space_end = int(columns[2])
            TEgenename = columns[3]
            col5_value = columns[4]  # Column 5 value to group by

            # Check for overlap and matching TEgenename with 40% threshold
            exclude_entry = False
            for enhancer_start, enhancer_end, enhancer_name in additional_data:
                if TEgenename == enhancer_name:
                    overlap_percentage = calculate_overlap_percentage(te_start, te_end, enhancer_start, enhancer_end)
                    if overlap_percentage > 0.4:  # Exclude if overlap is greater than 40%
                        exclude_entry = True
                        break  # Exit the loop as we found a match to exclude

            if exclude_entry:
                continue  # Skip this entry if it meets exclusion criteria

            # Calculate overlap percentage with enhancer space
            percent_overlap = calculate_overlap_percentage(te_start, te_end, space_start, space_end)
            percentages.append(percent_overlap)

            # Write to output file only if overlap > 50%
            if percent_overlap > 0.5:
                outfile.write(f"{line.strip()}\t{percent_overlap:.2f}%\n")


            # Accumulate the overlap percentage by column 5 identifier
            cumulative_overlaps[col5_value] += percent_overlap

# Write cumulative overlap results to a new file
with open(cumulative_output_file, 'w') as cum_outfile:
    cum_outfile.write("Col5\tCumulative_Overlap_Percentage\n")
    for col5_value, total_overlap in cumulative_overlaps.items():
        cumulative_percent = total_overlap * 100
        if cumulative_percent > 50:  # Only write if cumulative percent is greater than 50%
            cum_outfile.write(f"{col5_value}\t{cumulative_percent:.2f}%\n")


print("Cumulative overlap percentages written to", cumulative_output_file)

# Plot the distribution of cumulative overlap percentages
cumulative_values = [total_overlap * 100 for total_overlap in cumulative_overlaps.values()]

import matplotlib.pyplot as plt
import numpy as np

# Assuming percentages_scaled and cumulative_values have already been computed

# Calculate the median for non-additive and cumulative overlap percentages
percentages_scaled = [p * 100 for p in percentages]  # Convert to percentage format

median_non_additive = np.median(percentages_scaled)
median_cumulative = np.median(cumulative_values)

# Plotting the non-additive overlap percentages histogram with median line
plt.figure(figsize=(10, 6))
plt.hist(percentages_scaled, bins=20, edgecolor='black')
plt.axvline(median_non_additive, color='red', linestyle='dashed', linewidth=1.5, label=f'Median: {median_non_additive:.2f}%')
plt.title('Distribution of Non-Additive % Overlap')
plt.xlabel('% Overlap')
plt.ylabel('Frequency')
plt.legend()  # Show legend with median
plt.show()

# Plotting the cumulative overlap percentages histogram with median line
plt.figure(figsize=(10, 6))
plt.hist(cumulative_values, bins=20, edgecolor='black')
plt.axvline(median_cumulative, color='blue', linestyle='dashed', linewidth=1.5, label=f'Median: {median_cumulative:.2f}%')
plt.title('Distribution of Cumulative % Overlap by Column 5')
plt.xlabel('Cumulative % Overlap')
plt.ylabel('Frequency')
plt.legend()  # Show legend with median
plt.show()

# Counting unique col5 values (ignoring duplicates)
unique_col5_values = len(set(cumulative_overlaps.keys()))
print(f"Number of unique col5 values: {unique_col5_values}")


