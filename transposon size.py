import pandas as pd
import matplotlib.pyplot as plt

# Define the path to your BED file
bed_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_library_DM6/dm6_Final_combinedTEs/combined_TE.bed"

# Read the BED file assuming it's tab-delimited and has no headers
bed_df = pd.read_csv(bed_file_path, sep='\t', header=None)

# Calculate the differences between the end (column 3) and start (column 2)
# Columns are zero-indexed, so column 2 is start and column 3 is end
bed_df['difference'] = bed_df[2] - bed_df[1]

# Find rows with a difference greater than 15,000
large_differences = bed_df[bed_df['difference'] > 15000]
small_differences = bed_df[bed_df['difference'] < 20]
# Print the coordinates with differences greater than 15,000
print("Coordinates with differences greater than 15,000 bp:")
print(large_differences[[0, 1, 2, 'difference']])

# Calculate and print the range of sizes (minimum and maximum differences)
min_difference = bed_df['difference'].min()
max_difference = bed_df['difference'].max()
median_difference = bed_df['difference'].median()

print(f"\nRange of sizes:")
print(f"Minimum difference: {min_difference} bp")
print(f"Maximum difference: {max_difference} bp")
print(f"Median difference: {median_difference} bp")
print(f"small difference: {small_differences} bp")
# Plot the histogram of these differences
plt.figure(figsize=(10, 6))
plt.hist(bed_df['difference'], bins=50, edgecolor='black', alpha=0.7)
plt.axvline(median_difference, color='red', linestyle='--', linewidth=2, label=f'Median: {median_difference} bp')
plt.title('TE size distribution',fontsize=16)
plt.xlabel('Size (bp)',fontsize=15)
plt.ylabel('Frequency',fontsize=15)
plt.grid(True)
plt.legend()
plt.show()
