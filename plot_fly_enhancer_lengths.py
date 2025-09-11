import os
import statistics

import pandas as pd
import matplotlib.pyplot as plt

# Directory containing the .txt files
filepath = "/Users/jillianness/Desktop/mouse_analysis_031925/Mouse_enhancers/Final_enhancers_deduplicated.bed"


# List to store differences
differences = []
count=0
# Open and read the file line by line
with open(filepath, 'r') as file:
    for line in file:
        columns = line.strip().split()
        if len(columns) >= 3:
            try:
                start = int(columns[1])
                end = int(columns[2])
                difference = end - start
                differences.append(difference)

                if difference >= 40000:
                    print(line, difference)
                    count+=1
            except ValueError:
                # Skip lines where conversion to int fails
                continue
print(f"Number of entries:{count}")
# Print summary stats
print(f"Number of entries: {len(differences)}")
print(f"Median length: {statistics.median(differences)}")
print(f"Total sum of lengths: {sum(differences)}")

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(differences, bins=20, edgecolor='black')
plt.xlabel('Length (bps)')
plt.ylabel('Frequency')
plt.title('Histogram of Enhancer Sizes')
plt.tight_layout()
plt.show()