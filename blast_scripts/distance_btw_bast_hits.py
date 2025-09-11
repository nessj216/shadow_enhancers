import csv
import statistics
import matplotlib.pyplot as plt
import numpy as np

# Initialize lists to hold distances and enhancer lengths
distance = []
lengthenhancer = []

# Replace 'yourfile.csv' with the path to your actual CSV file
file_path = "/Users/jillianness/Desktop/all_mouse_seq0001_seq/coordinates_uniq.csv"
file_path2 = "/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv"

# Process the file to calculate distance and enhancer length based on 'qseqid' and 'sseqid' columns
with open(file_path2, newline='') as csvfile:
    reader = csv.DictReader(csvfile)  # Use DictReader to access columns by name
    for row in reader:
        # Access the 'qseqid' and 'sseqid' columns
        col1, col2 = row['qseqid'], row['sseqid']

        # Skip rows with empty values in 'col1' or 'col2'
        if not col1 or not col2:
            continue

        # Split the coordinates
        col1_numbers = col1.split(":")[1].split("-")
        col2_numbers = col2.split(":")[1].split("-")
        col1_first_num, col1_second_num = int(col1_numbers[0]), int(col1_numbers[1])
        col2_first_num, col2_second_num = int(col2_numbers[0]), int(col2_numbers[1])

        # Calculate max values for each column
        max_col1 = max(col1_first_num, col1_second_num)
        max_col2 = max(col2_first_num, col2_second_num)

        # Calculate the distance and store it in the distance list
        dist = abs(max_col2 - max_col1)
        distance.append(dist)

# Calculate median, minimum, and maximum distance if there are valid distances
if distance:
    median_distance = statistics.median(distance)
    min_distance = min(distance)
    max_distance = max(distance)

    print("Median Distance:", median_distance)
    print("Minimum Distance:", min_distance)
    print("Maximum Distance:", max_distance)
else:
    print("No valid distances calculated.")


