# import csv
#
# # Replace 'your_file.csv' with the path to your actual CSV file
# file_path = '/Users/jillianness/Desktop/all_mouse_seq0001/CollatedComparisons_0001.csv'
#
# # Open the CSV file and read the first column
# with open(file_path, 'r', newline='', encoding='utf-8-sig') as csvfile:
#     reader = csv.reader(csvfile)
#     # Extract the first column for each row in the CSV
#     first_column = [row[3] for row in reader]
#
# # Count the unique strings in the first column, excluding the header if it exists
# unique_strings = len(set(first_column)) - 1  # Assuming the first row is a header
#
# print(f'Number of unique strings in the first column: {unique_strings}')

########
import csv

# Replace 'your_file.csv' with the path to your actual CSV file
file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv'
# Open the CSV file and read the eighth column
with open(file_path, 'r', newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    # Extract the eighth column for each row in the CSV
    eighth_column = [row[0] for row in reader]

# Count the unique strings in the eighth column, excluding the header if it exists
unique_strings = len(set(eighth_column)) - 1  # Assuming the first row is a header

print(f'Number of unique strings in the 2nd column: {unique_strings}')

import csv
from collections import Counter

# Replace 'your_file.csv' with the path to your actual CSV file
file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv'


# Open the CSV file and read the second column
with open(file_path, 'r', newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    # Skip the header row and extract the second column
    second_column = [row[1] for i, row in enumerate(reader) if i > 0 and len(row) > 1]

# Count occurrences of each string in the second column
string_counts = Counter(second_column)

# Filter for strings that appear more than once
duplicate_strings = {string: count for string, count in string_counts.items() if count > 1}

print(f'Number of strings that appear more than once in the second column: {len(duplicate_strings)}')
print("Duplicate strings and their counts:", duplicate_strings)


