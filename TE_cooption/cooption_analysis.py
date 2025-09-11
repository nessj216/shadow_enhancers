from collections import Counter

# Input file path
input_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/single_enhancer_cooption/test/filtered_shadow_output_overlap.bed'
output_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/enhancer_overlabs_uniq.bed'
# Store unique lines
unique_lines = set()
col4_strings = []
col3_strings = []
counter=0
# Read the input file and store unique lines and column 4 values
with open(input_file, 'r') as infile:
    for line in infile:
        # Remove any extra whitespace and skip empty lines
        line = line.strip()
        counter+=1
        #print(line)
        if line:
            if line not in unique_lines:
                unique_lines.add(line)
                # Extract the value in column 4 (0-indexed: column index 3)
                columns = line.split('\t')
                if len(columns) > 3:
                    col4_strings.append(columns[3])

# Count occurrences of each unique string in column 4
col4_counts = Counter(col4_strings)

# Print the counts
print("Counts of strings in column 4:")
for string, count in col4_counts.items():
    print(f"{string}: {count}")

# Count the number of strings that appear more than once
more_than_once_count = sum(1 for count in col4_counts.values() if count > 1)
print(f"Number of unique column 4 strings that appear more than once: {more_than_once_count}")
# Write the unique lines to the output file
# with open(output_file, 'w') as outfile:
#     for line in unique_lines:
#         outfile.write(f"{line}\n")

print(f"Unique lines have been written to: {output_file}")
print("number of sets w a hut:",len(set(col4_strings)))
#print(set(col4_strings))
#print("number of sequences with TEs",len(unique_lines))
print(counter)

#input_file = "your_file.txt"  # Replace with your file path

# Initialize sets
unique_lines = set()
unique_col4_strings = set()
col_1_2_3=set()
# Open the file and process each line
with open(input_file, 'r') as infile:
    for line in infile:
        # Remove any extra whitespace and skip empty lines
        line = line.strip()
        if line:
            unique_lines.add(line)  # Add the entire line to unique lines

            # Split the line into columns (assuming tab-delimited)
            columns = line.split('\t')
            if len(columns) >= 4:  # Ensure column 4 exists
                unique_col4_strings.add(columns[3])  # Add column 4 to the set

                col_1_2_3.add(tuple(columns[:3]))
# Count the results
unique_line_count = len(unique_lines)
unique_col4_count = len(unique_col4_strings)

# Print results
print(f"Number of unique lines: {unique_line_count}")
print(f"Number of unique column 4 strings: {unique_col4_count}")
print(f"number of sequences col123:", len(col_1_2_3))