# Define a function to extract the desired column from a line
def get_column(line, column_number, delimiter='\t'):
    columns = line.strip().split(delimiter)
    return columns[column_number - 1]  # Adjusting to 0-based index

# Define the file paths
file1_path = '/Users/jillianness/Desktop/single_enhancer_cannavo/non_matching_lines_file2.txt'
file2_path = '/Users/jillianness/Desktop/single_enhancer_cannavo/cannavo_allenhancers.v3.txt'

# Load lines from file1 and store column 1 values in a set
with open(file1_path, 'r') as file1:
    column1_values = set(get_column(line, 1) for line in file1)

# Extract lines from file2 where column 4 matches column 1 of file1
matching_lines_file2 = []
with open(file2_path, 'r') as file2:
    for line in file2:
        if get_column(line, 4) in column1_values:
            matching_lines_file2.append(line.strip())

# Print matching lines from file2
'''print("Matching lines from file2:")
for line in matching_lines_file2:
    print(line)'''
print (len(matching_lines_file2))

with open('/Users/jillianness/Desktop/single_enhancer_cannavo/matching_lines_file2.txt', 'w') as output_file:
    for line in matching_lines_file2:
        output_file.write(line + '\n')