##this code goes through a directory and its subdirectories to concatenate 20 files at time into a new file...this is
##use REPEATMASKER browser tool which only takes 100 kb files at a time. the command line ReapeatMasker
#needs a linux O.S. and Perl, and a bunch of stuff that i dont want to bother with

import os

# Specify the directory where your subdirectories are located
base_directory = "/Users/jillianness/Desktop/mouse_original_files"

# Create a list to store the names of the files in each subdirectory
all_files = []

# Loop through all subdirectories and collect file names
for root, _, files in os.walk(base_directory):
    print(files)
    for file in files:
        if file.endswith('.txt'):
            all_files.append(os.path.join(root, file))

# Define the number of files to combine at a time
files_to_combine = 20

# Initialize a counter for total processed files
total_processed_files = 0

# Loop through the list of files and combine them in groups of 20
for i in range(0, len(all_files), files_to_combine):
    file_group = all_files[i:i + files_to_combine]

    # Create a new output file name based on the group
    output_file_name = f'combined_output_{i // files_to_combine}.txt'

    # Open the output file in write mode
    with open(output_file_name, 'w') as output_file:
        for file_name in file_group:
            with open(file_name, 'r') as input_file:
                output_file.write(input_file.read())
                total_processed_files += 1

    print(f'Combined {len(file_group)} files into {output_file_name}')

print(f'Total files processed and concatenated: {total_processed_files}')
