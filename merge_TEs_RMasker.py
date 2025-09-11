import pandas as pd

# File paths and column names
file_path = '/Users/jillianness/Desktop/sorting_cannavo_data/final_merged_libraries/mask_filtered_outputRM_ONECODEfile.bed'
columns = ['chr', 'start', 'end', 'column4', 'name', 'strand', 'family', 'column8', 'column9']

# Load the file into a pandas dataframe
df = pd.read_csv(file_path, sep='\t', names=columns)

# Convert 'start' and 'end' to numeric for calculation
df['start'] = pd.to_numeric(df['start'], errors='coerce')
df['end'] = pd.to_numeric(df['end'], errors='coerce')

# Helper function to extract string before underscore or use the full string if no underscore
def extract_prefix(name):
    return name.split('_')[0] if '_' in name else name

# Apply the extract_prefix function to the 'name' column
df['name'] = df['name'].apply(extract_prefix)

# Sort the dataframe by chromosome, strand, name, and start position
df_sorted = df.sort_values(by=['chr', 'strand', 'name', 'start']).reset_index(drop=True)

# Initialize a new DataFrame to store the merged results
merged_rows = []

# Initialize variables to store current merging information
current_chr = None
current_strand = None
current_name = None
current_start = None
current_end = None
current_family = None
# Counters to keep track of merges
lines_merged = 0  # Count the number of lines that were merged
merged_entries = 0  # Count the total number of lines in merged_rows

# Iterate through sorted rows and merge rows that overlap or are within the defined distance
for i, row in df_sorted.iterrows():
    space_between = 20
    if (current_chr == row['chr'] and current_strand == row['strand'] and current_name == row['name'] and row['start'] <= current_end + space_between):
        # Merge with current region
        current_end = max(current_end, row['end'])
        lines_merged += 1
    else:
        # If the current region is not None, append it to the results
        if current_chr is not None:
            # Append the merged row to the results
            merged_rows.append([current_chr, current_start, current_end,  current_name, current_strand, current_family])
            merged_entries += 1
            # Start a new region
        current_chr = row['chr']
        current_strand = row['strand']
        current_name = row['name']
        current_start = row['start']
        current_end = row['end']
        current_family = row['family']

# Append the last region
if current_chr is not None:
    merged_rows.append([current_chr, current_start, current_end,current_name, current_strand,current_family])
    merged_entries += 1  # Increment for the last merged entry

# Print the number of lines merged and the number of lines being sent to the file
print(f"Number of lines merged: {lines_merged}")
print(f"Number of lines being sent to the file: {merged_entries}")

# Convert merged results to DataFrame and save the necessary columns
merged_df = pd.DataFrame(merged_rows, columns=['chr', 'start', 'end', 'name', 'strand','family'])

# Save the merged dataframe to a file
output_file = '/Users/jillianness/Desktop/sorting_cannavo_data/final_merged_libraries/Jmerged_mask_filtered_outputRM_ONECODEfile.bed'
merged_df.to_csv(output_file, sep='\t', index=False)
