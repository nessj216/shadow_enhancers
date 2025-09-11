import pandas as pd
import zipfile
import os

# --- Step 1: Load files ---

# File paths
first_file = '/Users/jillianness/Desktop/mouse_analysis_031925/splitting_analysis/filter_out_cooptionhits/cooptinfile_50bp.bed'
second_file = '/Users/jillianness/Desktop/mouse_analysis_031925/splitting_analysis/space_TEoutput.bed'
#extract_folder = '/mnt/data/space_TEoutput_extracted'



# Read files
first_df = pd.read_csv(first_file, sep='\t', header=None)
second_df = pd.read_csv(second_file, sep='\t', header=None)

# --- Step 2: Filter based on matching col5 and col4 ---

# For each row in second_df, check:
#  - If any col5 from first_df is a substring of this row's col5
#  - AND col4 matches exactly

# Group substrings by col4 for fast matching
from collections import defaultdict

lookup = defaultdict(list)
for _, row in first_df.iterrows():
    lookup[row[3]].append(row[4])

# --- Step 3: Filtering without slow row-by-row apply ---

# Prepare a boolean mask
keep_rows = []

for idx, row in second_df.iterrows():
    col4 = row[3]
    col5 = row[4]
    substrings = lookup.get(col4, [])  # Get all substrings for this col4

    # Check if any substring is inside col5
    if any(substr in col5 for substr in substrings):
        keep_rows.append(False)  # REMOVE this row
    else:
        keep_rows.append(True)   # KEEP this row

# Apply the mask
filtered_second_df = second_df[keep_rows]

# --- Step 4: Save the output ---
output_path = '/Users/jillianness/Desktop/mouse_analysis_031925/splitting_analysis/filter_out_cooptionhits/splittingresults_filteredcooption.txt'
filtered_second_df.to_csv(output_path, sep='\t', header=False, index=False)

print(f"Filtered file saved to: {output_path}")