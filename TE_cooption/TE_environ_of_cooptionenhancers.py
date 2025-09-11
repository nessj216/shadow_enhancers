
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Load the two files
file1_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_envi_singles/v3_merged_single_env__TE_overlap.bed"
file2_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/single_enhancer_cooption/single_50bpFiltered_TE_Hits.csv"

# Load first file (no header, tab-separated)
file1_df = pd.read_csv(file1_path, sep="\t", header=None)

# Load second file (assuming it has a header)
file2_df = pd.read_csv(file2_path)

# Extract only the first three columns for comparison
file1_subset = file1_df.iloc[:, :3]
file2_subset = file2_df.iloc[:, :3]

# Find matching rows where the first three columns match
matching_rows = file1_df[file1_subset.apply(tuple, axis=1).isin(file2_subset.apply(tuple, axis=1))]

# Save the filtered file
filtered_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_envi_singles/TE_env_ofsinglecooptionhits.txt"
matching_rows.to_csv(filtered_file_path, sep="\t", index=False, header=False)

