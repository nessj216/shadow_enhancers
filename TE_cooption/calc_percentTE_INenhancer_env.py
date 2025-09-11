import pandas as pd

# Load the file (adjust the file_path as needed)
#file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_env_shadows/overlap_TE_shadows.bed"
import pandas as pd

# Load the file (adjust the file_path as needed)
file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_env_shadows/overlap_merged_shadowstes.txt"
df = pd.read_csv(file_path, sep="\t", header=None)

# Ensure that start (col 1) and end (col 2) are numeric
df[1] = pd.to_numeric(df[1])
df[2] = pd.to_numeric(df[2])
df[6] = pd.to_numeric(df[6])  # Ensure the column to be summed is numeric

# Calculate enhancer length from start and end positions
df['enhancer_length'] = df[2] - df[1]

# Group by chr (col 0), start (col 1), and end (col 2), summing column 4
grouped = df.groupby([0, 1, 2]).agg({
    3: 'first',            # Keep the first Enhancer ID
    6: 'sum',              # Sum column 4
    'enhancer_length': 'first'  # Keep the first enhancer length
}).reset_index()

# Rename the columns for clarity.
grouped = grouped.rename(columns={
    0: "chr",
    1: "start",
    2: "end",
    3: "Enhancer_ID",      # Enhancer ID from col 3
    6: "Total_Col4"
})

# Calculate the percentage of column 4 relative to enhancer length
grouped['percent_col4'] = (grouped['Total_Col4'] / grouped['enhancer_length']) * 100

# Select only the columns we want in the output
output_df = grouped[['chr', 'start', 'end', 'Enhancer_ID', 'percent_col4']]

# Save the final output to a new file
output_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_env_shadows/percentsums_merged_shadows.txt"
output_df.to_csv(output_path, sep="\t", header=True, index=False)

print(f"Modified file saved as: {output_path}")
