#calcualtes the percent total TE in each co-opted enhancer
#in this case TE total is in col 8
import pandas as pd

# Load file
file_path = "filtered_50bp_coptionenhancers.bed"
df = pd.read_csv(file_path, sep="\t", header=None)

# Ensure column 8 is numeric
df[7] = pd.to_numeric(df[7], errors='coerce')

# Create a unique identifier for each enhancer
df['enhancer_id'] = df[[0, 1, 2]].astype(str).agg('\t'.join, axis=1)

# Calculate enhancer length
df['length'] = df[2] - df[1]

# Group by enhancer and sum col 8
grouped = df.groupby('enhancer_id').agg({
    0: 'first',
    1: 'first',
    2: 'first',
    7: 'sum',
    'length': 'first'
})

# Compute percentage: sum(col 8) / enhancer length
grouped['percent_col8'] = (grouped[7] / grouped['length']) * 100

# Save to file
output_path = "enhancer_col8_percentages.tsv"
grouped[['percent_col8']].to_csv(output_path, sep="\t")

print(f"Saved enhancer-wise column 8 percentages to {output_path}")
