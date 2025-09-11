import pandas as pd

# Load the file into a DataFrame
file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/enhancer_overlabs.bed'

data = pd.read_csv(file_path, sep='\t', header=None)

# Extract column 4 (index 3) and count unique values
unique_strings = data[3].nunique()
# Count unique strings in the first column, counting duplicates only once
unique_gene_names = data['Gene Name'].nunique()
# Extract columns 1, 2, and 3 (index 0, 1, 2) and count unique sets
unique_col123_sets = data[[0, 1, 2]].drop_duplicates().shape[0]

print(f"Number of unique strings in column 4: {unique_strings}")
print(f"Number of unique (col1, col2, col3) sets: {unique_col123_sets}")
