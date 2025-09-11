import pandas as pd

# Load the two files into pandas DataFrames
file1_path = '/Users/jillianness/Desktop/Old_dont_useFly_SE_birth/cannavo comparison BLAST output/Comparisons_modorig001.csv'
file2_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv'

df1 = pd.read_csv(file1_path, header=None)
df2 = pd.read_csv(file2_path, header=None)

# Extract the second column (index 1) from both DataFrames
col1 = df1.iloc[:, 1]
col2 = df2.iloc[:, 1]

# Compare strings in the second column of both files
matches = col1[col1.isin(col2)]
non_matches_file1 = col1[~col1.isin(col2)]
non_matches_file2 = col2[~col2.isin(col1)]

# Display the results
print("Matching strings:")
print(matches)

print("\nNon-matching strings in file 1:")
print(non_matches_file1)

print("\nNon-matching strings in file 2:")
print(non_matches_file2)
