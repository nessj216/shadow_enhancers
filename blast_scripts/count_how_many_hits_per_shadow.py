import pandas as pd

# Load the CSV file
data = pd.read_csv('/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv')

# Group by column 1 (Gene name) and count unique occurrences in column 2 (number of pairs) for each value in column 1
unique_col2_counts = data.groupby(data.columns[0])[data.columns[1]].nunique()
print(unique_col2_counts)
# Categorize the counts of unique Column 2 values for each Column 1 value
unique_counts_summary = {
    '1 occurrence': (unique_col2_counts == 1).sum(),
    '2-4 occurrences': ((unique_col2_counts >= 2) & (unique_col2_counts <= 4)).sum(),
    '>=5 occurrences': (unique_col2_counts >= 5).sum()
}

print(unique_counts_summary)  #this plots the number of shadow sets (so 44 shadow sets have hits/409 shadowsets)
