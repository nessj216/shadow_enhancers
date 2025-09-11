'''Counts the hits/pair of hits
so counts the number of times a shadow pair appears in the Comparisons column '''


import pandas as pd

# Load the CSV file
data = pd.read_csv('/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/tfbs_scan_v2/CollatedComparisons_f.csv')

# Assuming the second column is at index 1
value_counts = data.iloc[:, 1].value_counts()

# Classify counts
counts_summary = {
    '1 occurrence': (value_counts == 1).sum(),
    '2-4 occurrences': ((value_counts >= 2) & (value_counts <= 4)).sum(),
    '>=5 occurrences': (value_counts >= 5).sum()
}

print(counts_summary)
