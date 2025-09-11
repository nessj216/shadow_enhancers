'''Analyzes CollatedComparison file to find different way to count number of hits'''
import pandas as pd
'''1. first just count number of genes with hits '''
# Load the CSV file
data = pd.read_csv('/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons_f.csv')
unique_gene_names = data.iloc[:, 0].nunique()
print('unqiue gene names/number of sets with hits:', unique_gene_names)
# Group by column 1 (Gene name) and count unique occurrences in column 2 (number of pairs) for each value in column 1
unique_col2_counts = data.groupby(data.columns[0])[data.columns[1]].nunique()




'''
2. next count unique pairs with hits/set 
Categorize the counts of unique Column 2 values for each Column 1 value



so for ex:

input
Col1	Col2
A	X
A	Y
A	X
B	Z
B	Z


output: 
Col1	Unique Counts in Col2
A	2
B	1
'''


unique_counts_summary = {
    '1 occurrence': (unique_col2_counts == 1).sum(),
    '2-4 occurrences': ((unique_col2_counts >= 2) & (unique_col2_counts <= 4)).sum(),
    '>=5 occurrences': (unique_col2_counts >= 5).sum()
}

print("number of pairs with hits/set:",unique_counts_summary)  #this plots the number of shadow sets (so 44 shadow sets have hits/409 shadowsets)




'''3. Counts the hits/pair
so counts the number of times a specific shadow pair appears in the Comparisons column of BLAST output CollatedComparison file'''


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

print('number of hits/pair count:',counts_summary)