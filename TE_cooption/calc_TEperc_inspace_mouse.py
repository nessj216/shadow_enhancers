import pandas as pd

# Load the file (tab-separated)
df = pd.read_csv("/Users/jillianness/Desktop/mouse_analysis_031925/splitting_analysis/filter_out_cooptionhits/splittingresults_filteredcooption.txt", sep="\t", header=None)

# Rename columns for clarity
df.columns = ['chr1', 'start1', 'end1', 'name1', 'name2', 'chr2', 'start2', 'end2', 'value']

# Group by the unique set of columns 1–4
grouped = df.groupby(['chr1', 'start1', 'end1', 'name1'])

# Calculate sum of col9 (value) and enhancer length (end1 - start1)
summary = grouped['value'].sum().reset_index(name='total_col9')
summary['enhancer_length'] = summary['end1'] - summary['start1']
summary['percent'] = summary['total_col9'] / summary['enhancer_length'] * 100

# Save to CSV
summary[['chr1', 'start1', 'end1', 'name1', 'total_col9', 'enhancer_length', 'percent']].to_csv("/Users/jillianness/Desktop/mouse_analysis_031925/splitting_analysis/filter_out_cooptionhits/enhancer_percent_perc.csv", index=False)

print(summary.head())

# Filter sets with percent > 60%
filtered = summary[summary['percent'] > 80]

# Save filtered output as BED file
filtered[['chr1', 'start1', 'end1', 'name1', 'percent' ]].to_csv("/Users/jillianness/Desktop/mouse_analysis_031925/splitting_analysis/filter_out_cooptionhits/filtered_enhancers_80percent.bed", sep='\t', header=False, index=False)

print(f"Saved {len(filtered)} enhancer sets with >60% coverage to 'filtered_enhancers_0percent.bed'")