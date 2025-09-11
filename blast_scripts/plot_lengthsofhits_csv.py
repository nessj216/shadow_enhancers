import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv'  # Replace with your CSV file path

data = pd.read_csv(file_path)

# Extract the 6th column (by index) and ensure it's treated as numeric
sixth_column = pd.to_numeric(data.iloc[:, 5], errors='coerce')

# Drop any NaN values that result from non-numeric entries
sixth_column = sixth_column.dropna()

# Plot the distribution
plt.figure(figsize=(10, 6))
plt.hist(sixth_column, bins=20, edgecolor='black')
plt.xlabel(' Hit Length (bps)',fontsize=18)
plt.ylabel('% Coverage',fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('% Hit Coverage in the enhancer', fontsize=20)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = '/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv'  # Replace with your CSV file path
data = pd.read_csv(file_path)

# Extract the 6th column (by index) and ensure it's treated as numeric
sixth_column = pd.to_numeric(data.iloc[:, 5], errors='coerce')

# Drop any NaN values that result from non-numeric entries
sixth_column = sixth_column.dropna()

# Calculate the median
median_value = sixth_column.median()

# Plot the distribution
plt.figure(figsize=(10, 6))
plt.hist(sixth_column, bins=20, edgecolor='black')
plt.axvline(median_value, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')
plt.xlabel('Hit Length (bps)', fontsize=20)
plt.ylabel('Frequency', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title('Distribution of BLAST hit lengths (bps)', fontsize=25)
plt.legend(fontsize=20)
plt.show()

print(f"The median hit length is: {median_value:.2f}")
