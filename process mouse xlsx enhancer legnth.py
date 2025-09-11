import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file
file_path = '/Users/jillianness/Downloads/Osterwalder_et_al_Supplementary_Table_11_modified2.csv'  # Update this with your file path

df = pd.read_csv(file_path, header=0)

# Convert the relevant columns to numeric, forcing errors to NaN
df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')

# Drop rows with NaN values in the relevant columns
df.dropna(subset=[df.columns[1], df.columns[2]], inplace=True)

# Subtract column 3 from column 2 and take the absolute value
df['abs_diff'] = (df.iloc[:, 1] - df.iloc[:, 2]).abs()

# Calculate the sum of all the absolute differences
total_diff = df['abs_diff'].sum()
print(f'The total sum of all absolute differences is: {total_diff}')





# Plot the result in a histogram
plt.figure(figsize=(10, 6))
plt.hist(df['abs_diff'], bins=30, edgecolor='black')
plt.xlabel('Absolute Difference (|Column 2 - Column 3|)')
plt.ylabel('Frequency')
plt.title('Histogram of Absolute Differences between Column 2 and Column 3')
plt.show()
