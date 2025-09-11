import csv
import matplotlib.pyplot as plt
import numpy as np
import math
# Replace with your CSV file path
csv_file_path = "/Users/jillianness/Desktop/Old_dont_useFly_SE_birth/KvonrandomBLAST/CollatedComparisons_final.csv"

# List to store the values from the 11th column
column_11_values = []
column_11_eval = []
# Read the CSV file and extract the 11th column values
with open(csv_file_path, "r") as csvfile:
    csvreader = csv.reader(csvfile)
    header = next(csvreader)  # Skip the header row
    for row in csvreader:
        if len(row) >= 10:
            column_11_value=row[10].strip('\n')
            #print(column_11_value)
            column_11_value = float(row[10])
            column_11_values.append(((column_11_value)))
            column_11_eval.append(((column_11_value)))

# Create a histogram
plt.hist(column_11_values, bins=100,color='blue')
plt.xlabel("evalue")
plt.ylabel("Frequency")
plt.title("Histogram of e-value")
plt.grid(True)

# Calculate the 90th percentile
percentile_90 = np.percentile(column_11_values, 10)
formatted_percentile_90 = "{:.4f}".format(percentile_90)
percentile_95 = np.percentile(column_11_values, 95)
print(10**-percentile_90,10**-percentile_95)
# Add a red vertical line at the 90th percentile
plt.axvline(x=percentile_90, color='red', linestyle='--', label='(90% -log(eval)'+formatted_percentile_90)
print(np.percentile(column_11_eval, 10))
print(len(column_11_values))
# Add a legend
plt.legend()

# Show the histogram
plt.show()