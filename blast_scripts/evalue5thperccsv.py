import csv
import statistics

import numpy as np


def calculate_percentile(filename, column_index):
    values = []

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header row

        for row in reader:
            if row[column_index]:  # Ignore blank entries
                values.append(float(row[column_index]))

    percentile = np.percentile(values, 10)
    print(percentile)
    index = range(len(values))

    values.sort()
    import matplotlib.pyplot as plt
    plt.plot(values, index, marker='o')
    plt.axhline(percentile, color='r', linestyle='--', label='5th Percentile')
    plt.xlabel('value')
    plt.ylabel('index')
    plt.title('Sorted Values')
    plt.grid(True)
    plt.show()
    return percentile
# Example usage
filename = "/Users/jillianness/Desktop/random_mouse/eavlues.csv"  # Replace with your CSV file path
column_index = 0  # Replace with the index of the desired column (zero-based)

percentile_value = calculate_percentile(filename, column_index)
print(f"The 5th percentile value is: {percentile_value}")
