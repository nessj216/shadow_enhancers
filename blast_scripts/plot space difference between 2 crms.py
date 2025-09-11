import pandas as pd
import matplotlib.pyplot as plt
import math
# Your Excel file pat
import numpy as np
import csv

import statistics

distance=[]
# Replace 'yourfile.csv' with the path to your actual CSV file
file_path = "/Users/jillianness/Desktop/all_mouse_seq0001_seq/coordinates_uniq.csv"
file_path2 = "/Users/jillianness/Desktop/SEbirthanalysis1024/Duplications_BLAST/ALL_shadowsets_dm3/CollatedComparisons.csv"
lengthenhancer=[]
with open(file_path2, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')  # Adjust delimiter if necessary
    next(reader)
    for row in reader:
        # Split the row into two columns
        col1, col2 = row[2], row[3]

        # Process the first column
        col1_numbers = col1.split(":")[1].split("-")

        # Process the second column if needed, similar to the first column
        col2_numbers = col2.split(":")[1].split("-")
        col1_first_num = int(col1_numbers[0])
        col2_first_num = int(col2_numbers[0])
        col1_second_num=int(col1_numbers[1])
        col2_second_num = int(col2_numbers[1])

        length1enhancer=col1_second_num-col1_first_num
        length2enhancer=col2_second_num-col2_first_num
        lengthenhancer.append(length1enhancer)
        lengthenhancer.append(length2enhancer)

        #print(col2_first_num,col2_first_num)
        if col2_first_num >= col1_first_num:

        # Subtract the first number of col2 from the first number of col1
            difference = col2_first_num - col1_second_num
            distance.append(difference)
        else:
            difference=col1_first_num - col2_second_num
            distance.append(difference)
filtered_numbers = sorted([num for num in distance if num >= 0])
enhancer_length = sorted([num for num in lengthenhancer if num >= 0])

enhancer_length = sorted([num for num in lengthenhancer if num >= 0])
print('median distance',statistics.median(filtered_numbers))
print('mean distance',statistics.mean(filtered_numbers))

print("Total length of enhancers:", sum(lengthenhancer))
print('median enhancer length',statistics.median(lengthenhancer))



#print(filtered_numbers)

        # Print the processed numbers; adjust this part as per your output needs
        #print(" ".join(col1_numbers))  # , "\t", " ".join(col2_numbers))
# Plotting the histogram

plt.hist(filtered_numbers, bins=20, edgecolor='black')
plt.title('Histogram of distance between crms')
plt.xlabel('distance (bps)')
plt.ylabel('Frequency')
plt.xticks(np.arange(0, 90000, 8000), rotation=90)
#plt.xticks(np.arange(0, 90000))
plt.yticks(np.arange(0, 32, 5))
# Apply these ticks to the x-axis
#plt.xticks(x_ticks)
plt.subplots_adjust(bottom=0.2)
plt.show()


count_less_than_5000 = sum(1 for num in filtered_numbers if num < 10000)
print(count_less_than_5000)

plt.hist(enhancer_length, bins=10, edgecolor='black')
plt.title('Histogram of enhancer length')
plt.xlabel('length (bps)')
plt.ylabel('Frequency')
#x_ticks = np.arange(0, np.max(filtered_numbers) + 200000, 200000)

# Apply these ticks to the x-axis
#plt.xticks(x_ticks)

plt.show()

import pandas as pd
import matplotlib.pyplot as plt

# Replace 'yourfile.csv' with the path to your actual CSV file




#plots length of hit
# Read the CSV file into a DataFrame
df = pd.read_csv(file_path2)

# Plot the 6th column, assuming the first column is indexed as 1
# If the first column is indexed as 0, then you would use df.iloc[:, 5]
plt.figure(figsize=(10, 6))  # Set the figure size as desired
plt.hist(df['length'], bins=10, edgecolor='black') # The 6th column has index 5
plt.title('Histogram of hit length')
plt.xlabel(' length (bps)')
plt.ylabel('frequency')
plt.show()
print(statistics.median(df['length']))
