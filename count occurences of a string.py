import csv
from collections import defaultdict

# Replace 'your_file.csv' with the path to your actual CSV file
file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_splitting/output_splitingmergedTE_enhancer.bed"
# Initialize a dictionary to count the occurrences of each string
string_counts = defaultdict(int)

# Open the CSV file and read the first column
with open(file_path, 'r', newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.reader(csvfile,delimiter='\t')
    next(reader)  # Skip the header row if there is one
    for row in reader:
        string_counts[row[4]] += 1

# Initialize counters for the frequency ranges
count_1 = 0
count_2 = 0
count_3 = 0
count_4 = 0
count_5_or_more = 0
print(string_counts.items())
# Count how many strings occur 1, 2-4, or more than 5 times
for count in string_counts.values():
    if count == 1:
        count_1 += 1
    elif count == 2:
        count_2 += 1
    elif count == 3:
        count_3 += 1
    elif count == 4:
        count_4 += 1
    elif count >= 5:
        count_5_or_more += 1

print(f'Strings that occur exactly once: {count_1}')
print(f'Strings that occur 2 times: {count_2}')
print(f'Strings that occur 3 times: {count_3}')
print(f'Strings that occur 4 times: {count_4}')
print(f'Strings that occur more than 5 times: {count_5_or_more}')
print(len(string_counts.values()))