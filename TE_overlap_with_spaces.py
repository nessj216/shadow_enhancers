input_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_splitting/output_splitting.bed"
output_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_splitting/output_splitting_60cov.bed"

import matplotlib.pyplot as plt


def calculate_overlap_percentage(te_start, te_end, space_start, space_end):
    overlap_start = max(te_start, space_start)
    overlap_end = min(te_end, space_end)
    overlap_length = max(0, overlap_end - overlap_start)
    space_length = space_end - space_start
    return (overlap_length / space_length)
 # * 100 if enhancer_length > 0 else 0



percentages = []
te_sizes = []
space_sizes = []
large_enhancer_ids = []
TEandpercet=[]
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.strip():  # Skip empty lines
            columns = line.strip().split('\t')
            te_start = int(columns[6])
            te_end = int(columns[7])
            space_start = int(columns[1])
            space_end = int(columns[2])
            TEgenename=(columns[3])

            # Calculate sizes
            te_size = te_end - te_start
            space_size = space_end - space_start

            # Store sizes
            te_sizes.append(te_size)
            space_sizes.append(space_size)

            # Calculate overlap percentage
            percent_overlap = calculate_overlap_percentage(te_start, te_end, space_start, space_end)
            percentages.append(percent_overlap)
            TEandpercet1=(te_size,percent_overlap)
            TEandpercet.append(TEandpercet1)
            # Write to output file
            outfile.write(f"{line.strip()}\t{percent_overlap:.2f}%\n")

            # # Check if enhancer space is greater than 150 kb
            # if enhancer_size > 150000:
            #     enhancer_id = columns[10]  # Column 11 contains the enhancer ID
            #     large_enhancer_ids.append(enhancer_id)
            # ##curious what if i plot overlap y axis and TE size on the x axis

print(percentages)
# Plotting the percentage overlap histogram
plt.figure(figsize=(10, 6))
plt.hist(percentages, bins=20, edgecolor='black')
plt.title('Distribution of TE Overlap Percentages')
plt.xlabel('Percentage Overlap')
plt.ylabel('Frequency')
plt.show()

# Plotting the TE size histogram
plt.figure(figsize=(10, 6))
plt.hist(te_sizes, bins=20, edgecolor='black')
plt.title('Distribution of TE Sizes')
plt.xlabel('TE Size (bp)')
plt.ylabel('Frequency')
plt.show()

# Plotting the enhancer space size histogram
plt.figure(figsize=(10, 6))
plt.hist(space_sizes, bins=20, edgecolor='black')
plt.title('Distribution of Enhancer Space Sizes')
plt.xlabel('Enhancer Space Size (bp)')
plt.ylabel('Frequency')
plt.show()

plt.figure(figsize=(10, 6))
plt.scatter(te_sizes, space_sizes)
plt.xlabel('te_sizes')
plt.ylabel('enhancer_space_sizes')
plt.show()

plt.figure(figsize=(10, 6))
plt.scatter(te_sizes,percentages)
plt.xlabel('te_sizes')
plt.ylabel('percentage of space')
plt.show()


# Count the number of overlaps greater than 50%
greater_than_50 = sum(1 for p in percentages if p > 50)
print(f"Number of TE overlaps greater than 50%: {greater_than_50}")

# Print enhancer IDs with spaces larger than 150 kb
if large_enhancer_ids:
    print(f"Enhancer IDs with spaces larger than 150 kb:")
    for enhancer_id in large_enhancer_ids:
        print(enhancer_id)
else:
    print("No enhancer spaces larger than 150 kb found.")
