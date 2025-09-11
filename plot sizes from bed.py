import matplotlib.pyplot as plt

# Read the data from your file
with open('/Users/jillianness/Desktop/fly_blast/repeatmask/canavo_seqs.bed', 'r') as file:
    lines = file.readlines()

# Extract the differences between column 3 and column 2
differences = []
for line in lines:
    columns = line.split('\t')
    if len(columns) >= 4:
        start = int(columns[1])
        end = int(columns[2])
        difference = end - start
        differences.append(difference)

# Create a histogram of the differences with a specified range
plt.hist(differences, bins=range(100, 1250, 20), edgecolor='k')
plt.xlabel('shadow sequence size (bps)')
plt.ylabel('Frequency')
plt.title('Distribution of Canavo et al shadow sequence size')

# Customize the y-axis ticks
#plt.yticks(range(25000, 180001, 25000))

plt.show()
