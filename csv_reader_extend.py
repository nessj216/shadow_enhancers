import os
import shutil
import csv
import numpy as np
column_11_values=[]
pathTillGene = "/Users/jillianness/Desktop/fly_SE_birth_dm3/SEoutput_blast_wseqhits"
file_list = os.listdir(pathTillGene)

onlyFolder = []

for item in file_list:
    full_path = os.path.join(pathTillGene, item)

    if os.path.isdir(full_path):
        onlyFolder.append(os.path.join(full_path, 'Comparisons_mod001'))

# Remove the old CSV file
csv_file_path = os.path.join(pathTillGene, 'Comparisons_mod001.csv')
if os.path.exists(csv_file_path):
    os.remove(csv_file_path)

with open(csv_file_path, mode='w', newline='') as csv_file:
    csvFileWriter = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csvFileWriter.writerow(
        ['Gene Name', 'Comparisons', 'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send',
         'evalue'])

    for compLoc in onlyFolder:
        onlyComp = os.listdir(compLoc)

        for compFile in onlyComp:
            with open(os.path.join(compLoc, compFile), "r") as f:
                data = f.readlines()

                for line in data:
                    line = line.strip()
                    if line:
                        splitData = line.split('\t')  # Split by tab character
                        geneName = compLoc.split("/")[-2]
                        splitData.insert(0, compFile)
                        splitData.insert(0, geneName)
                        csvFileWriter.writerow(splitData)
                    # Extract the 11th column value and store it in the list
                    if len(splitData) > 10:

                        column_11_value = float(splitData[10])
                        column_11_values.append(column_11_value) #list of e values from hits
print(column_11_values) #list of e values from hits
# Calculate the 90th percentile of the 11th column values using NumPy
percentile_90 = np.percentile(column_11_values, 5)
print(f"90th percentile of the 11th column: {percentile_90}")
print(len(column_11_values))