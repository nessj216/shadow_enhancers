import os
import re

# Define TF list
TF_list = ['cad', 'hb', 'bcd', 'kni', 'tll', 'Stat92E', 'gt', 'zelda', 'Kr']

# Define enhancer parameters
start_ = 20689644
end_ = 20690670
chr_ = 'chr3L'
f_name = 'kni2'

# Create a dictionary to store data for each key
data = {}

# Find matching files
for root, _, files in os.walk('.'):
    for file in files:
        if file.endswith("_bound.bed"):
            file_name = os.path.basename(file)

            # Find the matching TF in the file name
            tf_name = next((tf for tf in TF_list if tf in file_name), None)

            if tf_name:
                file_path = os.path.join(root, file)

                # Process the file
                with open(file_path, 'r') as f:
                    for line in f:
                        columns = line.strip().split('\t')
                        if (
                            columns[1] >= start_ and columns[1] <= end_
                            and columns[2] >= start_ and columns[2] <= end_
                            and columns[0] == chr_
                        ):
                            key = ''.join(columns[i] for i in [0, 1, 2, 4, 5, 6])
                            if key not in data:
                                data[key] = {
                                    'sum': 0,
                                    'count': 0,
                                    'lines': [],
                                }
                            if columns[5] == '+':
                                val = columns[1] - start_
                            elif columns[5] == '-':
                                val = columns[1] - start_
                            else:
                                val = ''
                            _val = columns[1] - columns[2]
                            abs_val = abs(_val)
                            data[key]['sum'] += float(columns[3])
                            data[key]['count'] += 1
                            data[key]['lines'].append(
                                f"binding_site {columns[1]} {abs_val} {columns[16]} {val} {tf_name} {columns[5]}",
                            )

# Write the output files
with open(f_name + "myconcatenated_files.txt", 'w') as output_file:
    for key in data:
        if data[key]['count'] >= 2:
            avg = data[key]['sum'] / data[key]['count']
            for line in data[key]['lines']:
                output_file.write(f"{avg} {line}\n")

print("Concatenation completed. Result saved in", f_name + "myconcatenated_files.txt")
