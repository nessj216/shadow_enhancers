import os
import statistics
import matplotlib.pyplot as plt
import numpy as np
directory = "/Users/jillianness/Desktop/TOBIAS091922/combd1d4"
enhancer_kni0 = {}
enhancer_kni2 = {}
TFscores = {}
tempor = []
tempor1 = {}
TFstats = {}
norm_scores = []

norm_dict1 = {}
norm_dict2 = {}
avg_norm_dict1 = {}
avg_norm_dict2 = {}
sum_norm_scores1 = {}
sum_norm_scores2 = {}
TFs = ['gt', 'cad', 'bcd', 'kni', 'zelda', 'Kr', 'Stat', 'hb', 'tll']

for path, subdir, files in os.walk(directory):
    for file1 in files:
        tempor = []
        filename = str(os.fsdecode(file1))
        filename1 = os.path.join(path, filename)

        if os.path.getsize(filename1) == 0:
            continue

        if 'sogprox' in filename:
            with open(filename1, "r") as t:
                tempor1[path] = []
                for line in t:
                    tempor1[path].append(float(line.strip('\n')))
            enhancer_kni0[path] = tempor1[path]

        if 'sogdis' in filename:
            with open(filename1, "r") as t:
                tempor1[path] = []
                for line in t:
                    tempor1[path].append(float(line.strip('\n')))
            enhancer_kni2[path] = tempor1[path]

        if 'all.bed.txt' in filename:
            with open(filename1, "r") as t:
                for line in t:
                    tempor.append(float(line.rstrip('\n')))
                TFscores[path] = tempor

for key, value in TFscores.items():
    avg = statistics.mean(value)
    sd = statistics.stdev(value)
    TFstats[key] = (avg, sd)

for key in enhancer_kni0:
    norm_scores = []
    for i in enhancer_kni0[key]:
        norm_scores.append((i - TFstats[key][0]) / TFstats[key][1])
    norm_dict1[key] = norm_scores

for key in enhancer_kni2:
    norm_scores = []
    for i in enhancer_kni2[key]:
        norm_scores.append((i - TFstats[key][0]) / TFstats[key][1])
    norm_dict2[key] = norm_scores

for key, value in norm_dict1.items():
    sum_norm_scores1[key] = sum(value)
for key, value in norm_dict2.items():
    sum_norm_scores2[key] = sum(value)
for TF in TFs:
    avg_norm_dict1[TF] = []
for TF in TFs:
    avg_norm_dict2[TF] = []
for key, value in sum_norm_scores1.items():
    for TF in TFs:
        if TF in key:
            avg_norm_dict1[TF].append(value)
for key, value in sum_norm_scores2.items():
    for TF in TFs:
        if TF in key:
            avg_norm_dict2[TF].append(value)
out1 = [{k: sum(v) / len(v)} for k, v in avg_norm_dict1.items() if v != []]
out2 = [{k: sum(v) / len(v)} for k, v in avg_norm_dict2.items() if v != []]
# Separate the data for "kni0" and "kni2"
print(out1)
print(out2)
# Extract TF names and average values
tf_names_kni0 = [list(d.keys())[0] for d in out1]
tf_avg_values_kni0 = [list(d.values())[0] for d in out1]

tf_names_kni2 = [list(d.keys())[0] for d in out2]
tf_avg_values_kni2 = [list(d.values())[0] for d in out2]
# Find the intersection of TF names between out1 and out2
common_tfs = set(avg_norm_dict1.keys()) & set(avg_norm_dict2.keys())

# Extract TF names and average values for the common TFs
tf_names_common = TFs
tf_avg_values_kni0 = [list(out1[tf_names_kni0.index(tf)].values())[0] if tf in tf_names_kni0 else 0 for tf in tf_names_common]
tf_avg_values_kni2 = [list(out2[tf_names_kni2.index(tf)].values())[0] if tf in tf_names_kni2 else 0 for tf in tf_names_common]

# Create positions for the bars
bar_positions_common = np.arange(len(tf_names_common))



# Define the width of the bars
bar_width = 0.4  # Adjust this value as needed

# Create a horizontal bar plot
plt.barh(bar_positions_common, tf_avg_values_kni0, height=bar_width, color='orange', label='prox')
plt.barh(bar_positions_common + bar_width, tf_avg_values_kni2, height=bar_width, color='blue', label='distal')

# Set the y-axis ticks and labels
plt.yticks(bar_positions_common + bar_width / 2, tf_names_common)

plt.xlabel('Average Value')
plt.ylabel('TF')
plt.title('Average Values of sog enhancer TFs')
plt.legend()

# Label each bar with its value
for tf_name, avg_value, position in zip(tf_names_common, tf_avg_values_kni0, bar_positions_common):
    plt.text(avg_value, position, f'{avg_value:.2f}', ha='center', va='center', fontsize=8, color='black')

for tf_name, avg_value, position in zip(tf_names_common, tf_avg_values_kni2, bar_positions_common + bar_width):
    plt.text(avg_value, position, f'{avg_value:.2f}', ha='center', va='center', fontsize=8, color='black')

plt.tight_layout()
plt.show()