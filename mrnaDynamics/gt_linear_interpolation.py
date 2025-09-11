import os
import numpy as np
import scipy.io
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.stats import sem  # For standard error calculation

# Define the construct list
construct_list = ['gtSE_5', 'gtsquish_0', 'gtsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Specify colors for each construct
construct_colors = {
    'gtSE_5': 'yellowgreen',
    'gtsquish_0': 'plum',
    'gtsquish_mini': 'slategrey'
}

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Forced linear interpolation function
def force_linear_interpolate(sorted_x, sorted_y, target_y):
    # Find the indices where the target_y falls between two data points
    for i in range(1, len(sorted_y)):
        if (sorted_y[i - 1] <= target_y <= sorted_y[i]) or (sorted_y[i - 1] >= target_y >= sorted_y[i]):
            x1, y1 = sorted_x[i - 1], sorted_y[i - 1]
            x2, y2 = sorted_x[i], sorted_y[i]
            return x1 + (target_y - y1) * (x2 - x1) / (y2 - y1)
    return None

# Dictionary to hold data for each embryo
construct_data = defaultdict(list)

# Iterate through the construct list and load the BurstPropertiesSlope.mat file for each
for construct in construct_list:
    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            construct_path = os.path.join(dropbox_folder, dir_name, 'BurstPropertiesSlope.mat')

            if os.path.exists(construct_path):
                try:
                    # Load the .mat file
                    mat_data = scipy.io.loadmat(construct_path)
                    burst_properties = mat_data.get('BurstProperties', None)

                    if burst_properties is not None:
                        apbin_totals = defaultdict(lambda: [0, 0])  # [sum, count]

                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                # Aggregate TotalmRNA by APBin
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna
                                    apbin_totals[ap_bin][1] += 1

                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if
                                          total[1] > 0}

                        # Define the ranges
                        range1_min, range1_max = 0.15, 0.27
                        range2_min, range2_max = 0.3, 0.4

                        # Find max in the range 0.15-0.27
                        max1_apbin = max((ap_bin for ap_bin in apbin_avg_mrna if range1_min <= ap_bin <= range1_max),
                                         key=apbin_avg_mrna.get, default=None)
                        max1_avg_mrna = apbin_avg_mrna[max1_apbin] if max1_apbin is not None else None

                        # Find max in the range 0.3-0.4
                        max2_apbin = max((ap_bin for ap_bin in apbin_avg_mrna if range2_min <= ap_bin <= range2_max),
                                         key=apbin_avg_mrna.get, default=None)
                        max2_avg_mrna = apbin_avg_mrna[max2_apbin] if max2_apbin is not None else None

                        # Calculate left half-max for range 0.15-0.27 with forced interpolation and restriction
                        left1_half_max_bin, left1_half_max_mrna = None, None
                        if max1_avg_mrna is not None:
                            half_max1 = max1_avg_mrna / 2
                            sorted_apbins = sorted(apbin_avg_mrna.keys())
                            sorted_avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in sorted_apbins]

                            # Restrict to AP bins less than max1_apbin for left half-max
                            left_apbins = [ap_bin for ap_bin in sorted_apbins if ap_bin < max1_apbin]
                            left_avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in left_apbins]

                            # Force interpolation if half-max is within the restricted range
                            left1_half_max_bin = force_linear_interpolate(left_apbins, left_avg_mrna, half_max1)
                            left1_half_max_mrna = half_max1 if left1_half_max_bin is not None else None

                        # Calculate right half-max for range 0.3-0.4 with forced interpolation and restriction
                        right2_half_max_bin, right2_half_max_mrna = None, None
                        if max2_avg_mrna is not None:
                            half_max2 = max2_avg_mrna / 2
                            sorted_apbins = sorted(apbin_avg_mrna.keys())
                            sorted_avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in sorted_apbins]

                            # Restrict to AP bins greater than max2_apbin for right half-max
                            right_apbins = [ap_bin for ap_bin in sorted_apbins if ap_bin > max2_apbin]
                            right_avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in right_apbins]

                            # Force interpolation if half-max is within the restricted range
                            right2_half_max_bin = force_linear_interpolate(right_apbins, right_avg_mrna, half_max2)
                            right2_half_max_mrna = half_max2 if right2_half_max_bin is not None else None

                        # Store results
                        construct_data[construct].append({
                            'max1_apbin': max1_apbin, 'max1_avg_mrna': max1_avg_mrna,
                            'left1_half_max_bin': left1_half_max_bin, 'left1_half_max_mrna': left1_half_max_mrna,
                            'max2_apbin': max2_apbin, 'max2_avg_mrna': max2_avg_mrna,
                            'right2_half_max_bin': right2_half_max_bin, 'right2_half_max_mrna': right2_half_max_mrna
                        })
                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Calculate and plot averages and standard errors for each construct
fig, ax = plt.subplots(figsize=(10, 6))

for construct, data in construct_data.items():
    max1_apbins = [d['max1_apbin'] for d in data if d['max1_apbin'] is not None]
    max1_mrna = [d['max1_avg_mrna'] for d in data if d['max1_avg_mrna'] is not None]
    left1_half_bins = [d['left1_half_max_bin'] for d in data if d['left1_half_max_bin'] is not None]
    left1_half_mrna = [d['left1_half_max_mrna'] for d in data if d['left1_half_max_mrna'] is not None]
    max2_apbins = [d['max2_apbin'] for d in data if d['max2_apbin'] is not None]
    max2_mrna = [d['max2_avg_mrna'] for d in data if d['max2_avg_mrna'] is not None]
    right2_half_bins = [d['right2_half_max_bin'] for d in data if d['right2_half_max_bin'] is not None]
    right2_half_mrna = [d['right2_half_max_mrna'] for d in data if d['right2_half_max_mrna'] is not None]

    avg_apbins = [np.mean(max1_apbins), np.mean(left1_half_bins), np.mean(max2_apbins), np.mean(right2_half_bins)]
    avg_mrna = [np.mean(max1_mrna), np.mean(left1_half_mrna), np.mean(max2_mrna), np.mean(right2_half_mrna)]
    sem_apbins = [sem(max1_apbins), sem(left1_half_bins), sem(max2_apbins), sem(right2_half_bins)]
    sem_mrna = [sem(max1_mrna), sem(left1_half_mrna), sem(max2_mrna), sem(right2_half_mrna)]

    # Plot average maxes and half-maxes with error bars for each construct
    ax.errorbar(
        avg_apbins, avg_mrna, xerr=sem_apbins, yerr=sem_mrna, fmt='o',
        capsize=5, label=construct, color=construct_colors[construct]
    )

ax.set_xlabel('AP Bin')
ax.set_ylabel('Total mRNA')
ax.set_title('Average Maxima and Half-Maxima for All Constructs')
ax.grid(True)
ax.legend()

plt.show()

import os
import numpy as np
import scipy.io
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.stats import sem  # For standard error calculation

# Define the construct list
construct_list = ['gtSE_5', 'gtsquish_0', 'gtsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Specify colors for each construct
construct_colors = {
    'gtSE_5': 'yellowgreen',
    'gtsquish_0': 'plum',
    'gtsquish_mini': 'slategrey'
}


# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)


# Forced linear interpolation function
def force_linear_interpolate(sorted_x, sorted_y, target_y):
    for i in range(1, len(sorted_y)):
        if (sorted_y[i - 1] <= target_y <= sorted_y[i]) or (sorted_y[i - 1] >= target_y >= sorted_y[i]):
            x1, y1 = sorted_x[i - 1], sorted_y[i - 1]
            x2, y2 = sorted_x[i], sorted_y[i]
            return x1 + (target_y - y1) * (x2 - x1) / (y2 - y1)
    return None


# Dictionary to hold data for each construct
construct_data = defaultdict(list)

# Iterate through the construct list and load the BurstPropertiesSlope.mat file for each
for construct in construct_list:
    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            construct_path = os.path.join(dropbox_folder, dir_name, 'BurstPropertiesSlope.mat')

            if os.path.exists(construct_path):
                try:
                    # Load the .mat file
                    mat_data = scipy.io.loadmat(construct_path)
                    burst_properties = mat_data.get('BurstProperties', None)

                    if burst_properties is not None:
                        apbin_totals = defaultdict(lambda: [0, 0])  # [sum, count]

                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                # Aggregate TotalmRNA by APBin
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna
                                    apbin_totals[ap_bin][1] += 1

                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if
                                          total[1] > 0}

                        # Store average total mRNA trace
                        sorted_apbins = sorted(apbin_avg_mrna.keys())
                        sorted_avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in sorted_apbins]

                        construct_data[construct].append({
                            'apbins': sorted_apbins,
                            'avg_mrna': sorted_avg_mrna
                        })

                        # Find maxes and half-maxes as before (same code for ranges and interpolation)

                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Plot averages, maxes, half-maxes, and total mRNA traces
fig, ax = plt.subplots(figsize=(10, 6))

for construct, data in construct_data.items():
    # Extract average total mRNA trace
    for trace_data in data:
        apbins = trace_data['apbins']
        avg_mrna_trace = trace_data['avg_mrna']

        # Plot total mRNA trace in background
        ax.plot(apbins, avg_mrna_trace, label=f'{construct} Trace', color=construct_colors[construct], alpha=0.5)

    max1_apbins = [d.get('max1_apbin') for d in data if 'max1_apbin' in d]
    max1_mrna = [d.get('max1_avg_mrna') for d in data if 'max1_avg_mrna' in d]
    left1_half_bins = [d.get('left1_half_max_bin') for d in data if 'left1_half_max_bin' in d]
    left1_half_mrna = [d.get('left1_half_max_mrna') for d in data if 'left1_half_max_mrna' in d]
    max2_apbins = [d.get('max2_apbin') for d in data if 'max2_apbin' in d]
    max2_mrna = [d.get('max2_avg_mrna') for d in data if 'max2_avg_mrna' in d]
    right2_half_bins = [d.get('right2_half_max_bin') for d in data if 'right2_half_max_bin' in d]
    right2_half_mrna = [d.get('right2_half_max_mrna') for d in data if 'right2_half_max_mrna' in d]

    # Plot maxes and half-maxes
    ax.scatter(max1_apbins, max1_mrna, label=f'{construct} Max1', color=construct_colors[construct])
    ax.scatter(left1_half_bins, left1_half_mrna, label=f'{construct} Left1 Half-Max', color=construct_colors[construct],
               marker='x')
    ax.scatter(max2_apbins, max2_mrna, label=f'{construct} Max2', color=construct_colors[construct])
    ax.scatter(right2_half_bins, right2_half_mrna, label=f'{construct} Right2 Half-Max',
               color=construct_colors[construct], marker='x')

ax.set_xlabel('AP Bin')
ax.set_ylabel('Total mRNA')
ax.set_title('Average Total mRNA Traces, Maxes, and Half-Maxima for All Constructs')
ax.grid(True)
#ax.legend()

plt.show()

import os
import numpy as np
import scipy.io
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.stats import sem  # For standard error calculation

# Define the construct list
construct_list = ['gtSE_5', 'gtsquish_0', 'gtsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Specify colors for each construct
construct_colors = {
    'gtSE_5': 'yellowgreen',
    'gtsquish_0': 'plum',
    'gtsquish_mini': 'slategrey'
}

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Forced linear interpolation function
def force_linear_interpolate(sorted_x, sorted_y, target_y):
    for i in range(1, len(sorted_y)):
        if (sorted_y[i - 1] <= target_y <= sorted_y[i]) or (sorted_y[i - 1] >= target_y >= sorted_y[i]):
            x1, y1 = sorted_x[i - 1], sorted_y[i - 1]
            x2, y2 = sorted_x[i], sorted_y[i]
            return x1 + (target_y - y1) * (x2 - x1) / (y2 - y1)
    return None

# Dictionary to hold aggregated total mRNA data for each construct
construct_totals = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # {construct: {ap_bin: [sum, count]}}

# Iterate through the construct list and load the BurstPropertiesSlope.mat file for each
for construct in construct_list:
    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            construct_path = os.path.join(dropbox_folder, dir_name, 'BurstPropertiesSlope.mat')

            if os.path.exists(construct_path):
                try:
                    # Load the .mat file
                    mat_data = scipy.io.loadmat(construct_path)
                    burst_properties = mat_data.get('BurstProperties', None)

                    if burst_properties is not None:
                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                # Aggregate TotalmRNA by APBin for averaging across all embryos
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    construct_totals[construct][ap_bin][0] += mrna
                                    construct_totals[construct][ap_bin][1] += 1

                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Calculate average total mRNA across all embryos for each construct
avg_mrna_traces = {}

for construct, apbin_totals in construct_totals.items():
    avg_mrna_traces[construct] = {
        ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0
    }

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

for construct, avg_mrna_trace in avg_mrna_traces.items():
    # Sort AP bins for plotting
    sorted_apbins = sorted(avg_mrna_trace.keys())
    sorted_avg_mrna = [avg_mrna_trace[ap_bin] for ap_bin in sorted_apbins]

    # Plot average total mRNA trace for each construct
    ax.plot(sorted_apbins, sorted_avg_mrna, label=f'{construct} Avg Trace', color=construct_colors[construct], alpha=0.7)

    # Find maxes and half-maxes
    range1_min, range1_max = 0.15, 0.27
    range2_min, range2_max = 0.3, 0.4

    # Max1 and Max2
    max1_apbin = max((ap_bin for ap_bin in avg_mrna_trace if range1_min <= ap_bin <= range1_max),
                     key=avg_mrna_trace.get, default=None)
    max1_avg_mrna = avg_mrna_trace[max1_apbin] if max1_apbin is not None else None

    max2_apbin = max((ap_bin for ap_bin in avg_mrna_trace if range2_min <= ap_bin <= range2_max),
                     key=avg_mrna_trace.get, default=None)
    max2_avg_mrna = avg_mrna_trace[max2_apbin] if max2_apbin is not None else None

    # Left half-max for range1
    left1_half_max_bin, left1_half_max_mrna = None, None
    if max1_avg_mrna is not None:
        half_max1 = max1_avg_mrna / 2
        left_apbins = [ap_bin for ap_bin in sorted_apbins if ap_bin < max1_apbin]
        left_avg_mrna = [avg_mrna_trace[ap_bin] for ap_bin in left_apbins]

        # Force interpolation if half-max is within the restricted range
        left1_half_max_bin = force_linear_interpolate(left_apbins, left_avg_mrna, half_max1)
        left1_half_max_mrna = half_max1 if left1_half_max_bin is not None else None

    # Right half-max for range2
    right2_half_max_bin, right2_half_max_mrna = None, None
    if max2_avg_mrna is not None:
        half_max2 = max2_avg_mrna / 2
        right_apbins = [ap_bin for ap_bin in sorted_apbins if ap_bin > max2_apbin]
        right_avg_mrna = [avg_mrna_trace[ap_bin] for ap_bin in right_apbins]

        # Force interpolation if half-max is within the restricted range
        right2_half_max_bin = force_linear_interpolate(right_apbins, right_avg_mrna, half_max2)
        right2_half_max_mrna = half_max2 if right2_half_max_bin is not None else None

    # Plot maxes and half-maxes
    if max1_apbin is not None:
        ax.scatter(max1_apbin, max1_avg_mrna, color=construct_colors[construct], label=f'{construct} Max1')
    if left1_half_max_bin is not None:
        ax.scatter(left1_half_max_bin, left1_half_max_mrna, color=construct_colors[construct], marker='x', label=f'{construct} Left1 Half-Max')
    if max2_apbin is not None:
        ax.scatter(max2_apbin, max2_avg_mrna, color=construct_colors[construct], label=f'{construct} Max2')
    if right2_half_max_bin is not None:
        ax.scatter(right2_half_max_bin, right2_half_max_mrna, color=construct_colors[construct], marker='x', label=f'{construct} Right2 Half-Max')

ax.set_xlabel('AP Bin')
ax.set_ylabel('Average Total mRNA')
ax.set_title('Average Total mRNA Traces, Maxes, and Half-Maxima Across All Embryos')
ax.grid(True)
#ax.legend()

#plt.show()





