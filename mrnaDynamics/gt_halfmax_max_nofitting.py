import os
import numpy as np
import scipy.io
from collections import defaultdict
import matplotlib.pyplot as plt

# Define the construct list
construct_list = ['gtsquish_32C']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

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

                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0}

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

                        # Calculate left half-max for range 0.15-0.27
                        left1_half_max_bin, left1_half_max_mrna = None, None
                        if max1_avg_mrna is not None:
                            half_max1 = max1_avg_mrna / 2
                            # Calculate the absolute difference from the half-max for each AP bin less than max1_apbin
                            left_candidates = [(ap_bin, abs(apbin_avg_mrna[ap_bin] - half_max1))
                                               for ap_bin in sorted(apbin_avg_mrna.keys())
                                               if ap_bin < max1_apbin]
                            if left_candidates:
                                left1_half_max_bin, _ = min(left_candidates, key=lambda x: x[1])
                                left1_half_max_mrna = apbin_avg_mrna[left1_half_max_bin]

                        # Calculate right half-max for range 0.3-0.4
                        right2_half_max_bin, right2_half_max_mrna = None, None
                        if max2_avg_mrna is not None:
                            half_max2 = max2_avg_mrna / 2
                            # Calculate the absolute difference from the half-max for each AP bin greater than max2_apbin
                            right_candidates = [(ap_bin, abs(apbin_avg_mrna[ap_bin] - half_max2))
                                                for ap_bin in sorted(apbin_avg_mrna.keys())
                                                if ap_bin > max2_apbin]
                            if right_candidates:
                                right2_half_max_bin, _ = min(right_candidates, key=lambda x: x[1])
                                right2_half_max_mrna = apbin_avg_mrna[right2_half_max_bin]

                        # Store results
                        construct_data[construct].append({
                            'apbin_avg_mrna': apbin_avg_mrna,
                            'max1_apbin': max1_apbin, 'max1_avg_mrna': max1_avg_mrna,
                            'left1_half_max_bin': left1_half_max_bin, 'left1_half_max_mrna': left1_half_max_mrna,
                            'max2_apbin': max2_apbin, 'max2_avg_mrna': max2_avg_mrna,
                            'right2_half_max_bin': right2_half_max_bin, 'right2_half_max_mrna': right2_half_max_mrna
                        })
                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Plotting the results
for construct, data in construct_data.items():
    for i, d in enumerate(data):
        apbins = list(d['apbin_avg_mrna'].keys())
        avg_mrna = list(d['apbin_avg_mrna'].values())

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(apbins, avg_mrna, label='Avg Total mRNA', marker='o', linestyle='-')

        # Plot max and half-max points for range 0.15-0.27
        if d['max1_apbin'] is not None:
            ax.plot(d['max1_apbin'], d['max1_avg_mrna'], 'ro', label='Max 1')
            if d['left1_half_max_bin'] is not None:
                ax.plot(d['left1_half_max_bin'], d['left1_half_max_mrna'], 'rx', label='Left Half-Max 1')

        # Plot max and half-max points for range 0.3-0.4
        if d['max2_apbin'] is not None:
            ax.plot(d['max2_apbin'], d['max2_avg_mrna'], 'bo', label='Max 2')
            if d['right2_half_max_bin'] is not None:
                ax.plot(d['right2_half_max_bin'], d['right2_half_max_mrna'], 'bx', label='Right Half-Max 2')

        ax.set_xlabel('APBin')
        ax.set_ylabel('Total mRNA')
        ax.set_title(f'Embryo {i + 1} of {construct}: Max and Half-Max')
        ax.grid(True)
        ax.legend()

        plt.show()