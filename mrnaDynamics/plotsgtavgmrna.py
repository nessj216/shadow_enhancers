import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from collections import defaultdict

# Define the construct list
construct_list = ['gtSE_5-m5']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    if isinstance(nested_array, (np.ndarray, list)) and nested_array.size > 0:
        return np.ravel(nested_array).astype(float)
    return np.array([])

# Dictionary to hold the interpolated data for each embryo
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
                            if 'APBin' in entry.dtype.names and 'TotalmRNA' in entry.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                if ap_bin_values.size == 0 or total_mrna_values.size == 0:
                                    continue  # Skip if no valid data remains

                                # Update totals
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna
                                    apbin_totals[ap_bin][1] += 1

                        # Calculate average TotalmRNA for each APBin
                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0}
                        ap_bins = sorted(apbin_avg_mrna.keys())
                        avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in ap_bins]

                        if len(ap_bins) > 0 and len(avg_mrna) > 0:
                            # Plot avg mRNA for valid data points
                            plt.figure(figsize=(10, 6))
                            plt.plot(ap_bins, avg_mrna, 'o-', label='Avg mRNA (non-NaN)')
                            plt.xlabel('APBin')
                            plt.ylabel('Total mRNA')
                            plt.title(f'{construct} - {dir_name}')
                            plt.grid(True)
                            plt.legend()
                            plt.show()

                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")
