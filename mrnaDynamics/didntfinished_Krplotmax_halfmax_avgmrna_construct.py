import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

# Define the Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Define the construct list
construct_list = ['gtsquish_0']

# Set the path to the Dropbox folder
dropbox_folder = '/Users/jillianness/Downloads/livemRNA/Data/DynamicsResults'

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Helper function to find max and half-max within a specified AP bin range
def find_max_and_half_max(ap_bins, avg_mrna, min_range, max_range):
    # Filter AP bins based on the given floating-point range
    ap_bins_filtered = [ap_bins[i] for i in range(len(ap_bins)) if min_range <= ap_bins[i] <= max_range]
    avg_mrna_filtered = [avg_mrna[i] for i in range(len(ap_bins)) if min_range <= ap_bins[i] <= max_range]

    if not ap_bins_filtered or not avg_mrna_filtered:
        print(f"No data found in AP range {min_range} to {max_range}")
        return None, None, None, None

    try:
        popt, _ = curve_fit(gaussian, ap_bins_filtered, avg_mrna_filtered,
                            p0=[max(avg_mrna_filtered), np.mean(ap_bins_filtered), np.std(ap_bins_filtered)])
        fitted_max = popt[1]
        half_max_value = popt[0] / 2

        # Calculate the interpolated half-max
        interpolated_left_half_max = fitted_max - np.sqrt(2 * (np.log(popt[0] / half_max_value)) * popt[2]**2)
        interpolated_right_half_max = fitted_max + np.sqrt(2 * (np.log(popt[0] / half_max_value)) * popt[2]**2)

        return fitted_max, interpolated_left_half_max, interpolated_right_half_max, half_max_value
    except RuntimeError:
        print("Gaussian fit failed.")
        return None, None, None, None

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

                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna
                                    apbin_totals[ap_bin][1] += 1

                        # Calculate average TotalmRNA for each APBin
                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0}
                        ap_bins = sorted(apbin_avg_mrna.keys())
                        avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in ap_bins]

                        # Debug: Check ap_bins and avg_mrna
                        print(f"AP Bins for {construct}: {ap_bins}")
                        print(f"Avg mRNA for {construct}: {avg_mrna}")

                        # Find max and half-max for AP bins within 0.0-0.27 (left of max)
                        left_results = find_max_and_half_max(ap_bins, avg_mrna, 0.0, 0.27)
                        left_max, left_half_max, _, left_half_max_value = left_results

                        # Debug: Check left_results
                        print(f"Left results for {construct}: {left_results}")

                        # Find max and half-max for AP bins within 0.29-0.40 (right of max)
                        right_results = find_max_and_half_max(ap_bins, avg_mrna, 0.29, 0.40)
                        right_max, _, right_half_max, right_half_max_value = right_results

                        # Debug: Check right_results
                        print(f"Right results for {construct}: {right_results}")

                        # Check if the returned results are valid
                        if all(v is not None for v in left_results) and all(v is not None for v in right_results):
                            construct_data[construct].append({
                                'left_max': left_max,
                                'left_half_max': left_half_max,
                                'left_half_max_value': left_half_max_value,
                                'right_max': right_max,
                                'right_half_max': right_half_max,
                                'right_half_max_value': right_half_max_value
                            })
                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Prepare data for calculating averages and SE
construct_means = {}

for construct, data in construct_data.items():
    left_maxs = [d['left_max'] for d in data if d['left_max'] is not None]
    left_half_maxs = [d['left_half_max'] for d in data if d['left_half_max'] is not None]
    left_half_max_values = [d['left_half_max_value'] for d in data if d['left_half_max_value'] is not None]
    right_maxs = [d['right_max'] for d in data if d['right_max'] is not None]
    right_half_maxs = [d['right_half_max'] for d in data if d['right_half_max'] is not None]
    right_half_max_values = [d['right_half_max_value'] for d in data if d['right_half_max_value'] is not None]

    # Calculate mean and SE for left and right max and half-max
    construct_means[construct] = {
        'mean_left_max': np.mean(left_maxs) if left_maxs else None,
        'se_left_max': np.std(left_maxs) / np.sqrt(len(left_maxs)) if left_maxs else None,
        'mean_left_half_max': np.mean(left_half_maxs) if left_half_maxs else None,
        'se_left_half_max': np.std(left_half_maxs) / np.sqrt(len(left_half_maxs)) if left_half_maxs else None,
        'mean_left_half_max_value': np.mean(left_half_max_values) if left_half_max_values else None,
        'se_left_half_max_value': np.std(left_half_max_values) / np.sqrt(len(left_half_max_values)) if left_half_max_values else None,
        'mean_right_max': np.mean(right_maxs) if right_maxs else None,
        'se_right_max': np.std(right_maxs) / np.sqrt(len(right_maxs)) if right_maxs else None,
        'mean_right_half_max': np.mean(right_half_maxs) if right_half_maxs else None,
        'se_right_half_max': np.std(right_half_maxs) / np.sqrt(len(right_half_maxs)) if right_half_maxs else None,
        'mean_right_half_max_value': np.mean(right_half_max_values) if right_half_max_values else None,
        'se_right_half_max_value': np.std(right_half_max_values) / np.sqrt(len(right_half_max_values)) if right_half_max_values else None,
    }

# Debug: Check construct_means
print(f"Construct means: {construct_means}")

# Plotting the results
fig, ax = plt.subplots(figsize=(12, 8))

for construct, mean_data in construct_means.items():
    # Plot left max with error bars
    if mean_data['mean_left_max'] is not None:
        ax.errorbar(
            mean_data['mean_left_max'], mean_data['mean_left_half_max_value'],
            xerr=mean_data['se_left_max'], yerr=mean_data['se_left_half_max_value'],
            fmt='o', capsize=5, label=f'{construct} Left Max'
        )

    # Plot left half-max with error bars
    if mean_data['mean_left_half_max'] is not None:
        ax.errorbar(
            mean_data['mean_left_half_max'], mean_data['mean_left_half_max_value'],
            xerr=mean_data['se_left_half_max'], yerr=mean_data['se_left_half_max_value'],
            fmt='o', capsize=5, label=f'{construct} Left Half-Max'
        )

    # Plot right max with error bars
    if mean_data['mean_right_max'] is not None:
        ax.errorbar(
            mean_data['mean_right_max'], mean_data['mean_right_half_max_value'],
            xerr=mean_data['se_right_max'], yerr=mean_data['se_right_half_max_value'],
            fmt='o', capsize=5, label=f'{construct} Right Max'
        )

    # Plot right half-max with error bars
    if mean_data['mean_right_half_max'] is not None:
        ax.errorbar(
            mean_data['mean_right_half_max'], mean_data['mean_right_half_max_value'],
            xerr=mean_data['se_right_half_max'], yerr=mean_data['se_right_half_max_value'],
            fmt='o', capsize=5, label=f'{construct} Right Half-Max'
        )

# Final plot settings
ax.set_xlabel('APBin')
ax.set_ylabel('TotalmRNA')
ax.set_title('Max and Half-Max TotalmRNA for Left and Right AP Bins Across Embryos')
ax.grid(True)
ax.legend()

plt.show()
