##plots curve but not happy with max and half max extraction

"""import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

# Define the construct list
construct_list = ['gtsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'  # Update this with your Dropbox folder path

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Two-component Gaussian function
def two_gaussian(x, a1, mu1, sigma1, a2, mu2, sigma2):
    return (a1 * np.exp(-((x - mu1) ** 2) / (2 * sigma1 ** 2)) +
            a2 * np.exp(-((x - mu2) ** 2) / (2 * sigma2 ** 2)))

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

                    if burst_properties is not None and len(burst_properties) > 0:
                        apbin_totals = defaultdict(lambda: [0, 0])  # [sum, count]

                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                if len(ap_bin_values) == 0 or len(total_mrna_values) == 0:
                                    print(f"Skipping {construct_path} due to empty data")
                                    continue

                                # Aggregate TotalmRNA by APBin
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna
                                    apbin_totals[ap_bin][1] += 1

                        # Ensure there is valid data to process
                        if len(apbin_totals) == 0:
                            print(f"No valid data for {construct_path}")
                            continue

                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0}
                        sorted_apbins = np.array(sorted(apbin_avg_mrna.keys()))
                        sorted_avg_mrna = np.array([apbin_avg_mrna[ap_bin] for ap_bin in sorted_apbins])

                        # Initial guesses for Gaussian parameters, with reasonable bounds
                        initial_guess = [max(sorted_avg_mrna), 0.2, 0.05, max(sorted_avg_mrna), 0.35, 0.05]
                        bounds = ([0, 0.15, 0, 0, 0.3, 0], [np.inf, 0.27, 0.1, np.inf, 0.4, 0.1])

                        # Fit the two-component Gaussian to the average mRNA trace
                        try:
                            popt, _ = curve_fit(two_gaussian, sorted_apbins, sorted_avg_mrna, p0=initial_guess, bounds=bounds)
                            fitted_curve = two_gaussian(sorted_apbins, *popt)

                            # Calculate half-max points from the fitted curve
                            half_max1 = popt[0] / 2  # Half of the first Gaussian's amplitude
                            half_max2 = popt[3] / 2  # Half of the second Gaussian's amplitude

                            # Find left half-max for the first Gaussian
                            left_indices = np.where((fitted_curve <= half_max1) & (sorted_apbins < popt[1]))[0]
                            if len(left_indices) > 0:
                                left_half_max1_bin = sorted_apbins[left_indices[-1]]
                                left_half_max1_mrna = fitted_curve[left_indices[-1]]
                            else:
                                left_half_max1_bin, left_half_max1_mrna = None, None

                            # Find right half-max for the second Gaussian
                            right_indices = np.where((fitted_curve <= half_max2) & (sorted_apbins > popt[4]))[0]
                            if len(right_indices) > 0:
                                right_half_max2_bin = sorted_apbins[right_indices[0]]
                                right_half_max2_mrna = fitted_curve[right_indices[0]]
                            else:
                                right_half_max2_bin, right_half_max2_mrna = None, None

                            # Store results
                            construct_data[construct].append({
                                'apbin_avg_mrna': apbin_avg_mrna,
                                'sorted_apbins': sorted_apbins,
                                'sorted_avg_mrna': sorted_avg_mrna,
                                'fitted_curve': fitted_curve,
                                'max1_apbin': popt[1], 'max1_mrna': popt[0],
                                'max2_apbin': popt[4], 'max2_mrna': popt[3],
                                'left_half_max1_bin': left_half_max1_bin, 'left_half_max1_mrna': left_half_max1_mrna,
                                'right_half_max2_bin': right_half_max2_bin, 'right_half_max2_mrna': right_half_max2_mrna
                            })

                        except RuntimeError:
                            print(f"Fit failed for construct {construct} in {construct_path}")

                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Plotting the original trace and fitted Gaussian for each embryo
for construct, data in construct_data.items():
    for i, d in enumerate(data):
        if len(d['sorted_apbins']) == 0 or len(d['fitted_curve']) == 0:
            print(f"Skipping plot for embryo {i+1} due to insufficient data")
            continue

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot original average mRNA trace
        ax.plot(d['sorted_apbins'], d['sorted_avg_mrna'], label='Avg Total mRNA (Original)', marker='o', linestyle='-', color='gray')

        # Plot fitted Gaussian curve
        ax.plot(d['sorted_apbins'], d['fitted_curve'], label='Fitted Gaussian', linestyle='--', color='blue')

        # Plot maxima and half-maxima from Gaussian fit
        if d['max1_apbin'] is not None:
            ax.plot(d['max1_apbin'], d['max1_mrna'], 'ro', label='Max 1')
        if d['left_half_max1_bin'] is not None:
            ax.plot(d['left_half_max1_bin'], d['left_half_max1_mrna'], 'rx', label='Left Half-Max 1')
        if d['max2_apbin'] is not None:
            ax.plot(d['max2_apbin'], d['max2_mrna'], 'bo', label='Max 2')
        if d['right_half_max2_bin'] is not None:
            ax.plot(d['right_half_max2_bin'], d['right_half_max2_mrna'], 'bx', label='Right Half-Max 2')

        # Set labels and title
        ax.set_xlabel('APBin')
        ax.set_ylabel('Total mRNA')
        ax.set_title(f'Embryo {i+1} of {construct}: Gaussian Fit and Half-Max Calculation')
        ax.grid(True)
        ax.legend()

        plt.show()"""


import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

# Define the construct list
construct_list = ['gtSE_5', 'gtsquish_0', 'gtsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Two-component Gaussian function
def two_gaussian(x, a1, mu1, sigma1, a2, mu2, sigma2):
    return (a1 * np.exp(-((x - mu1) ** 2) / (2 * sigma1 ** 2)) +
            a2 * np.exp(-((x - mu2) ** 2) / (2 * sigma2 ** 2)))

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

                    if burst_properties is not None and len(burst_properties) > 0:
                        apbin_totals = defaultdict(lambda: [0, 0])  # [sum, count]

                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                if len(ap_bin_values) == 0 or len(total_mrna_values) == 0:
                                    continue

                                # Aggregate TotalmRNA by APBin
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna
                                    apbin_totals[ap_bin][1] += 1

                        # Ensure there is valid data to process
                        if len(apbin_totals) == 0:
                            continue

                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0}
                        sorted_apbins = np.array(sorted(apbin_avg_mrna.keys()))
                        sorted_avg_mrna = np.array([apbin_avg_mrna[ap_bin] for ap_bin in sorted_apbins])

                        # Initial guesses for Gaussian parameters
                        initial_guess = [max(sorted_avg_mrna), 0.2, 0.05, max(sorted_avg_mrna), 0.35, 0.05]
                        bounds = ([0, 0.15, 0, 0, 0.3, 0], [np.inf, 0.27, 0.1, np.inf, 0.4, 0.1])

                        # Fit the two-component Gaussian to the average mRNA trace
                        try:
                            popt, _ = curve_fit(two_gaussian, sorted_apbins, sorted_avg_mrna, p0=initial_guess, bounds=bounds)
                            fitted_curve = two_gaussian(sorted_apbins, *popt)

                            # Find max1 within AP bin range 0.1 to 0.25
                            max1_range_indices = np.where((sorted_apbins >= 0.1) & (sorted_apbins <= 0.25))
                            max1_idx = max1_range_indices[0][np.argmax(fitted_curve[max1_range_indices])]
                            max1_apbin = sorted_apbins[max1_idx]
                            max1_mrna = fitted_curve[max1_idx]

                            # Calculate half-max1 (should be lower than max1 AP bin)
                            half_max1 = max1_mrna / 2
                            left_half_max1_idx = np.where((fitted_curve[:max1_idx] <= half_max1))[0]
                            if len(left_half_max1_idx) > 0:
                                left_half_max1_bin = sorted_apbins[left_half_max1_idx[-1]]
                                left_half_max1_mrna = fitted_curve[left_half_max1_idx[-1]]
                            else:
                                left_half_max1_bin, left_half_max1_mrna = None, None

                            # Find max2 within AP bin range 0.28 to 0.4
                            max2_range_indices = np.where((sorted_apbins >= 0.28) & (sorted_apbins <= 0.4))
                            max2_idx = max2_range_indices[0][np.argmax(fitted_curve[max2_range_indices])]
                            max2_apbin = sorted_apbins[max2_idx]
                            max2_mrna = fitted_curve[max2_idx]

                            # Calculate half-max2 (should be higher than max2 AP bin)
                            half_max2 = max2_mrna / 2
                            right_half_max2_idx = np.where((fitted_curve[max2_idx:] <= half_max2))[0]
                            if len(right_half_max2_idx) > 0:
                                right_half_max2_bin = sorted_apbins[max2_idx:][right_half_max2_idx[0]]
                                right_half_max2_mrna = fitted_curve[max2_idx:][right_half_max2_idx[0]]
                            else:
                                right_half_max2_bin, right_half_max2_mrna = None, None

                            # Store results
                            construct_data[construct].append({
                                'sorted_apbins': sorted_apbins,
                                'sorted_avg_mrna': sorted_avg_mrna,
                                'fitted_curve': fitted_curve,
                                'max1_apbin': max1_apbin, 'max1_mrna': max1_mrna,
                                'max2_apbin': max2_apbin, 'max2_mrna': max2_mrna,
                                'left_half_max1_bin': left_half_max1_bin, 'left_half_max1_mrna': left_half_max1_mrna,
                                'right_half_max2_bin': right_half_max2_bin, 'right_half_max2_mrna': right_half_max2_mrna
                            })

                        except RuntimeError:
                            print(f"Fit failed for construct {construct} in {construct_path}")

                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Plotting the original trace and fitted Gaussian for each embryo
for construct, data in construct_data.items():
    for i, d in enumerate(data):
        if len(d['sorted_apbins']) == 0 or len(d['fitted_curve']) == 0:
            continue

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot original average mRNA trace
        ax.plot(d['sorted_apbins'], d['sorted_avg_mrna'], label='Avg Total mRNA (Original)', marker='o', linestyle='-', color='gray')

        # Plot fitted Gaussian curve
        ax.plot(d['sorted_apbins'], d['fitted_curve'], label='Fitted Gaussian', linestyle='--', color='blue')

        # Plot maxima and half-maxima from Gaussian fit
        if d['max1_apbin'] is not None:
            ax.plot(d['max1_apbin'], d['max1_mrna'], 'ro', label='Max 1')
        if d['left_half_max1_bin'] is not None:
            ax.plot(d['left_half_max1_bin'], d['left_half_max1_mrna'], 'rx', label='Left Half-Max 1')
        if d['max2_apbin'] is not None:
            ax.plot(d['max2_apbin'], d['max2_mrna'], 'bo', label='Max 2')
        if d['right_half_max2_bin'] is not None:
            ax.plot(d['right_half_max2_bin'], d['right_half_max2_mrna'], 'bx', label='Right Half-Max 2')

        # Set labels and title
        ax.set_xlabel('APBin')
        ax.set_ylabel('Total mRNA')
        ax.set_title(f'Embryo {i+1} of {construct}: Gaussian Fit and Half-Max Calculation')
        ax.grid(True)
        ax.legend()

        plt.show()

