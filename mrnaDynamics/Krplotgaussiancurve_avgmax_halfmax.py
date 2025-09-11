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
construct_list = ['KrSE_0', 'Krsquish_0', 'Krsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'  # Update this with your Dropbox folder path

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Dictionary to hold the interpolated half-max data for each embryo
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

                    # Process 'BurstProperties' if it exists
                    if burst_properties is not None:
                        # Dictionary to hold the sum and count of 'TotalmRNA' for each APBin within an embryo
                        apbin_totals = defaultdict(lambda: [0, 0])  # [sum, count]

                        for entry in burst_properties[0]:
                            # Extract 'APBin' and 'TotalmRNA' from each entry
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values from TotalmRNA
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                # Aggregate TotalmRNA by APBin for this embryo
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna  # Sum of TotalmRNA
                                    apbin_totals[ap_bin][1] += 1      # Count for averaging

                        # Calculate the average TotalmRNA for each APBin in the current embryo
                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0}
                        ap_bins = sorted(apbin_avg_mrna.keys())
                        avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in ap_bins]

                        # Fit a Gaussian curve to the average TotalmRNA over the AP length
                        try:
                            popt, _ = curve_fit(gaussian, ap_bins, avg_mrna, p0=[max(avg_mrna), np.mean(ap_bins), np.std(ap_bins)])
                            fitted_max = popt[1]
                            half_max_value = popt[0] / 2

                            # Interpolate the AP bins at half-max
                            interpolated_left_half_max = fitted_max - np.sqrt(2 * (np.log(popt[0] / half_max_value)) * popt[2]**2)
                            interpolated_right_half_max = fitted_max + np.sqrt(2 * (np.log(popt[0] / half_max_value)) * popt[2]**2)

                            construct_data[construct].append({
                                'ap_bins': ap_bins,
                                'avg_mrna': avg_mrna,
                                'popt': popt,
                                'interpolated_left_half_max': interpolated_left_half_max,
                                'interpolated_right_half_max': interpolated_right_half_max
                            })
                        except RuntimeError:
                            print(f"Gaussian fit failed for embryo: {dir_name}")
                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Plotting the results for each embryo
for construct, data in construct_data.items():
    for i, embryo_data in enumerate(data):
        ap_bins = embryo_data['ap_bins']
        avg_mrna = embryo_data['avg_mrna']
        popt = embryo_data['popt']

        # Generate Gaussian curve for fitted values
        fitted_curve = gaussian(np.array(ap_bins), *popt)

        # Plot the average mRNA trace and the fitted Gaussian curve
        plt.figure(figsize=(8, 6))
        plt.plot(ap_bins, avg_mrna, 'o-', label='Average TotalmRNA')
        plt.plot(ap_bins, fitted_curve, '-', label='Fitted Gaussian')
        plt.axhline(y=popt[0] / 2, color='r', linestyle='--', label='Half-Max')
        plt.xlabel('AP Bin')
        plt.ylabel('Average TotalmRNA')
        plt.title(f'Embryo {i + 1} - {construct}: Avg TotalmRNA and Gaussian Fit')
        plt.legend()
        plt.grid(True)
        plt.show()

