import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict
import matplotlib.cm as cm

# Define the Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Define the construct list
construct_list = ['KrSE_0', 'Krsquish_0', 'Krsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Dictionary to hold the data for each construct
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
                            fitted_curve = gaussian(np.array(ap_bins), *popt)

                            # Store data for plotting, including the directory name
                            construct_data[construct].append({
                                'ap_bins': ap_bins,
                                'avg_mrna': avg_mrna,
                                'fitted_curve': fitted_curve,
                                'dir_name': dir_name  # Add directory name for labeling
                            })
                        except RuntimeError:
                            print(f"Gaussian fit failed for embryo: {dir_name}")
                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Plotting mRNA curves and Gaussian fits for each construct
for construct, data in construct_data.items():
    fig, ax = plt.subplots(figsize=(10, 6))

    # Create a colormap for unique colors
    colors = cm.viridis(np.linspace(0, 1, len(data)))

    # Plot each embryo's mRNA data and Gaussian fit with a unique color
    for idx, embryo_data in enumerate(data):
        ap_bins = embryo_data['ap_bins']
        avg_mrna = embryo_data['avg_mrna']
        fitted_curve = embryo_data['fitted_curve']
        dir_name = embryo_data['dir_name']

        # Plot the mRNA curve with a unique color
        ax.plot(ap_bins, avg_mrna, 'o-', label=f'{dir_name} mRNA Data', color=colors[idx])

        # Plot the Gaussian fit with the same color
        ax.plot(ap_bins, fitted_curve, '-', label=f'{dir_name} Gaussian Fit', color=colors[idx])

    # Set plot details
    ax.set_xlabel('APBin')
    ax.set_ylabel('TotalmRNA')
    ax.set_title(f'{construct}: mRNA Curves and Gaussian Fits')
    ax.set_ylim(bottom=0)  # Set y-axis to start at 0
    ax.grid(True)
    ax.legend(loc='upper right', fontsize='small')

    plt.show()
