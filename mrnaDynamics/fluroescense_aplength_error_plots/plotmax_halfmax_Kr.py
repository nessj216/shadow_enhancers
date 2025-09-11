import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

# Define the Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Define the construct list and their corresponding colors
construct_list = ['Krsquish_mini', 'Kr_SEmini_endog']
#construct_list = ['KrSE_0', 'Krsquish_0','Krsquish_mini']
construct_colors = {
    'KrSE_0': 'yellowgreen',
    'Krsquish_0': 'plum',
    'Krsquish_mini': 'slategrey',
    'Kr_SEmini_endog': 'darkgreen'
}
#'yellowgreen', 'plum', 'darkgreen', 'purple'
# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

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
                                'max_apbin': fitted_max,
                                'max_avg_mrna': popt[0],
                                'interpolated_left_half_max': interpolated_left_half_max,
                                'interpolated_right_half_max': interpolated_right_half_max,
                                'half_max_value': half_max_value
                            })
                        except RuntimeError:
                            print(f"Gaussian fit failed for embryo: {dir_name}")
                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Prepare data for calculating averages and SE
construct_means = {}

for construct, data in construct_data.items():
    max_apbins = [d['max_apbin'] for d in data]
    max_avg_mrnas = [d['max_avg_mrna'] for d in data]
    left_half_bins = [d['interpolated_left_half_max'] for d in data]
    right_half_bins = [d['interpolated_right_half_max'] for d in data]
    half_max_values = [d['half_max_value'] for d in data]

    # Calculate mean and SE for max APBin and TotalmRNA
    mean_max_apbin = np.mean(max_apbins)
    se_max_apbin = np.std(max_apbins) / np.sqrt(len(max_apbins))

    mean_max_avg_mrna = np.mean(max_avg_mrnas)
    se_max_avg_mrna = np.std(max_avg_mrnas) / np.sqrt(len(max_avg_mrnas))

    # Calculate mean and SE for left half-max APBin and TotalmRNA
    mean_left_half_bin = np.mean(left_half_bins)
    se_left_half_bin = np.std(left_half_bins) / np.sqrt(len(left_half_bins))

    mean_half_max_value = np.mean(half_max_values)
    se_half_max_value = np.std(half_max_values) / np.sqrt(len(half_max_values))

    # Calculate mean and SE for right half-max APBin
    mean_right_half_bin = np.mean(right_half_bins)
    se_right_half_bin = np.std(right_half_bins) / np.sqrt(len(right_half_bins))

    construct_means[construct] = {
        'mean_max_apbin': mean_max_apbin,
        'mean_max_avg_mrna': mean_max_avg_mrna,
        'se_max_apbin': se_max_apbin,
        'se_max_avg_mrna': se_max_avg_mrna,
        'mean_left_half_bin': mean_left_half_bin,
        'se_left_half_bin': se_left_half_bin,
        'mean_right_half_bin': mean_right_half_bin,
        'se_right_half_bin': se_right_half_bin,
        'mean_half_max_value': mean_half_max_value,
        'se_half_max_value': se_half_max_value
    }

# Plotting the results
fig, ax = plt.subplots(figsize=(5, 3))

for construct, mean_data in construct_means.items():
    color = construct_colors.get(construct, 'black')  # Default to black if construct not found

    # Plot max APBin with error bars
    ax.errorbar(
        mean_data['mean_max_apbin'], mean_data['mean_max_avg_mrna'],
        xerr=mean_data['se_max_apbin'], yerr=mean_data['se_max_avg_mrna'],
        fmt='o', capsize=5, color=color, label=f'{construct} Max'
    )

    # Plot left half-max APBin with error bars
    ax.errorbar(
        mean_data['mean_left_half_bin'], mean_data['mean_half_max_value'],
        xerr=mean_data['se_left_half_bin'], yerr=mean_data['se_half_max_value'],
        fmt='o', capsize=5, color=color, label=f'{construct} Left Half-Max'
    )

    # Plot right half-max APBin with error bars
    ax.errorbar(
        mean_data['mean_right_half_bin'], mean_data['mean_half_max_value'],
        xerr=mean_data['se_right_half_bin'], yerr=mean_data['se_half_max_value'],
        fmt='o', capsize=5, color=color, label=f'{construct} Right Half-Max'
    )

ax.set_xlabel('% egg length',fontsize=14)
ax.set_ylabel('Integrated fluorescence (AU)',fontsize=14)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
#ax.set_title('Average Max and Half-Max Fluorescsence',fontsize=15)
ax.grid(True)
#ax.legend()
#ax.set_ylim(bottom=0)
ax.set_xlim(.35,.65)
plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/Krembryoaperror.png', dpi=300, bbox_inches='tight')

plt.show()



#









######next script
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

# Define the Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Define the construct list and their corresponding colors
construct_list = ['KrSE_0', 'Krsquish_0', 'Krsquish_mini']
construct_colors = {
    'KrSE_0': 'yellowgreen',
    'Krsquish_0': 'plum',
    'Krsquish_mini': 'slategrey'
}

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Dictionary to hold the interpolated half-max data for each embryo
construct_data = defaultdict(list)
construct_avg_traces = defaultdict(lambda: defaultdict(list))

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
                        apbin_totals = defaultdict(lambda: [0, 0])  # [sum, count]

                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    apbin_totals[ap_bin][0] += mrna  # Sum of TotalmRNA
                                    apbin_totals[ap_bin][1] += 1      # Count for averaging

                        apbin_avg_mrna = {ap_bin: total[0] / total[1] for ap_bin, total in apbin_totals.items() if total[1] > 0}
                        ap_bins = sorted(apbin_avg_mrna.keys())
                        avg_mrna = [apbin_avg_mrna[ap_bin] for ap_bin in ap_bins]

                        for ap_bin, avg in apbin_avg_mrna.items():
                            construct_avg_traces[construct][ap_bin].append(avg)

                        try:
                            popt, _ = curve_fit(gaussian, ap_bins, avg_mrna, p0=[max(avg_mrna), np.mean(ap_bins), np.std(ap_bins)])
                            fitted_max = popt[1]
                            half_max_value = popt[0] / 2

                            interpolated_left_half_max = fitted_max - np.sqrt(2 * (np.log(popt[0] / half_max_value)) * popt[2]**2)
                            interpolated_right_half_max = fitted_max + np.sqrt(2 * (np.log(popt[0] / half_max_value)) * popt[2]**2)

                            construct_data[construct].append({
                                'max_apbin': fitted_max,
                                'max_avg_mrna': popt[0],
                                'interpolated_left_half_max': interpolated_left_half_max,
                                'interpolated_right_half_max': interpolated_right_half_max,
                                'half_max_value': half_max_value
                            })
                        except RuntimeError:
                            print(f"Gaussian fit failed for embryo: {dir_name}")
                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Calculate average traces and 95% CI for each construct
for construct, ap_data in construct_avg_traces.items():
    for ap_bin, mrna_values in ap_data.items():
        mean_mrna = np.mean(mrna_values)
        se_mrna = np.std(mrna_values) / np.sqrt(len(mrna_values))
        ci_mrna = 1.96 * se_mrna  # 95% CI

        ap_data[ap_bin] = {'mean': mean_mrna, 'ci': ci_mrna}

# Plotting the results
fig, ax = plt.subplots(figsize=(5, 3))

for construct, mean_data in construct_means.items():
    color = construct_colors.get(construct, 'black')  # Default to black if construct not found

    # Plot max APBin with error bars, scaling the x-axis by 100
    ax.errorbar(
        mean_data['mean_max_apbin'] * 100, mean_data['mean_max_avg_mrna'],
        xerr=mean_data['se_max_apbin'] * 100, yerr=mean_data['se_max_avg_mrna'],
        fmt='o', capsize=5, color=color, label=f'{construct} Max'
    )

    # Plot left half-max APBin with error bars, scaling the x-axis by 100
    ax.errorbar(
        mean_data['mean_left_half_bin'] * 100, mean_data['mean_half_max_value'],
        xerr=mean_data['se_left_half_bin'] * 100, yerr=mean_data['se_half_max_value'],
        fmt='o', capsize=5, color=color, label=f'{construct} Left Half-Max'
    )

    # Plot right half-max APBin with error bars, scaling the x-axis by 100
    ax.errorbar(
        mean_data['mean_right_half_bin'] * 100, mean_data['mean_half_max_value'],
        xerr=mean_data['se_right_half_bin'] * 100, yerr=mean_data['se_half_max_value'],
        fmt='o', capsize=5, color=color, label=f'{construct} Right Half-Max'
    )

ax.set_xlabel('% egg length', fontsize=15)
ax.set_ylabel('integrated fluorescence (au)', fontsize=15)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
#ax.set_title('Average Max and Half-Max Fluorescence', fontsize=15)
ax.grid(True)
# ax.legend()
yticks = np.arange(1e6, 4e6, 0.5e6)
ax.set_yticks(yticks)

ax.set_xlim([40, 65])
#ax.set_xlim(35, 65)  # Adjusted for percentage scale
plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/Kr_maxhalfmax.png', dpi=300, bbox_inches='tight')
plt.show()


'''qq Plot section '''

'''import scipy.stats as stats

# Define the APBin range you want to include (example: 0.4 to 0.6)
apbin_min = 0.35
apbin_max = 0.55

# Collect filtered fluorescence values for each construct
construct_raw_values = defaultdict(list)

for construct in construct_list:
    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            construct_path = os.path.join(dropbox_folder, dir_name, 'BurstPropertiesSlope.mat')
            if os.path.exists(construct_path):
                try:
                    mat_data = scipy.io.loadmat(construct_path)
                    burst_properties = mat_data.get('BurstProperties', None)
                    if burst_properties is not None:
                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaNs and apply APBin range filter
                                valid_indices = (~np.isnan(total_mrna_values)) & \
                                                (ap_bin_values >= apbin_min) & \
                                                (ap_bin_values <= apbin_max)

                                filtered_mrna = total_mrna_values[valid_indices]
                                construct_raw_values[construct].extend(filtered_mrna.tolist())
                except Exception as e:
                    print(f"Error reading {construct_path} for QQ plot: {e}")

# QQ plot for each construct vs KrSE_0
reference_data = np.sort(construct_raw_values['KrSE_0'])

fig, axes = plt.subplots(1, len(construct_list)-1, figsize=(15, 5))

for idx, construct in enumerate([c for c in construct_list if c != 'KrSE_0']):
    test_data = np.sort(construct_raw_values[construct])

    min_len = min(len(reference_data), len(test_data))
    ref_quantiles = reference_data[:min_len]
    test_quantiles = test_data[:min_len]
    # Reference line
    ax.plot(ref_quantiles, ref_quantiles, 'k--', linewidth=1, label='y = x')

    # 95% confidence ribbon for QQ plot reference line
    # Use standard error of the mean or a fixed margin (adjust as needed)
    se = np.std(ref_quantiles - test_quantiles) / np.sqrt(min_len)
    ci_upper = ref_quantiles + 1.96 * se
    ci_lower = ref_quantiles - 1.96 * se

    ax.fill_between(ref_quantiles, ci_lower, ci_upper, color='gray', alpha=0.2, label='95% CI')

    ax = axes[idx]
    ax.scatter(ref_quantiles, test_quantiles, color=construct_colors[construct], alpha=0.6)
    ax.plot(ref_quantiles, ref_quantiles, 'k--', linewidth=1)  # reference line y = x
    ax.set_title(f'QQ Plot: {construct} vs KrSE_0\n(APBin {apbin_min:.2f}–{apbin_max:.2f})', fontsize=13)
    ax.set_xlabel('KrSE_0 Quantiles', fontsize=12)
    ax.set_ylabel(f'{construct} Quantiles', fontsize=12)
    ax.grid(True)

plt.tight_layout()
plt.show()'''
