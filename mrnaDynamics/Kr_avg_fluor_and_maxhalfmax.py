import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict

# Gaussian function for curve fitting
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Define constructs and colors
construct_list = ['KrSE_0', 'Krsquish_0', 'Krsquish_mini']
construct_colors = {
    'KrSE_0': 'yellowgreen',
    'Krsquish_0': 'plum',
    'Krsquish_mini': 'slategrey'
}

# Path to data
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Helper to extract values
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)

# Storage
construct_data = defaultdict(list)
construct_avg_traces = defaultdict(dict)

# Process each construct
for construct in construct_list:
    apbin_aggregate = defaultdict(list)

    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            construct_path = os.path.join(dropbox_folder, dir_name, 'BurstPropertiesSlope.mat')
            if os.path.exists(construct_path):
                try:
                    mat_data = scipy.io.loadmat(construct_path)
                    burst_properties = mat_data.get('BurstProperties', None)
                    if burst_properties is not None:
                        apbin_totals = defaultdict(lambda: [0, 0])
                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                valid = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid]
                                total_mrna_values = total_mrna_values[valid]

                                for ap, val in zip(ap_bin_values, total_mrna_values):
                                    apbin_aggregate[ap].append(val)
                                    apbin_totals[ap][0] += val
                                    apbin_totals[ap][1] += 1

                        # Compute avg mRNA per APBin per embryo
                        apbin_avg_mrna = {ap: total[0] / total[1] for ap, total in apbin_totals.items() if total[1] > 0}
                        ap_bins = sorted(apbin_avg_mrna.keys())
                        avg_mrna = [apbin_avg_mrna[ap] for ap in ap_bins]

                        # Fit Gaussian to avg trace
                        try:
                            popt, _ = curve_fit(gaussian, ap_bins, avg_mrna,
                                                p0=[max(avg_mrna), np.mean(ap_bins), np.std(ap_bins)])
                            fitted_max = popt[1]
                            half_max = popt[0] / 2

                            left_half = fitted_max - np.sqrt(2 * (np.log(popt[0] / half_max)) * popt[2]**2)
                            right_half = fitted_max + np.sqrt(2 * (np.log(popt[0] / half_max)) * popt[2]**2)

                            construct_data[construct].append({
                                'max_apbin': fitted_max,
                                'max_avg_mrna': popt[0],
                                'interpolated_left_half_max': left_half,
                                'interpolated_right_half_max': right_half,
                                'half_max_value': half_max
                            })

                        except RuntimeError:
                            print(f"Fit failed for: {dir_name}")
                except Exception as e:
                    print(f"Error loading: {construct_path} – {e}")

    # Compute avg trace + 95% CI across embryos
    for ap in sorted(apbin_aggregate.keys()):
        vals = np.array(apbin_aggregate[ap])
        mean = np.mean(vals)
        se = np.std(vals) / np.sqrt(len(vals))
        ci = 1.96 * se
        construct_avg_traces[construct][ap] = {'mean': mean, 'ci': ci}

# Compute mean ± SE for max and half-max values
construct_means = {}
for construct, data in construct_data.items():
    max_apbins = [d['max_apbin'] for d in data]
    max_mrnas = [d['max_avg_mrna'] for d in data]
    left_bins = [d['interpolated_left_half_max'] for d in data]
    right_bins = [d['interpolated_right_half_max'] for d in data]
    half_vals = [d['half_max_value'] for d in data]

    construct_means[construct] = {
        'mean_max_apbin': np.mean(max_apbins),
        'se_max_apbin': np.std(max_apbins) / np.sqrt(len(max_apbins)),
        'mean_max_avg_mrna': np.mean(max_mrnas),
        'se_max_avg_mrna': np.std(max_mrnas) / np.sqrt(len(max_mrnas)),
        'mean_left_half_bin': np.mean(left_bins),
        'se_left_half_bin': np.std(left_bins) / np.sqrt(len(left_bins)),
        'mean_right_half_bin': np.mean(right_bins),
        'se_right_half_bin': np.std(right_bins) / np.sqrt(len(right_bins)),
        'mean_half_max_value': np.mean(half_vals),
        'se_half_max_value': np.std(half_vals) / np.sqrt(len(half_vals))
    }
fig, ax = plt.subplots(figsize=(10, 6))

for construct in construct_list:
    color = construct_colors.get(construct, 'black')
    mean_data = construct_means[construct]

    # Plot avg trace with CI
    ap_bins_sorted = sorted(construct_avg_traces[construct].keys())
    means = [construct_avg_traces[construct][ap]['mean'] for ap in ap_bins_sorted]
    cis = [construct_avg_traces[construct][ap]['ci'] for ap in ap_bins_sorted]
    ap_bins_percent = [ap * 100 for ap in ap_bins_sorted]

    ax.plot(ap_bins_percent, means, label=f'{construct} avg trace', color=color, linewidth=1.8)
    ax.fill_between(ap_bins_percent,
                    np.array(means) - np.array(cis),
                    np.array(means) + np.array(cis),
                    color=color, alpha=0.3)

    # Plot max and half-max points
    ax.errorbar(mean_data['mean_max_apbin'] * 100, mean_data['mean_max_avg_mrna'],
                xerr=mean_data['se_max_apbin'] * 100, yerr=mean_data['se_max_avg_mrna'],
                fmt='o', capsize=5, color=color, label=f'{construct} Max')

    ax.errorbar(mean_data['mean_left_half_bin'] * 100, mean_data['mean_half_max_value'],
                xerr=mean_data['se_left_half_bin'] * 100, yerr=mean_data['se_half_max_value'],
                fmt='o', capsize=5, color=color, label=f'{construct} Left Half-Max')

    ax.errorbar(mean_data['mean_right_half_bin'] * 100, mean_data['mean_half_max_value'],
                xerr=mean_data['se_right_half_bin'] * 100, yerr=mean_data['se_half_max_value'],
                fmt='o', capsize=5, color=color, label=f'{construct} Right Half-Max')

ax.set_xlabel('% egg length', fontsize=14)
ax.set_ylabel('average fluorescence', fontsize=14)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.set_title('Average Fluorescence Trace with Max and Half-Max', fontsize=15)
ax.grid(True)
ax.set_xlim(35, 70)
ax.legend()
plt.tight_layout()
plt.show()
