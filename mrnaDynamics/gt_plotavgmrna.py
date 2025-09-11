import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from collections import defaultdict



# Define the construct list
#construct_list = ['gtSE_5', 'gtsquish_0','gtsquish_mini']

# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'
construct_list = ['gtSE_5', 'gtsquish_0','gtsquish_mini']
#construct_list = ['gtsquish_mini','gtsquish_0']
# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)


# Dictionary to hold the aggregated mRNA data for each construct
construct_data = {construct: defaultdict(lambda: [0, 0, 0]) for construct in construct_list}

# Colors for each construct plot
'''colors = {
    'KrSE_0': 'yellowgreen',
    'Krsquish_0': 'plum',
    'Krsquish_mini': 'slategrey'
}
'''

colors = {
    'gtSE_5': 'yellowgreen',
    'gtsquish_0': 'plum',
    'gtsquish_mini': 'slategrey',
    'gt_SEmini_endog': 'darkgreen'
}
#
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
                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bin_values = safe_extract_float(entry['APBin'])
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaN values from TotalmRNA
                                valid_indices = ~np.isnan(total_mrna_values)
                                ap_bin_values = ap_bin_values[valid_indices]
                                total_mrna_values = total_mrna_values[valid_indices]

                                # Aggregate TotalmRNA by APBin across embryos
                                for ap_bin, mrna in zip(ap_bin_values, total_mrna_values):
                                    construct_data[construct][ap_bin][0] += mrna  # Sum of TotalmRNA
                                    construct_data[construct][ap_bin][1] += 1  # Count for averaging
                                    construct_data[construct][ap_bin][2] += mrna ** 2  # Sum of squared TotalmRNA for CI

                except Exception as e:
                    print(f"Error loading {construct_path}: {e}")

# Plotting the average mRNA trace with 95% CI ribbon for each construct
plt.figure(figsize=(5, 3))

for construct, apbin_totals in construct_data.items():
    ap_bins = sorted(apbin_totals.keys())
    avg_mrna = [apbin_totals[ap_bin][0] / apbin_totals[ap_bin][1] for ap_bin in ap_bins]

    # Calculate 95% CI
    ci_95 = [
        1.96 * np.sqrt(
            (apbin_totals[ap_bin][2] / apbin_totals[ap_bin][1] - (avg_mrna[i]) ** 2) / apbin_totals[ap_bin][1])
        for i, ap_bin in enumerate(ap_bins)
    ]

    # Plot average mRNA traces with 95% CI ribbon
    plt.plot(ap_bins, avg_mrna, label=f'{construct}', color=colors[construct])
    plt.fill_between(ap_bins, np.array(avg_mrna) - np.array(ci_95), np.array(avg_mrna) + np.array(ci_95),
                     alpha=0.3, color=colors[construct])

plt.xlabel('% egg length',fontsize=15)
plt.ylabel('Integrated fluoresence intensity (AU)', fontsize=15)
#plt.title('Average mRNA Traces per Construct with 95% CI')

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(plt.xticks()[0], labels=[f'{int(x*100)}' for x in plt.xticks()[0]])

plt.xlim([.15, .45])

plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/minisdarkgreen colors_average_gtmrna_traces.png', dpi=300, bbox_inches='tight')

#plt.legend()
#plt.grid(True)
plt.show()


#raw nuclei values
import seaborn as sns
import pandas as pd

# Collect all raw TotalmRNA values by construct
violin_data = []

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
                            if 'TotalmRNA' in burst_properties.dtype.names:
                                total_mrna_values = safe_extract_float(entry['TotalmRNA'])

                                # Filter out NaNs
                                total_mrna_values = total_mrna_values[~np.isnan(total_mrna_values)]

                                # Append to violin_data
                                for val in total_mrna_values:
                                    violin_data.append({
                                        'Construct': construct,
                                        'TotalmRNA': val
                                    })

                except Exception as e:
                    print(f"Error loading for violin plot: {construct_path}: {e}")

# Create DataFrame for plotting
violin_df = pd.DataFrame(violin_data)

# Plot
plt.figure(figsize=(5, 4))  # Smaller figure for violin plot
sns.violinplot(x='Construct', y='TotalmRNA', data=violin_df, palette=colors, cut=0)

plt.xlabel('')
plt.ylabel('Integrated fluorescence intensity (AU)', fontsize=13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Optional: Limit y-axis or set scale
# plt.ylim([0, 3e6])
# plt.yscale("log")

plt.tight_layout()
plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/violin_gt_totalmRNA.png', dpi=300)
plt.show()









#Avg nuclei/ AP bin mrna distribution
import os
import numpy as np
import scipy.io
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define construct list and paths
construct_list = ['gtSE_5', 'gtsquish_0','gtsquish_mini']
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'
colors = {
    'gtSE_5': 'yellowgreen',
    'gtsquish_0': 'plum',
    'gtsquish_mini': 'slategrey'
}

# Helper to flatten safely
def safe_extract_float(nested_array):
    return np.ravel(nested_array).astype(float)

# Store AP bin means per construct
apbin_means_data = []

for construct in construct_list:
    apbin_totals = {}  # {ap_bin: [list of mRNA values]}

    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            file_path = os.path.join(dropbox_folder, dir_name, 'BurstPropertiesSlope.mat')
            if os.path.exists(file_path):
                try:
                    mat_data = scipy.io.loadmat(file_path)
                    burst_properties = mat_data.get('BurstProperties', None)

                    if burst_properties is not None:
                        for entry in burst_properties[0]:
                            if 'APBin' in burst_properties.dtype.names and 'TotalmRNA' in burst_properties.dtype.names:
                                ap_bins = safe_extract_float(entry['APBin'])
                                total_mrnas = safe_extract_float(entry['TotalmRNA'])

                                valid = ~np.isnan(ap_bins) & ~np.isnan(total_mrnas)
                                ap_bins = ap_bins[valid]
                                total_mrnas = total_mrnas[valid]

                                for ap_bin, mrna in zip(ap_bins, total_mrnas):
                                    apbin_totals.setdefault(ap_bin, []).append(mrna)

                except Exception as e:
                    print(f"Error processing {file_path}: {e}")

    # Now compute the mean TotalmRNA per AP bin and store it
    for ap_bin, mrna_vals in apbin_totals.items():
        if len(mrna_vals) > 0:
            mean_mrna = np.mean(mrna_vals)
            apbin_means_data.append({
                'Construct': construct,
                'AvgTotalmRNA_perAPbin': mean_mrna
            })
print(np.unique(list(apbin_totals.keys())))

# Convert to DataFrame
violin_df = pd.DataFrame(apbin_means_data)

# Plot violin of averaged per-AP-bin mRNA values
plt.figure(figsize=(5, 4))
sns.violinplot(
    x='Construct',
    y='AvgTotalmRNA_perAPbin',
    data=violin_df,
    palette=colors,
    cut=0
)

plt.xlabel('')
plt.ylabel('Avg nuclear mRNA per AP bin (AU)', fontsize=13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.tight_layout()
plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/violin_APbin_averages.png', dpi=300)
plt.show()

