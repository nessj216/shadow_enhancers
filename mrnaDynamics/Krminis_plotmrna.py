import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from collections import defaultdict



# Define the construct list


# Set the path to the Dropbox folder
dropbox_folder = r'/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'
#construct_list = ['gt23_0','5-m2']

#construct_list = ['KrSE_0','Krsquish_0','Krsquish_mini']
construct_list = ['Krsquish_mini', 'Kr_SEmini_endog']
#construct_list = ['Kr1_0','Kr_proxmini']
# Helper function to safely extract float values from nested arrays
def safe_extract_float(nested_array):
    flat_list = []
    if isinstance(nested_array, (np.ndarray, list)):
        flat_list.extend(np.ravel(nested_array).astype(float))
    return np.array(flat_list)


# Dictionary to hold the aggregated mRNA data for each construct
construct_data = {construct: defaultdict(lambda: [0, 0, 0]) for construct in construct_list}

# Colors for each construct plot
colors = {
    'KrSE_0': 'yellowgreen',
    'Krsquish_0': 'plum',
    'Krsquish_mini': 'slategrey',
    'Kr_SEmini_endog': 'darkgreen'
}


'''colors = {
    'Kr_dist_': 'royalblue',
     'Kr1_0': 'darkorange',
    'Kr_distalmini': 'deepskyblue',
    'Kr_proxmini' : "gold"
}'''
#
# Iterate through the construct list and load the BurstPropertiesSlope.mat file for each
for construct in construct_list:
    print(f"\nLoading data for construct: {construct}")
    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            construct_path = os.path.join(dropbox_folder, dir_name, 'BurstPropertiesSlope.mat')

            if os.path.exists(construct_path):
                print(f"  -> Found and loading file: {construct_path}")  # <-- Add this line

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
fig, ax = plt.subplots(figsize=(5, 3))

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

plt.xlim([.35, .7])
yticks = np.arange(0, 4e6, 0.5e6)
ax.set_yticks(yticks)
plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/KRmrna_traces.png', dpi=300, bbox_inches='tight')

#plt.legend()
#plt.grid(True)
plt.show()
