import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt

# Define the construct list
construct_list = ['KrSE_0', 'Krsquish_0']  # Add more constructs as needed

# Set the path to the main Dropbox folder (update this path)
dropbox_folder = '/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/data_for_pythonplotting'

# Define number of AP bins
num_ap_bins = 41
ap_bin_edges = np.linspace(0, 1, num_ap_bins + 1)  # AP bins range from 0 to 1

# Define color mapping for each construct
construct_colors = {
    'KrSE_0': 'green',
    'Krsquish_0': 'blue',
    'gtsquish_mini': 'orange'
}

# Initialize storage for each construct
construct_data = {construct: [] for construct in construct_list}

# Iterate over each construct
for construct in construct_list:
    # Search for directories containing the construct name
    for dir_name in os.listdir(dropbox_folder):
        if construct in dir_name:
            prefix_path = os.path.join(dropbox_folder, dir_name, 'SpotCorrelationAdj.mat')
            if not os.path.exists(prefix_path):
                continue

            # Load the .mat file
            mat_data = scipy.io.loadmat(prefix_path)
            spot_diff = mat_data.get('SpotDiff', None)
            if spot_diff is None:
                continue

            # Initialize temporary storage for CVs in this embryo
            embryo_cvs = [[] for _ in range(num_ap_bins)]  # List of lists for each AP bin

            # Iterate over each row in SpotDiff
            for row in spot_diff[0]:  # Assuming spot_diff is a list of rows
                smooth_spot_one = row['SmoothSpotOne'].flatten()
                ap_bin = float(row['APBin'])  # Extract the AP bin as a float

                # Ensure AP bin is within the range of 0 to 1
                if 0 <= ap_bin <= 1:
                    # Find the corresponding AP bin index
                    ap_bin_idx = np.digitize(ap_bin, ap_bin_edges) - 1

                    # Calculate CV for SmoothSpotOne
                    valid_frames = ~np.isnan(smooth_spot_one)
                    if valid_frames.sum() > 1:  # Only consider spots with >1 frame
                        avg_fluo = np.nanmean(smooth_spot_one)
                        std_fluo = np.nanstd(smooth_spot_one)
                        cv_fluo = std_fluo / avg_fluo if avg_fluo else np.nan

                        # Append CV to the corresponding AP bin
                        embryo_cvs[ap_bin_idx].append(cv_fluo)

            # Append the embryo's CVs to the construct's data
            construct_data[construct].append(embryo_cvs)

# Prepare the plot
plt.figure(figsize=(12, 8))

# Calculate the average CV and 95% CI for each construct across embryos
for construct in construct_list:
    # Prepare to store the average CV and 95% CI across AP bins
    avg_cv = np.full(num_ap_bins, np.nan)
    ci_95 = np.full(num_ap_bins, np.nan)

    # Convert the list of lists to an array
    for ap_bin_idx in range(num_ap_bins):
        all_cvs = [cv for embryo in construct_data[construct] for cv in embryo[ap_bin_idx] if not np.isnan(cv)]

        if all_cvs:
            avg_cv[ap_bin_idx] = np.nanmean(all_cvs)
            ci_95[ap_bin_idx] = 1.96 * np.nanstd(all_cvs) / np.sqrt(len(all_cvs))  # 95% CI

    # Define the AP bin centers for plotting
    ap_bin_centers = (ap_bin_edges[:-1] + ap_bin_edges[1:]) / 2

    # Plot the continuous line for average CV
    plt.plot(ap_bin_centers, avg_cv, label=construct, color=construct_colors.get(construct, 'black'), linewidth=2)

    # Add a shaded region for the 95% CI ribbon
    plt.fill_between(ap_bin_centers, avg_cv - ci_95, avg_cv + ci_95, color=construct_colors.get(construct, 'black'), alpha=0.2)

# Customize the plot
plt.xlabel('% egg length', fontsize=15)
plt.ylabel('Coefficient of Variation (CV)', fontsize=15)
#plt.title('Average CV Across AP Bins for All Constructs with 95% CI', fontsize=16)
plt.xlim([.38, .68])
plt.xticks(plt.xticks()[0], labels=[f'{int(x*100)}' for x in plt.xticks()[0]])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0, 4])  # Adjust based on expected CV range
#plt.grid(True, linestyle='--', alpha=0.5)
#plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/1_KravCV_traces.png', dpi=300, bbox_inches='tight')

plt.show()
