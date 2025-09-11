import scipy.io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the .mat file
mat_file_path = '/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/KrViolinPlot_CV_Data_AllBins.mat'
data = scipy.io.loadmat(mat_file_path)

# Access the main structure
combined_data = data['combinedDataStruct']

# Initialize lists to store AP bin and CV data
construct_names = []
cv_values = []
ap_bins = []

# Loop through each construct within the combinedDataStruct
for construct_data in combined_data[0, 0]:
    # Extract and convert construct name to a string
    construct_name = str(construct_data['ConstructName'][0])  # Convert construct name to string

    # Fully unpack and flatten AP bin IDs and CV values
    ap_bin_id = np.array(construct_data['APbinID'][0, 0]).flatten()
    cv_value = np.array(construct_data['CVValue'][0, 0]).flatten()

    # Filter to include only AP bins between 0.22 and 0.32
    mask = np.logical_and(ap_bin_id >= 0.3, ap_bin_id <= 0.37)  # Mask for specified range
    filtered_cv = cv_value[mask]
    filtered_ap = ap_bin_id[mask]

    # Extend lists with filtered values
    cv_values.extend(filtered_cv)
    ap_bins.extend(filtered_ap)
    construct_names.extend([construct_name] * len(filtered_cv))

# Create a DataFrame for easier plotting
data_dict = {
    'Construct Name': construct_names,
    'AP bin ID': ap_bins,
    'CV Value': cv_values
}
df = pd.DataFrame(data_dict)

# Filter the DataFrame to include only CV values up to 4
filtered_df = df[df['CV Value'] <= 4]

# Define a custom color palette using a list of colors
custom_palette = ['skyblue', 'salmon', 'lightgreen', 'orange', 'purple']

# Create violin plots with inner box plots for each construct with filtered CV values
plt.figure(figsize=(8, 4))
sns.violinplot(x='Construct Name', y='CV Value', data=filtered_df, inner='box', linewidth=1.2, palette=custom_palette, width=0.7)
plt.xlabel('Construct Name')
plt.ylabel('CV Value')
plt.title('Violin Plot of CV Values for Each Construct (CV <= 4) with Box Plot and Median')
plt.xticks(rotation=45)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()

