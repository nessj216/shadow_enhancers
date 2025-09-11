import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import itertools

# Load the .mat file
mat = scipy.io.loadmat('/Users/jillianness/Downloads/myCellArrays110324.mat')

# Access the saved cell array
cell_array = mat['data_to_plot']

# Define the AP range for filtering
ap_min = 0.10
ap_max = 0.45

# Clean, filter by AP range, and prepare the data
filtered_data_list = []
for i in range(cell_array.shape[1]):
    cell_content = cell_array[0, i].squeeze()  # Squeeze to remove unnecessary dimensions
    # Ensure the content is an ndarray and contains numeric data
    if isinstance(cell_content, np.ndarray) and np.issubdtype(cell_content.dtype, np.number):
        # Remove NaN values and filter by AP range
        filtered_content = cell_content[~np.isnan(cell_content)]
        filtered_content = filtered_content[(filtered_content >= ap_min) & (filtered_content <= ap_max)]
        filtered_data_list.append(filtered_content)

# Set up the figure size to decrease space between violins
plt.figure(figsize=(5, 5))  # Adjust width for closer plots

# Plotting the violin plots with reduced space between them
sns.violinplot(data=filtered_data_list, palette=['blue', 'dodgerblue', 'orangered', 'goldenrod', 'green'], cut=0)

# Add x-axis labels with larger font size
labels = ['distal', '2x distal', 'proximal', '2xprox', 'Shadow']
plt.xticks(ticks=range(len(labels)), labels=labels, fontsize=15)  # Set larger font size for x-axis labels

# Adjust y-axis tick label size
plt.yticks(fontsize=15)

# Adjust subplots to minimize space between violins
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.9, wspace=0.02)

# Add median bars
for i, data in enumerate(filtered_data_list):
    median = np.median(data)
    plt.hlines(median, i - 0.25, i + 0.25, color='black', linewidth=2)

# Perform Mann-Whitney U tests for all combinations
combinations = list(itertools.combinations(range(len(filtered_data_list)), 2))

# Determine the maximum y-value for annotations
y_max = max([item for sublist in filtered_data_list for item in sublist]) + 2
y_text = y_max * 1.05  # Starting point for text annotations
bar_height = y_max * 0.0001  # Height of the significance bar

# Number of tests for Bonferroni correction
num_tests = len(combinations)

# Loop over all combinations and print p-values
for (idx1, idx2) in combinations:
    _, p_value = mannwhitneyu(filtered_data_list[idx1], filtered_data_list[idx2])
    corrected_p_value = p_value * num_tests
    print(f"Comparing {labels[idx1]} vs {labels[idx2]}: p-value (Bonferroni corrected) = {corrected_p_value:.10f}")

# Set y-axis label and title with larger font size
plt.ylabel('Coefficient of Variation', fontsize=14)
plt.title('Comparison of Groups with Mann-Whitney U Test (Bonferroni corrected)', fontsize=16)

# Adjust layout and save the plot
plt.tight_layout()
plt.ylim(0, 2.5)

plt.savefig('/Users/jillianness/Desktop/comittee_meeting_figures_2024/CV_GTall_violoinplot.png', dpi=300, bbox_inches='tight')
plt.show()
