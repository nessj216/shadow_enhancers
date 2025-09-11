import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import itertools

# Load the .mat file
mat = scipy.io.loadmat('/Users/jillianness/Downloads/myCellArrays0213.mat')

# Access the saved cell array
cell_array = mat['data_to_plot']

# Clean and prepare the data
for i in range(cell_array.shape[1]):
    cell_content = cell_array[0, i].squeeze()  # Squeeze to remove unnecessary dimensions
    # Ensure the content is an ndarray and contains numeric data
    if isinstance(cell_content, np.ndarray) and np.issubdtype(cell_content.dtype, np.number):
        # Remove NaN values
        cell_array[0, i] = cell_content[~np.isnan(cell_content)]

# Convert to a list of numpy arrays
data_list = [cell_array[0, i] for i in range(cell_array.shape[1])]

# Plotting
sns.violinplot(data=data_list, palette=['yellowgreen', 'plum', 'darkgreen', 'purple'])
labels = ['gt_SE_25C', 'gtsquish_25', 'gt_SE_32C', 'gtsquish_32']
plt.xticks(ticks=range(len(labels)), labels=labels)

# Add median bars
for i, data in enumerate(data_list):
    median = np.median(data)
    plt.hlines(median, i - 0.25, i + 0.25, color='black', linewidth=2)


# Perform Mann-Whitney U tests for all combinations
combinations = list(itertools.combinations(range(len(data_list)), 2))

# Determine the maximum y-value for annotations
y_max = max([item for sublist in data_list for item in sublist]) + 2
y_text = y_max * 1.05  # Starting point for text annotations

bar_height = y_max * 0.001  # Height of the significance bar

# Loop over all combinations
for (idx1, idx2) in combinations:
    _, p_value = mannwhitneyu(data_list[idx1], data_list[idx2])
    print(f"Comparing {labels[idx1]} vs {labels[idx2]}: p-value = {p_value:.4f}")

  # Annotate if p-value is less than 0.05 (common significance level)
    '''if p_value < 0.05:
        # Draw a line between the pairs and a star for significance
        plt.plot([idx1, idx2], [y_text, y_text], color='black')
        sig_symbol = '*' if p_value < 0.05 else 'n.s.'  # Could add more symbols for different significance levels
        plt.text((idx1 + idx2) / 2, y_text + bar_height / 2, sig_symbol, ha='center', va='bottom')
        y_text += bar_height * 3  # Increase y_text for next annotation to avoid overlap'''

plt.ylabel('Coefficient of Variation')
plt.title('Comparison of Groups with Mann-Whitney U Test')
plt.tight_layout()  # Adjust layout to not cut off labels
plt.show()
