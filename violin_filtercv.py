import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import itertools

# Load the .mat file
mat = scipy.io.loadmat('/Users/jillianness/Downloads/myCellArraysgt0227.mat')

# Access the saved cell array
cell_array = mat['data_to_plot']

# Clean and prepare the data
for i in range(cell_array.shape[1]):
    cell_content = cell_array[0, i].squeeze()  # Squeeze to remove unnecessary dimensions
    # Ensure the content is an ndarray and contains numeric data
    if isinstance(cell_content, np.ndarray) and np.issubdtype(cell_content.dtype, np.number):
        # Remove NaN values and values over 4
        cell_array[0, i] = cell_content[~np.isnan(cell_content) & (cell_content <= 4)]

# Convert to a list of numpy arrays
data_list = [cell_array[0, i] for i in range(cell_array.shape[1]) if np.any(cell_array[0, i]) <= 4]

# Plotting
sns.violinplot(data=data_list, palette=['yellowgreen', 'lightcoral', 'darkgreen', 'darkred'])
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
y_text = y_max * 1.05  # Starting point f
