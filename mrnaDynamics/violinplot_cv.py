
'''the purpose of the script is take a matlab dataframe of a construct noise data (some of the output from CV_final.m)
and plot in a violin plot and perform mann whitney test '''




import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import numpy
# Load the .mat file
mat = scipy.io.loadmat('/Users/jillianness/Downloads/myCellArrays031124.mat')

# The loaded .mat files use a dictionary-like structure in Python.
# To access the saved cell array:
cell_array = mat['data_to_plot']
'''
#check if there is data in array
for index, content in np.ndenumerate(cell_array):
    print(f"Cell {index}:")
    print(content)
    print('-' * 50)  # separator for clarity
    
print(type(cell_array[0][0].squeeze()))
print(cell_array[0][0].squeeze().dtype)
'''


for i in range(cell_array.shape[1]):
    cell_content = cell_array[0, i]
    if isinstance(cell_content, (np.ndarray, np.float64, np.int64)) and np.issubdtype(cell_content.dtype, np.number):
        cell_array[0, i] = cell_content[~np.isnan(cell_content)]


# Convert the MATLAB cell array to a list of numpy arrays and filter out NaN values
data_list =  [cell_array[0, i] for i in range(cell_array.shape[1])]

# Plot the violin plot
sns.violinplot(data=data_list)


# Define custom colors and labels
colors = ['dodgerblue', 'goldenrod','green']
labels = ['duplicate distal','duplicate proximal', 'shadow enhancer']

# Plot the violin plot with custom colors
sns.violinplot(data=data_list, palette=colors, fontsize=15)

# Add labels to the x-axis
plt.xticks(ticks=range(len(labels)), labels=labels,fontsize=15)

# Add a legend
patch_list = [plt.Line2D([0], [0], color=color, lw=4, label=label) for color, label in zip(colors, labels)]
#plt.legend(handles=patch_list)

'''
# Mann-Whitney U tests and annotate significance
comparisons = [(0, 1), (1, 2), (0, 2)]  # Comparisons: 1 vs 2, 2 vs 3, 1 vs 3
y_max = max([item for sublist in data_list for item in sublist]) + 2  # get an offset for placing text
y_text = 5.3

bar_height = 0.1
for (idx1, idx2) in comparisons:
    _, p_value = mannwhitneyu(data_list[idx1], data_list[idx2])
    print(idx1, idx2, p_value)
    if p_value < .0001:  # Choose a significance level
        plt.plot([idx1, idx2], [y_text, y_text], color='black')
        # Place the asterisk above the bar
        plt.text((idx1 + idx2) / 2, y_text + bar_height / 3, "***", ha='center', va='bottom')
        y_text += .25  # Adjust for the next text placement
'''
# Mann-Whitney U tests and annotate significance
comparisons = [(0, 1), (1, 2), (0, 2)]  # Comparisons: 1 vs 2, 2 vs 3, 1 vs 3
y_max = max([item for sublist in data_list for item in sublist]) + 2  # get an offset for placing text
y_text = 5.3

bar_height = 0.1

for (idx1, idx2) in comparisons:
    _, p_value = mannwhitneyu(data_list[idx1], data_list[idx2])
    print(idx1, idx2, p_value)
    if p_value < .0001:  # Choose a significance level
        plt.plot([idx1, idx2], [y_text, y_text], color='black')
        # Place the asterisk above the bar
        plt.text((idx1 + idx2) / 2, y_text + bar_height / 3, "***", ha='center', va='bottom')
        y_text += .25  # Adjust for the next text placement

plt.ylim(0, 3.5)
plt.ylabel('coefficient of variation',fontsize=15)
plt.title('')
plt.show()
