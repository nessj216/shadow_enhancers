import matplotlib.pyplot as plt
import numpy as np

# Data for the heatmap
buckets = ['2 shadows', '3-4 shadows', '>=5 shadows']
te_categories = ['1 TE', '2 TE', '3+ TE']
heatmap_data = [
    [10, 5, 0],   # '2 shadows' bucket: '1 TE', '2 TE', '3+ TE'
    [15, 10, 5],  # '3-4 shadows' bucket: '1 TE', '2 TE', '3+ TE'
    [20, 15, 10]  # '>=5 shadows' bucket: '1 TE', '2 TE', '3+ TE'
]

# Create the figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Create the heatmap with a white-to-blue gradient
cmap = plt.cm.Blues
im = ax.imshow(heatmap_data, cmap=cmap, aspect='auto')

# Add category labels
ax.set_xticks(np.arange(len(buckets)))
ax.set_yticks(np.arange(len(te_categories)))
ax.set_xticklabels(buckets)
ax.set_yticklabels(te_categories)
ax.set_xlabel("Total Shadows in Set")
ax.set_ylabel("TE-Derived Shadow Categories")
ax.set_title("Heatmap of TE-Derived Enhancer Counts")

# Add color bar for reference
cbar = ax.figure.colorbar(im, ax=ax)
cbar.set_label("Count of Sets")

# Show the plot
plt.tight_layout()
plt.show()
