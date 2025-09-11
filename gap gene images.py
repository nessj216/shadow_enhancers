import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Initialize the figure
fig, ax = plt.subplots(figsize=(6, 4))

# Draw the embryo (represented as an ellipse)
embryo = patches.Ellipse((0.5, 0.5), width=1, height=0.7, fc='lightgrey')
ax.add_patch(embryo)

# Dummy gap gene expression patterns
# These are just simple rectangles for illustration purposes
# The coordinates are (x, y) for the lower-left corner, and then width and height

# Gene A (red)
gene_a = patches.Rectangle((0.1, 0.3), 0.2, 0.4, fc='red', alpha=0.5, label='Gene A')
ax.add_patch(gene_a)

# Gene B (green)
gene_b = patches.Rectangle((0.3, 0.25), 0.2, 0.5, fc='green', alpha=0.5, label='Gene B')
ax.add_patch(gene_b)

# Gene C (blue)
gene_c = patches.Rectangle((0.6, 0.3), 0.2, 0.4, fc='blue', alpha=0.5, label='Gene C')
ax.add_patch(gene_c)

# Adjust the plot limits and show
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.axis('off')
ax.legend(loc='upper right')

plt.title("Simplified Gap Gene Expression in Drosophila Embryo")
plt.show()
