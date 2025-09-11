import matplotlib.pyplot as plt
import numpy as np

# Data
genes = ["Bcd", "Zld", "Stat", "Cad", "Hb", "Kr", "Tll", "Kni*"]
proximal_values = [0, 2.7375, 0.65, 0, 0, 8.739, 12.525, 0]
distal_values = [5.335, 0.4975, 0, 4.395, 2.067, 0.65, 0, 1.745]
genes.reverse()
distal_values.reverse()
proximal_values.reverse()
# Create an array of indices for the genes
y = np.arange(len(genes))

# Create the figure and axes
fig, ax = plt.subplots()
height = 0.8
# Create the horizontal bars for proximal and distal values (stacked)
proximal_bars = ax.barh(y, proximal_values, height, label='Proximal', color='darkorange', edgecolor='lightgrey')
distal_bars = ax.barh(y, distal_values, height, label='Distal', color='royalblue', edgecolor='lightgrey')

# Set the labels and title
ax.set_xlabel('normalized additive score', fontsize='12')

ax.set_title('$\it{gt}$ shadow enhancers', fontsize='14')
ax.set_yticks(y)
ax.set_yticklabels(genes, fontsize='12')
#ax.legend(fontsize='11')

# Add a dotted line between "Hb" bars
hb_index = genes.index("Hb")
ax.axhline(y=hb_index +.48, color='black', linestyle='dotted')
nu=hb_index +.48
ax.text(10.2, nu + .8, 'Activator', fontsize=11)
ax.text(10.15, nu-1, 'Repressor', fontsize=11)

# Show the plot
plt.tight_layout()
plt.show()
