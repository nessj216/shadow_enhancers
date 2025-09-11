import matplotlib.pyplot as plt
import numpy as np

# Data
# Reversed lists
genes = ["Bcd", "Zld", "Stat", "Hb", "Cad", "Kni", "Gt"]
distal_values = [0, 0.340399232, 0.580323357, 4.860569534, 0, 1.629262956, 0]
proximal_values = [8.600575623, 10.64281329, 2.206528512, 10.82168052, 7.350592832, 8.629823406, 8.68190809]


# Reverse the order of genes and scores
genes.reverse()
distal_values.reverse()
proximal_values.reverse()

# Create an array of indices for the genes
y = np.arange(len(genes))

# Height of each bar
height = 0.45

# Create the figure and axes
fig, ax = plt.subplots()

# Create the bars for distal and proximal values
distal_bars = ax.barh(y - height/2, distal_values, height, label='Distal', color='royalblue')
proximal_bars = ax.barh(y + height/2, proximal_values, height, label='Proximal', color='darkorange',edgecolor='lightgrey')

# Set the labels and title
#ax.set_ylabel('Genes')
ax.set_xlabel('Normalized Additive Score',fontsize='12')
ax.set_title('$\it{Kr}$ Shadow Enhancers ',fontsize='14')
ax.set_yticks(y)
ax.set_yticklabels(genes,fontsize='12')
ax.legend(fontsize='11')
ax.tick_params(axis='x', labelsize=12)
# Show the plot
plt.tight_layout()
plt.show()
