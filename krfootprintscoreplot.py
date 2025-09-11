import matplotlib.pyplot as plt
import numpy as np


genes = ["Bcd", "Zld", "Stat", "Hb", "Cad", "Kni", "Gt",'Kr']
distal_values = [0, 0.340399232, 0.580323357, 4.860569534, 0, 1.629262956, 0]
proximal_values = [8.600575623, 10.64281329, 2.206528512, 10.82168052, 7.350592832, 8.629823406, 8.68190809]

genes.reverse()
distal_values.reverse()
proximal_values.reverse()
# Create an array of indices for the genes
y = np.arange(len(genes))

# Create the figure and axes
fig, ax = plt.subplots()
height = 0.7
# Create the horizontal bars for proximal and distal values (stacked)
proximal_bars = ax.barh(y, proximal_values, height, label='Proximal', color='darkorange', edgecolor='lightgrey')
distal_bars = ax.barh(y, distal_values, height, label='Distal', color='royalblue', edgecolor='lightgrey')

# Set the labels and title
ax.set_xlabel('normalized additive score', fontsize='12')

ax.set_title('$\it{Kr}$ shadow enhancers', fontsize='14')
ax.set_yticks(y)
ax.set_yticklabels(genes, fontsize='12')
#ax.legend(fontsize='11')



# Add a dotted line between "Hb" bars
hb_index = genes.index("Hb")
ax.axhline(y=hb_index+.6, color='black', linestyle='dotted')

ax.axhline(y=hb_index-.6, color='black', linestyle='dotted')

ax.text(1.8, hb_index-.1 , 'Activator', fontsize=11)
#ax.text(9.2, hb_index + 2.8, 'Activator', fontsize=11)

ax.text(9.5, hb_index + 3, 'Activator', fontsize=11)
ax.text(9.5, hb_index-3.2, 'Repressor', fontsize=11)

ax.text(7.15, hb_index-.1, 'Repressor', fontsize=11)
#ax.text(9.15, hb_index-3, 'Repressor', fontsize=11)
# Show the plot
plt.tight_layout()
plt.show()
