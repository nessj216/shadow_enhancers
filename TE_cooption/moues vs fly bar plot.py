import matplotlib.pyplot as plt

# Data
categories = ['Fly', 'Mouse']
subcategories = ['TE co-option', 'Duplication', 'Enhancer splitting', 'Other? De novo?']
values = [
    [10, 14, 1, 75],  # Fly percentages
    [70, 18, 3, 9]    # Mouse percentages
]

# Softer pastel colors (white for De novo)
colors = ['#a1d99b', '#ffe0b2', '#fcbba1', 'white']

# Create horizontal stacked bar chart with wider figure
fig, ax = plt.subplots(figsize=(12, 4))

# Add major grid lines on the x-axis
ax.grid(True, axis='x', linestyle='--', color='lightgrey', linewidth=0.7)

# Stack bars
left = [0, 0]
for idx, (subcat, color) in enumerate(zip(subcategories, colors)):
    bars = ax.barh(
        categories,
        [v[idx] for v in values],
        left=left,
        color=color,
        edgecolor='black',
        height=0.5,
        label=subcat
    )
    # If De novo, use thicker dotted border
    if subcat == 'De novo':
        for bar in bars:
            bar.set_linestyle(':')
            bar.set_linewidth(2.5)  # increased from 1.5 to 2.5 for a bolder dotted outline
    left = [left[i] + values[i][idx] for i in range(len(categories))]

# Add percentage labels inside each segment
for i, _ in enumerate(categories):
    start = 0
    for j, _ in enumerate(subcategories):
        pct = values[i][j]
        if pct > 0:
            ax.text(
                start + pct / 2,
                i,
                f'{pct}%',
                va='center',
                ha='center',
                color='black',
                fontsize=12,
                fontweight='bold'
            )
        start += pct

# Axes settings
ax.set_xlim(0, 100)
ax.set_xlabel('Percentage', fontsize=14)

# Customize y-axis labels to include sample sizes
ytick_labels = ['Fly\nn=1,122', 'Mouse\nn=22,755']
ax.set_yticks([0, 1])
ax.set_yticklabels(ytick_labels, fontsize=14)

ax.invert_yaxis()  # Fly on top

# Place legend outside on the right
legend = ax.legend(
    loc='center left',
    bbox_to_anchor=(1, 0.5),
    fontsize=12,
    frameon=True
)
legend.get_frame().set_edgecolor('black')

# Show all spines to create a perimeter
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_edgecolor('black')
    spine.set_linewidth(1)

# Remove x-axis tick lines
ax.tick_params(axis='x', which='both', length=0)
plt.xticks(fontsize=12)

# Adjust layout to fit legend
plt.tight_layout(rect=[0, 0, 0.75, 1])
plt.show()
