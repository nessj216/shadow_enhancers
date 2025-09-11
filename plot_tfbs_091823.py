'''Plots TFBS along the length of an enhancer. the height of the TFBS bar corresponds to a TOBIAS-specific
strength metric

Inputs required for the following variables:
1) 'enhancer_length': length of enhancer
2) a path to your file
3) 'okaylist' is the names of your TFs from the file. this allows same TFBS to be plotted as the same color
4)color_map enter the names of TFs and the colors you want them'''





import matplotlib.pyplot as plt

# Read data from the text file
with open("path to your file", "r") as file:
    lines = file.readlines()
enhancer_length = 1585 #user input required
# Process each line and create the data list
data = []
okaylist = ['zelda','bcd', 'gt', 'Kr', 'hb', 'cad', 'tll', 'kni']
#okaylist = ['zelda','bcd', 'gt', 'Kr', 'hb']  # specific to evestripe2
#['zelda','bcd', 'gt', 'Kr', 'hb', 'cad', 'tll', 'kni']
for line in lines:
    line = line.strip()  # Remove leading/trailing whitespace
    elements = line.split("\t")  # Split line by comma

    site = {
        "start": int(elements[4]), #format of file
        "length": int(elements[2]),
        "score": float(elements[3]),
        "strand": elements[6],
        "type": elements[5]
    }
    data.append(site)



max_score = max(site["score"] for site in data)

color_map = {
    "zelda":"lightgreen",
    "gt": "lightblue",
    "hb": "orange",
    "bcd": "red",
    "Kr" : "midnightblue",
    "tll":"yellow",
    "kni":"mediumpurple",
    "cad": "slategrey",
    "Stat92E":"rosybrown"

}
# Set up the figure
fig, ax = plt.subplots()

# Configure the plot settings
bar_width = 0.8

# Iterate over the data
for site in data:
    start = site["start"]
    length = site["length"]
    score = site["score"]
    strand = site["strand"]
    type = site["type"]

    # Calculate the normalized score
    normalized_score = score / max_score

    # Calculate the x-coordinate range for the bar
    x = start
    x_range = [x, x + length]

    # Calculate the height of the bar based on the normalized score and strand
    if strand == '+' :
        y = normalized_score
    else:
        y = -normalized_score

    print(type, y)
    # Plot the bar

    if type in okaylist:
        print(type)
        ax.bar(x=x_range, height=y, width=length, align="edge", color=color_map[type])

# Set the x-axis limits
ax.set_xlim(0,enhancer_length)

# Set the y-axis limits
ax.set_ylim(-1, 1)

# Remove the y-axis
ax.yaxis.set_visible(False)

# Set the x-axis label
ax.set_xlabel("Position",  labelpad=-50,  fontsize=14)

# Set the title
#ax.set_title("Tobias Tool TFBS Predictions",  fontsize=18)

# Move the x-axis to y=0
ax.spines["bottom"].set_position("zero")
ax.spines["top"].set_visible(False)

# Remove the top spine
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

legend_handles = []
for type, color in color_map.items():
    if type in okaylist:
        legend_handles.append(plt.Line2D([0], [0], color=color, lw=4, label=type))

ax.legend(handles=legend_handles, loc="upper left")
plt.savefig("/Users/jillianness/Desktop/TFBS_tobias_figures/krprox_fixed", dpi=300, bbox_inches='tight')
plt.show()
# Show the plot

