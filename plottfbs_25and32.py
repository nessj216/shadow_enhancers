import matplotlib.pyplot as plt

def read_data(file_path):
    # Read data from the text file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Process each line and create the data list
    data = []
    #okaylist = ['zelda','bcd', 'gt', 'Kr', 'hb', 'cad', 'tll', 'kni']
    okaylist = ['sna', 'twi', 'dl']
    for line in lines:
        line = line.strip()  # Remove leading/trailing whitespace
        elements = line.split("\t")  # Split line by comma

        site = {
            "start": int(elements[4]),
            "length": int(elements[2]),
            "score": float(elements[3]),
            "strand": elements[6],
            "type": elements[5]
        }
        data.append(site)

    return data, okaylist

def plot_data(ax, data, okaylist, enhancer_length, color_map, title, show_legend=True):
    max_score = max(site["score"] for site in data)
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
        if strand == '+':
            y = normalized_score
        else:
            y = -normalized_score

        # Plot the bar if the type is in the okaylist
        if type in okaylist:
            ax.bar(x=x_range, height=y, width=length, align="edge", color=color_map[type])

    # Configure the plot settings
    ax.set_xlim(0, enhancer_length)
    ax.set_ylim(-1, 1)
    ax.yaxis.set_visible(False)
    ax.set_xlabel("Position", labelpad=-50, fontsize=16)  # Increased font size
    ax.set_title(title, fontsize=18)  # Increased font size
    ax.tick_params(axis='x', labelsize=14)  # Increased tick label font size
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Ensure tick labels are shown
    ax.set_xticks(range(0, enhancer_length, 200))  # Adjust tick positions as needed
    ax.set_xticklabels(range(0, enhancer_length, 200), fontsize=14)

    # Add legend if needed
    if show_legend:
        legend_handles = []
        for type, color in color_map.items():
            if type in okaylist:
                legend_handles.append(plt.Line2D([0], [0], color=color, lw=4, label=type))
        ax.legend(handles=legend_handles, loc="upper left", fontsize=12)  # Adjust legend font size as needed

# Paths to the files for 32°C and 25°C conditions
file_32C = "/Users/jillianness/Downloads/TF_Binding_4-15-24_Sample_TOBIAS_chr/BINDetect_4-15_4-16_picard_Dl_Sna_Twist/snadistal_32_filesall.txt"
file_25C = "/Users/jillianness/Downloads/TF_Binding_4-15-24_Sample_TOBIAS_chr/BINDetect_4-15_4-16_picard_Dl_Sna_Twist/snadistal_25_filesall.txt"

# Read data for both conditions
data_32C, okaylist = read_data(file_32C)
data_25C, _ = read_data(file_25C)

# Define color map and enhancer length
color_map = {
    "zelda":"lightgreen",
    "gt": "lightskyblue",
    "hb": "orange",
    "bcd": "red",
    "Kr" : "midnightblue",
    "tll":"yellow",
    "kni":"mediumpurple",
    "cad": "slategrey",
    "Stat92E":"rosybrown",
    "twi": "darkslategray",
    "dl": 'thistle',
    "sna": "mediumaquamarine"
}

enhancer_length = 1137

# Set up the figure with 2 subplots without shared x-axis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot data for 32°C and 25°C conditions
plot_data(ax1, data_32C, okaylist, enhancer_length, color_map, "TFBS Predictions at 32°C", show_legend=False)
plot_data(ax2, data_25C, okaylist, enhancer_length, color_map, "TFBS Predictions at 25°C", show_legend=False)

# Adjust layout and save the plot
plt.tight_layout()
plt.savefig("/Users/jillianness/Desktop/TFBS_tobias_figures/snaildistal_32_25_comparison.png", dpi=300, bbox_inches='tight')
plt.show()
