import matplotlib.pyplot as plt

# Read data from the text file
with open("/Users/jillianness/Desktop/Tobias1e2/D2/uniqueall.txt", "r") as file:
    lines = file.readlines()

# Process each line and create the data list
data = []
okaylist = ['Gt', 'Hb', 'Bcd', 'Kr']  # specific to evestripe2
#okaylist = ['Gt', 'Hb', 'Bcd', 'Kr', 'Tll', 'Kni', 'Cad',  'Stat92E']  # specific to evestripe2
for line in lines:
    line = line.strip()  # Remove leading/trailing whitespace
    elements = line.split("\t")  # Split line by tab

    site = {
        "start": int(elements[4]),
        "length": int(elements[2]),
        "score": float(elements[3]),
        "strand": elements[6],
        "type": elements[5]
    }
    data.append(site)

enhancer_length = 1483
max_score = max(site["score"] for site in data)

color_map = {
    "Gt": "skyblue",
    "Hb": "orange",
    "Bcd": "red",
    "Kr": "midnightblue",
    "Tll": "yellow",
    "Kni": "mediumpurple",
    "Cad": "slategrey",
    "Stat92E": "rosybrown",
    "zelda": "blue"


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
    if type == '+' :
        y = 3  # normalized_score
    else:
        y = -3  # -normalized_score

    print(type, y)
    # Plot the bar

    if type in okaylist:
        print(type)
        ax.bar(x=x_range, height=y, width=length, align="edge", color=color_map[type], alpha=0.7)

# Set the x-axis limits
ax.set_xlim(0, enhancer_length)

# Set the y-axis limits
ax.set_ylim(-1, 1)

# Remove the y-axis
ax.yaxis.set_visible(False)

# Set the x-axis label
ax.set_xlabel("Position", labelpad=-50, fontsize=12)

# Set the title
ax.set_title("Tobias Tool TFBS Predictions", fontsize=18)

# Move the x-axis to y=0
ax.spines["bottom"].set_position("zero")
ax.spines["top"].set_visible(False)

# Remove the top spine
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

# Show the plot


legend_labels = ['Gt', 'Hb', 'Bcd', 'Kr', 'Tll', 'Kni', 'Cad',  'Stat92E']
legend_colors = [color_map[label] for label in legend_labels]

# Plot the bars and add grey outlines as before

# Add a legend to the plot with correct labels and colors
ax.legend(legend_labels, loc='upper left', bbox_to_anchor=(1, 1),
          handles=[plt.Rectangle((0,0),1,1, color=color) for color in legend_colors],
          labels=legend_labels)

plt.show()