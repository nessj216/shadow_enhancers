import matplotlib.pyplot as plt
import pandas as pd

# Load your data (adjust the path to your file)
tfbs_data = pd.read_csv('/Users/jillianness/Desktop/zeba_snail_Kr/trl400.tsv', sep='\t')
figname="/Users/jillianness/Desktop/zeba_snail_Kr/trl400.png"
# Extract start, end positions, and strand information for TFBS
start_positions = tfbs_data['start']
end_positions = tfbs_data['stop']
strand_info = tfbs_data['strand']

# Define the actual enhancer length and the largest enhancer length
enhancer_length_actual = 400 # replace with the actual size of this enhancer
largest_enhancer_length = 2043

# Proportionally scale the enhancer
scaling_factor = 10 / largest_enhancer_length  # Scale the plot width to a manageable size
scaled_enhancer_length = enhancer_length_actual * scaling_factor

# Create the plot
plt.figure(figsize=(scaled_enhancer_length, 2))

# Plot enhancer as a rectangle, scaled proportionally
plt.gca().add_patch(plt.Rectangle((0, 0.5), scaled_enhancer_length, 0.5, color='orange', alpha=0.3, label="Enhancer"))

# Plot each TFBS as smaller rectangles, considering the strand orientation
for start, end, strand in zip(start_positions, end_positions, strand_info):
    scaled_start = start * scaling_factor
    scaled_end = end * scaling_factor

    if strand == '+':
        # Plot normally for '+' strand
        plt.gca().add_patch(plt.Rectangle((scaled_start, 0.5), scaled_end - scaled_start, 0.5, color='grey', alpha=0.8))
    else:
        # Reverse plotting for '-' strand
        plt.gca().add_patch(
            plt.Rectangle((scaled_enhancer_length - scaled_end, 0.5), scaled_end - scaled_start, 0.5, color='grey',
                          alpha=0.8))

# Customize plot
plt.title(f'TFBS Locations along the Enhancer ({enhancer_length_actual} bp), Strands Considered')
plt.xlim(0, scaled_enhancer_length)
plt.ylim(0, 1.5)
plt.gca().axis('off')  # Remove both x and y axes

# Save the plot as a high-resolution image (300 DPI)
plt.savefig(figname, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
