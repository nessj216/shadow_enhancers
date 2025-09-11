
import pandas as pd
import matplotlib.pyplot as plt

# Use a clean style with default fonts and update font size globally
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({
    'font.size': 20,               # used by default for text, axis titles, etc.
    'xtick.labelsize': 18,
    'ytick.labelsize': 18
})

def analyze_enhancers(file_path1, file_path2):

    def calculate_percentages(df):
        total_lines = len(df)
        rna_types = ['LTR', 'SINE', 'LINE']
        dna_type = 'DNA'
        grouped = df.groupby([0, 1, 2])[7].unique().reset_index()

        only_rna = grouped[grouped[7].apply(
            lambda x: all(any(rna in item for rna in rna_types) for item in x)
                      and not any(dna_type in item for item in x))
        ]
        only_dna = grouped[grouped[7].apply(
            lambda x: all(dna_type in item for item in x)
                      and not any(any(rna in item for rna in rna_types) for item in x))
        ]
        both = grouped[grouped[7].apply(
            lambda x: any(any(rna in item for rna in rna_types) for item in x)
                      and any(dna_type in item for item in x))
        ]

        percent_rna = (len(only_rna) / len(grouped)) * 100
        percent_dna = (len(only_dna) / len(grouped)) * 100
        percent_both = (len(both) / len(grouped)) * 100

        percent_ltr = (df[7].str.contains('LTR').sum() / total_lines) * 100
        percent_line = (df[7].str.contains('LINE').sum() / total_lines) * 100
        percent_sine = (df[7].str.contains('SINE').sum() / total_lines) * 100
        percent_dna_rc = (df[7].str.contains('DNA/RC').sum() / total_lines) * 100
        percent_dna_excl_rc = ((df[7].str.contains('DNA').sum() - df[7].str.contains('DNA/RC').sum()) / total_lines) * 100

        return [percent_rna, percent_dna, percent_both, percent_ltr, percent_line, percent_sine, percent_dna_rc, percent_dna_excl_rc]

    df1 = pd.read_csv(file_path1, sep='\t', header=None)
    df2 = pd.read_csv(file_path2, sep='\t', header=None)

    percentages1 = calculate_percentages(df1)
    percentages2 = calculate_percentages(df2)

    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    for a in ax:
        a.grid(False, which='major')
        a.grid(False, which='minor')

    # Define new x-axis labels: Shadows for file1 and All enhancers for file2.
    new_categories = ["Shadows", "All enhancers"]
    x_positions = [0, 1]

    #### Plot 1: RNA, DNA, Both Comparison ####
    file1_values = percentages1[:3]
    file2_values = percentages2[:3]
    # Overall color scheme: RNA (blue: "#9ecae1"), DNA (orange: "#fdae6b"), Both (gray: "#7f7f7f")
    # Plot for file 1 ("Shadows")
    # Overall color scheme: RNA (blue: "#9ecae1"), DNA (orange: "#fdae6b"), Both (gray: "#7f7f7f")
    # Plot for file 1 ("Shadows")
    colors = ["#9ecae1", "#fdae6b", "#7f7f7f"]
    labels = ['RNA Class+', 'DNA Class+', 'Both+']
    cumulative = 0
    # Overall color scheme: RNA (blue: "#9ecae1"), DNA (orange: "#fdae6b"), Both (gray: "#7f7f7f")
    # Plot for file 1 ("Shadows")
    colors = ["#9ecae1", "#fdae6b", "#7f7f7f"]
    labels = ['RNA Class+', 'DNA Class+', 'Both+']
    cumulative = 0
    for val, color, label in zip(file1_values[::-1], colors[::-1], labels[::-1]):
        ax[0].bar(x_positions[0], val, width=0.5, bottom=cumulative, label=label, color=color, edgecolor='black')
        if val < 5:
            ax[0].text(x_positions[0] + 0.3, cumulative + val / 2, f'{val:.1f}%', ha='left', va='center', fontsize=18)
        else:
            ax[0].text(x_positions[0], cumulative + val / 2, f'{val:.1f}%', ha='center', va='center', fontsize=18)
        cumulative += val

    # Plot for file 2 ("All enhancers")
    cumulative = 0
    for val, color in zip(file2_values[::-1], colors[::-1]):
        ax[0].bar(x_positions[1], val, width=0.5, bottom=cumulative, color=color, edgecolor='black')
        if val < 5:
            ax[0].text(x_positions[1] + 0.3, cumulative + val / 2, f'{val:.1f}%', ha='left', va='center', fontsize=18)
        else:
            ax[0].text(x_positions[1], cumulative + val / 2, f'{val:.1f}%', ha='center', va='center', fontsize=18)
        cumulative += val

    #ax[0].set_title('TE Class of TE+ enhancers (RNA Class, DNA Class, Both)')
    ax[0].set_ylabel('Percentage')
    ax[0].set_xticks(x_positions)
    ax[0].set_xticklabels(new_categories, fontsize=18)
    ax[0].tick_params(axis='both')
    ax[0].legend(frameon=True)

    #### Plot 2: Detailed TE Type Breakdown for Both Files ####
    te_labels = ['LTR', 'LINE', 'SINE', 'DNA/RC', 'DNA']
    file1_te_values = percentages1[3:]
    file2_te_values = percentages2[3:]
    # Colors for detailed TE breakdown:
    # Blue shades for LTR, LINE, SINE and orange shades for DNA/RC, DNA.
    te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#fdd49e', '#f16913']


    # Plot for file 1 ("Shadows")
    # Plot for file 1 ("Shadows")
    bottom_file1 = 0
    for val, label, color in zip(file1_te_values, te_labels, te_colors):
        ax[1].bar(x_positions[0], val, width=0.5, bottom=bottom_file1, label=label, color=color, edgecolor='black')
        if val < 5:
            ax[1].text(x_positions[0] + 0.3, bottom_file1 + val / 2, f'{val:.1f}%', ha='left', va='center', fontsize=18)
        else:
            ax[1].text(x_positions[0], bottom_file1 + val / 2, f'{val:.1f}%', ha='center', va='center', fontsize=18)
        bottom_file1 += val

    # Plot for file 2 ("All enhancers")
    bottom_file2 = 0
    for val, color in zip(file2_te_values, te_colors):
        ax[1].bar(x_positions[1], val, width=0.5, bottom=bottom_file2, color=color, edgecolor='black')
        if val < 5:
            ax[1].text(x_positions[1] + 0.3, bottom_file2 + val / 2, f'{val:.1f}%', ha='left', va='center', fontsize=18)
        else:
            ax[1].text(x_positions[1], bottom_file2 + val / 2, f'{val:.1f}%', ha='center', va='center', fontsize=18)
        bottom_file2 += val

    #ax[1].set_title('TE Family of TE+ enhancers')
    ax[1].set_ylabel('Percentage')
    ax[1].set_xticks(x_positions)
    ax[1].set_xticklabels(new_categories, fontsize=18)
    ax[1].tick_params(axis='both')
    ax[1].legend(frameon=True)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    # Create header proxy items with no visible marker (using Line2D with color 'none')
    dummy_rna = Line2D([], [], color='none', label=r'$\mathbf{RNA\ Class}$')
    dummy_dna = Line2D([], [], color='none', label=r'$\mathbf{DNA\ Class}$')

    # Create patches for each of your groups using the colors chosen earlier.
    rna_handles = [
        Patch(facecolor='#c6dbef', edgecolor='black', label='LTR'),
        Patch(facecolor='#6baed6', edgecolor='black', label='LINE'),
        Patch(facecolor='#1f78b4', edgecolor='black', label='SINE')
    ]
    dna_handles = [
        Patch(facecolor='#fdd49e', edgecolor='black', label='RC'),
        Patch(facecolor='#f16913', edgecolor='black', label='Other')
    ]

    # Combine headers and handles into one list.
    handles = [dummy_rna] + rna_handles + [dummy_dna] + dna_handles

    ax[0].legend(loc='upper left', bbox_to_anchor=(1.13, 1), frameon=True)
    ax[1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1.15, 1), frameon=True)
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.tight_layout()
    # Optionally save the figure: plt.savefig('enhancer_analysis_dual.png', dpi=300, bbox_inches='tight')
    #plt.show()
# Save figure as high-resolution PNG
    plt.savefig("/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/TEfamily/TE_class_and_family_breakdown.png", dpi=600, bbox_inches='tight', pad_inches=0.1)

# Example usage (update the file paths as needed)

analyze_enhancers(
    '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/FINAL_TE_shadow_ouput3.bed',
    '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/singels/filtered50bp_single_TEoutput.bed'
)

## log graph
plt.style.use('seaborn-whitegrid')
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

shadow_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/FINAL_TE_shadow_ouput3.bed'
all_enhancers_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/singels/filtered50bp_single_TEoutput.bed'
genome_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/final_merged_cleaned_TE.bed'

shadow_df = pd.read_csv(shadow_file, sep='\t', header=None)
all_enhancers_df = pd.read_csv(all_enhancers_file, sep='\t', header=None)
genome_df = pd.read_csv(genome_file, sep='\t', header=None)
from scipy.stats import fisher_exact

# prepare counts for Fisher's tests
shadow_total = len(shadow_df)
single_total = len(all_enhancers_df)

# raw counts for each category
counts = {
    'LTR': (
        shadow_df[7].str.contains('LTR').sum(),
        all_enhancers_df[7].str.contains('LTR').sum()
    ),
    'LINE': (
        shadow_df[7].str.contains('LINE').sum(),
        all_enhancers_df[7].str.contains('LINE').sum()
    ),
    'SINE': (
        shadow_df[7].str.contains('SINE').sum(),
        all_enhancers_df[7].str.contains('SINE').sum()
    ),
    'DNA/RC': (
        shadow_df[7].str.contains('DNA/RC').sum(),
        all_enhancers_df[7].str.contains('DNA/RC').sum()
    ),
    'DNA/Other': (
        # total "DNA" occurrences minus the DNA/RC ones
        shadow_df[7].str.contains('DNA').sum() - shadow_df[7].str.contains('DNA/RC').sum(),
        all_enhancers_df[7].str.contains('DNA').sum() - all_enhancers_df[7].str.contains('DNA/RC').sum()
    )
}

print("Fisher’s exact test (Shadow vs Single) for each TE category:")
for label, (shadow_count, single_count) in counts.items():
    # build 2×2 table:
    #               TE present      TE absent
    # Shadow          a                 b
    # Single          c                 d
    table = [
        [shadow_count, shadow_total - shadow_count],
        [single_count, single_total - single_count]
    ]
    oddsratio, pvalue = fisher_exact(table)
    print(f"{label:10s}  Shadow: {shadow_count:4d}/{shadow_total:<4d}   "
          f"Single: {single_count:4d}/{single_total:<4d}   "
          f"OR={oddsratio:.2f}   p={pvalue:.2e}")

# Function to calculate TE family percentages
def get_te_family_percentages(df, te_col):
    total = len(df)
    ltr = df[te_col].str.contains('LTR').sum() / total * 100
    line = df[te_col].str.contains('LINE').sum() / total * 100
    sine = df[te_col].str.contains('SINE').sum() / total * 100
    dna_rc = df[te_col].str.contains('DNA/RC').sum() / total * 100
    dna = (df[te_col].str.contains('DNA').sum() - df[te_col].str.contains('DNA/RC').sum()) / total * 100
    return [ltr, line, sine, dna_rc, dna]

# Compute observed and expected
shadow_obs = get_te_family_percentages(shadow_df, 7)
all_obs = get_te_family_percentages(all_enhancers_df, 7)
genome_exp = get_te_family_percentages(genome_df, 3)

# Observed/expected ratios (log2-transformed)
shadow_oe = [np.log2(obs / exp) if exp > 0 and obs > 0 else 0 for obs, exp in zip(shadow_obs, genome_exp)]
all_oe = [np.log2(obs / exp) if exp > 0 and obs > 0 else 0 for obs, exp in zip(all_obs, genome_exp)]





# Plot config
labels = ['LTR', 'LINE', 'SINE', 'DNA/RC', 'DNA/Other']
te_colors = {
    'LTR': '#c6dbef',
    'LINE': '#6baed6',
    'SINE': '#1f78b4',
    'DNA/RC': '#fdd49e',
    'DNA/Other': '#f16913'
}
x = range(len(labels))
bar_width = 0.46

# Set up single plot
fig, ax = plt.subplots(figsize=(8, 6))

bars1 = ax.bar([i - bar_width/2 for i in x], shadow_oe, width=bar_width,
               color=[te_colors[label] for label in labels], edgecolor='black', hatch='//')
bars2 = ax.bar([i + bar_width/2 for i in x], all_oe, width=bar_width,
               color=[te_colors[label] for label in labels], edgecolor='black')

# Axis settings
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=20)
ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Expected}}\right)$', fontsize=18)
#ax.set_title('log₂ Observed/Expected TE Family Enrichment', fontsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(axis='x', labelsize=18)
ax.set_ylim(-2, 4)
ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
ax.grid(axis='x', which='both', linestyle='', linewidth=0)

# Add value labels on top of bars
# Add value labels clearly above each bar, even for negative bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        offset = 0.1 if height >= 0 else -0.3  # push label up for positive, down for negative
        ax.text(bar.get_x() + bar.get_width() / 2, height + offset, f'{height:.2f}',
                ha='center', va='bottom' if height >= 0 else 'top', fontsize=18)


# Add legend
shadow_patch = mpatches.Patch(facecolor='white', edgecolor='black', hatch='//', label='Shadow enhancers')
all_patch = mpatches.Patch(facecolor='white', edgecolor='black', label='Single enhancers')
ax.legend(handles=[shadow_patch, all_patch], fontsize=20, loc='upper right')

plt.tight_layout()
plt.savefig("/Users/jillianness/Desktop/Figures_shadowbirth/v2TE_class_and_family_genome_log2_single.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.show()

