import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Load data
shadow_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/FINAL_TE_shadow_ouput3.bed'
all_enhancers_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/singels/filtered50bp_single_TEoutput.bed'
genome_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/final_merged_cleaned_TE.bed'

shadow_df = pd.read_csv(shadow_file, sep='\t', header=None)
all_enhancers_df = pd.read_csv(all_enhancers_file, sep='\t', header=None)
genome_df = pd.read_csv(genome_file, sep='\t', header=None)

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

# Observed/expected ratios
shadow_oe = [obs / exp if exp > 0 else 0 for obs, exp in zip(shadow_obs, genome_exp)]
all_oe = [obs / exp if exp > 0 else 0 for obs, exp in zip(all_obs, genome_exp)]

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
bar_width = 0.35

# Plot setup with broken y-axis
fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, figsize=(10, 9),
    gridspec_kw={'height_ratios': [0.3, 2.7]}
)

# Top axis (7.5–8)
ax1.bar([i - bar_width/2 for i in x], shadow_oe, width=bar_width,
        color=[te_colors[label] for label in labels], edgecolor='black', hatch='//')
ax1.bar([i + bar_width/2 for i in x], all_oe, width=bar_width,
        color=[te_colors[label] for label in labels], edgecolor='black')
ax1.set_ylim(7.5, 8)
ax1.set_yticks([7.5, 8])
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(labeltop=False)

# Bottom axis (0–2)
ax2.bar([i - bar_width/2 for i in x], shadow_oe, width=bar_width,
        color=[te_colors[label] for label in labels], edgecolor='black', hatch='//')
ax2.bar([i + bar_width/2 for i in x], all_oe, width=bar_width,
        color=[te_colors[label] for label in labels], edgecolor='black')
ax2.set_ylim(0, 2)
ax2.set_yticks([0, 0.5, 1, 1.5, 2])
ax2.spines['top'].set_visible(False)
ax2.xaxis.tick_bottom()

# Add labels above bars in bottom axis
for i, (shadow, allval) in enumerate(zip(shadow_oe, all_oe)):
    ax2.text(i - bar_width/2, shadow + 0.05, f'{shadow:.2f}', ha='center', va='bottom', fontsize=16)
    ax2.text(i + bar_width/2, allval + 0.05, f'{allval:.2f}', ha='center', va='bottom', fontsize=16)

# Break marks
d = .01
kwargs = dict(marker=[(-1, -1), (1, 1)], markersize=10,
              linestyle='none', color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Labels and legend
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=18)
ax2.set_ylabel('Observed / Expected', fontsize=18)
ax1.set_title('Observed/Expected TE Type Enrichment', fontsize=18)

shadow_patch = mpatches.Patch(facecolor='white', edgecolor='black', hatch='//', label='Shadows')
all_patch = mpatches.Patch(facecolor='white', edgecolor='black', label='All enhancers')
ax1.legend(handles=[shadow_patch, all_patch], fontsize=14, loc='upper right')

# Save high-res figure
plt.tight_layout()
plt.savefig("/Users/jillianness/Desktop/Figures_shadowbirth/TE_class_and_family_genome.png", dpi=600,
            bbox_inches='tight', pad_inches=0.1)
plt.show()

shadow_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/FINAL_TE_shadow_ouput3.bed'
all_enhancers_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/singels/filtered50bp_single_TEoutput.bed'
genome_file = '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/final_merged_cleaned_TE.bed'

shadow_df = pd.read_csv(shadow_file, sep='\t', header=None)
all_enhancers_df = pd.read_csv(all_enhancers_file, sep='\t', header=None)
genome_df = pd.read_csv(genome_file, sep='\t', header=None)

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
fig, ax = plt.subplots(figsize=(9, 7))

bars1 = ax.bar([i - bar_width/2 for i in x], shadow_oe, width=bar_width,
               color=[te_colors[label] for label in labels], edgecolor='black', hatch='//')
bars2 = ax.bar([i + bar_width/2 for i in x], all_oe, width=bar_width,
               color=[te_colors[label] for label in labels], edgecolor='black')

# Axis settings
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=18)
ax.set_ylabel(r'$\log_2\left(\frac{\mathrm{Observed}}{\mathrm{Expected}}\right)$', fontsize=18)
#ax.set_title('log₂ Observed/Expected TE Family Enrichment', fontsize=18)
ax.tick_params(axis='y', labelsize=16)
ax.tick_params(axis='x', labelsize=16)
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
                ha='center', va='bottom' if height >= 0 else 'top', fontsize=16)


# Add legend
shadow_patch = mpatches.Patch(facecolor='white', edgecolor='black', hatch='//', label='Shadows')
all_patch = mpatches.Patch(facecolor='white', edgecolor='black', label='All enhancers')
ax.legend(handles=[shadow_patch, all_patch], fontsize=18, loc='upper right')

plt.tight_layout()
plt.savefig("/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/TEfamily/TE_class_and_family_genome_log2_single.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.show()

