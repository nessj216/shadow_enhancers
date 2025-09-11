import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('seaborn-whitegrid')
plt.rcParams.update({'font.size': 16})

def analyze_enhancers(file_path1, te_col):

    def calculate_percentages(df, te_col):
        total_lines = len(df)
        rna_types = ['LTR', 'SINE', 'LINE']
        dna_type = 'DNA'
        grouped = df.groupby([0, 1, 2])[te_col].unique().reset_index()

        only_rna = grouped[grouped[te_col].apply(
            lambda x: all(any(rna in item for rna in rna_types) for item in x)
                      and not any(dna_type in item for item in x))
        ]
        only_dna = grouped[grouped[te_col].apply(
            lambda x: all(dna_type in item for item in x)
                      and not any(any(rna in item for rna in rna_types) for item in x))
        ]
        both = grouped[grouped[te_col].apply(
            lambda x: any(any(rna in item for rna in rna_types) for item in x)
                      and any(dna_type in item for item in x))
        ]

        percent_rna = (len(only_rna) / len(grouped)) * 100
        percent_dna = (len(only_dna) / len(grouped)) * 100
        percent_both = (len(both) / len(grouped)) * 100

        percent_ltr = (df[te_col].str.contains('LTR').sum() / total_lines) * 100
        percent_line = (df[te_col].str.contains('LINE').sum() / total_lines) * 100
        percent_sine = (df[te_col].str.contains('SINE').sum() / total_lines) * 100
        percent_dna_rc = (df[te_col].str.contains('DNA/RC').sum() / total_lines) * 100
        percent_dna_excl_rc = ((df[te_col].str.contains('DNA').sum() - df[te_col].str.contains('DNA/RC').sum()) / total_lines) * 100

        return [percent_rna, percent_dna, percent_both, percent_ltr, percent_line, percent_sine, percent_dna_rc, percent_dna_excl_rc]

    # Load the file
    df1 = pd.read_csv(file_path1, sep='\t', header=None)
    percentages1 = calculate_percentages(df1, te_col)

    # Make the plot
    fig, ax = plt.subplots(1, 1, figsize=(6, 7))

    x_pos = 0
    colors = ["#9ecae1", "#fdae6b", "#7f7f7f"]
    labels = ['RNA Class+', 'DNA Class+', 'Both+']
    cumulative = 0

    for val, color, label in zip(percentages1[:3][::-1], colors[::-1], labels[::-1]):
        ax.bar(x_pos, val, width=0.5, bottom=cumulative, label=label, color=color, edgecolor='black')
        if val < 5:
            ax.text(x_pos + 0.3, cumulative + val / 2, f'{val:.1f}%', ha='left', va='center', fontsize=16)
        else:
            ax.text(x_pos, cumulative + val / 2, f'{val:.1f}%', ha='center', va='center', fontsize=16)
        cumulative += val

    ax.set_title('TE Class of Shadow TE+ Enhancers')
    ax.set_ylabel('Percentage', fontsize=18)
    ax.set_xticks([x_pos])
    ax.set_xticklabels(['Shadows'])
    ax.tick_params(axis='both', labelsize=18)
    ax.legend(frameon=True)

    plt.tight_layout()
    plt.savefig("/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/TEfamily/50percTEcomp/50peroverlap_additive_familyTE_enhancer.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
    plt.show()

# Example of how to call the updated function:
analyze_enhancers('/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/TEfamily/50percTEcomp/50peroverlap_additive_familyTE_enhancer.bed', te_col=9)


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Load only shadow and genome files
shadow_file = '/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/TEfamily/50percTEcomp/50peroverlap_additive_familyTE_enhancer.bed'
genome_file = '/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/filtered_mm10.fa.bed'

shadow_df = pd.read_csv(shadow_file, sep='\t', header=None)
genome_df = pd.read_csv(genome_file, sep='\t', header=None)

# TE family percent calculation
def get_te_family_percentages(df, te_col):
    total = len(df)
    ltr = df[te_col].str.contains('LTR').sum() / total * 100
    line = df[te_col].str.contains('LINE').sum() / total * 100
    sine = df[te_col].str.contains('SINE').sum() / total * 100
    dna_rc = df[te_col].str.contains('DNA/RC').sum() / total * 100
    dna = (df[te_col].str.contains('DNA').sum() - df[te_col].str.contains('DNA/RC').sum()) / total * 100
    return [ltr, line, sine, dna_rc, dna]

# Compute observed/expected for shadows only # specific column
shadow_obs = get_te_family_percentages(shadow_df, 9)
genome_exp = get_te_family_percentages(genome_df, 4)
shadow_oe = [obs / exp if exp > 0 else 0 for obs, exp in zip(shadow_obs, genome_exp)]

# Plot setup
labels = ['LTR', 'LINE', 'SINE', 'DNA/RC', 'DNA/Other']
te_colors = {
    'LTR': '#c6dbef',
    'LINE': '#6baed6',
    'SINE': '#1f78b4',
    'DNA/RC': '#fdd49e',
    'DNA/Other': '#f16913'
}
x = range(len(labels))
bar_width = 0.5

# Axis breaks
upper_ylim = (9, 9.5)
lower_ylim = (0, 2)

plt.style.use('seaborn-whitegrid')
plt.rcParams.update({'font.size': 18, 'font.family': 'DejaVu Sans'})

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, figsize=(8, 7),
    gridspec_kw={'height_ratios': [0.3, 2.7]}
)

# Top plot for outliers
bars1 = ax1.bar(x, shadow_oe, width=bar_width,
                color=[te_colors[label] for label in labels],
                edgecolor='black', hatch='//')
ax1.set_ylim(upper_ylim)
ax1.set_yticks([9, 9.5])
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(labeltop=False)

# Bottom plot for normal range
bars2 = ax2.bar(x, shadow_oe, width=bar_width,
                color=[te_colors[label] for label in labels],
                edgecolor='black', hatch='//')
ax2.set_ylim(lower_ylim)
ax2.set_yticks([0, 0.5, 1, 1.5, 2])
ax2.spines['top'].set_visible(False)
ax2.xaxis.tick_bottom()

# Annotate bars
for bar in bars2:
    height = bar.get_height()
    if height < 2:
        ax2.text(bar.get_x() + bar.get_width() / 2, height + 0.05, f'{height:.2f}', ha='center', va='bottom', fontsize=16)

for bar in bars1:
    height = bar.get_height()
    if height >= upper_ylim[0]:
        ax1.text(bar.get_x() + bar.get_width() / 2, height + 0.02, f'{height:.2f}', ha='center', va='bottom', fontsize=16)

# Break marks
d = .01
kwargs = dict(marker=[(-1, -1), (1, 1)], markersize=10,
              linestyle='none', color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Labels
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=18)
ax2.set_ylabel('Observed / Expected', fontsize=18)
ax1.set_title('Observed/Expected TE Family Enrichment (Shadows)', fontsize=18)

# Legend
shadow_patch = mpatches.Patch(facecolor='white', edgecolor='black', hatch='//', label='Shadows')
ax1.legend(handles=[shadow_patch], fontsize=16, loc='upper right')

print("Shadow percentages (Observed):", shadow_obs)
print("Genome percentages (Expected):", genome_exp)
print("Observed/Expected ratios:", shadow_oe)


# Save and show
plt.tight_layout()
plt.savefig("/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/TEfamily/50percTEcomp/50peroverlap_additive_familyTE_enhancer.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
plt.show()
