import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Global style
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({
    'font.size': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18
})

def analyze_enhancers(file_path1, file_path2, genome_file_path):
    def calculate_percentages(df, te_col):
        total_lines = len(df)
        rna_types = ['LTR', 'SINE', 'LINE']
        dna_type = 'DNA'

        # Only for enhancer files (not full genome)
        grouped = df.groupby([0, 1, 2])[te_col].unique().reset_index() if te_col == 7 else None

        if grouped is not None:
            only_rna = grouped[grouped[te_col].apply(
                lambda x: all(any(rna in item for rna in rna_types) for item in x)
                          and not any(dna_type in item for item in x))]
            only_dna = grouped[grouped[te_col].apply(
                lambda x: all(dna_type in item for item in x)
                          and not any(any(rna in item for rna in rna_types) for item in x))]
            both = grouped[grouped[te_col].apply(
                lambda x: any(any(rna in item for rna in rna_types) for item in x)
                          and any(dna_type in item for item in x))]

            percent_rna = (len(only_rna) / len(grouped)) * 100
            percent_dna = (len(only_dna) / len(grouped)) * 100
            percent_both = (len(both) / len(grouped)) * 100
        else:
            percent_rna = percent_dna = percent_both = 0

        percent_ltr = (df[te_col].str.contains('LTR').sum() / total_lines) * 100
        percent_line = (df[te_col].str.contains('LINE').sum() / total_lines) * 100
        percent_sine = (df[te_col].str.contains('SINE').sum() / total_lines) * 100
        percent_dna_rc = (df[te_col].str.contains('DNA/RC').sum() / total_lines) * 100
        percent_dna_excl_rc = ((df[te_col].str.contains('DNA').sum() - df[te_col].str.contains('DNA/RC').sum()) / total_lines) * 100

        return [percent_rna, percent_dna, percent_both, percent_ltr, percent_line, percent_sine, percent_dna_rc, percent_dna_excl_rc]

    # Load files
    df1 = pd.read_csv(file_path1, sep='\t', header=None)
    df2 = pd.read_csv(file_path2, sep='\t', header=None)
    df3 = pd.read_csv(genome_file_path, sep='\t', header=None)

    percentages1 = calculate_percentages(df1, te_col=7)
    percentages2 = calculate_percentages(df2, te_col=7)
    percentages3 = calculate_percentages(df3, te_col=3)

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    for a in ax:
        a.grid(False)

    ########### Subplot 1: RNA/DNA/Both (original code) ###########
    file1_values = percentages1[:3]
    file2_values = percentages2[:3]
    new_categories = ["Shadows", "Singles"]
    x_positions = [0, 1]
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

    cumulative = 0
    for val, color in zip(file2_values[::-1], colors[::-1]):
        ax[0].bar(x_positions[1], val, width=0.5, bottom=cumulative, color=color, edgecolor='black')
        if val < 5:
            ax[0].text(x_positions[1] + 0.3, cumulative + val / 2, f'{val:.1f}%', ha='left', va='center', fontsize=18)
        else:
            ax[0].text(x_positions[1], cumulative + val / 2, f'{val:.1f}%', ha='center', va='center', fontsize=18)
        cumulative += val

    ax[0].set_ylabel('Percentage')
    ax[0].set_xticks(x_positions)
    ax[0].set_xticklabels(new_categories, fontsize=20)
    ax[0].tick_params(axis='both')
    ax[0].legend(frameon=True)

    ########### Subplot 2: TE Family Breakdown with Full Genome ###########
    te_labels = ['LTR', 'LINE', 'SINE', 'DNA/RC', 'DNA']
    te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#fdd49e', '#f16913']
    all_te_values = [percentages1[3:], percentages2[3:], percentages3[3:]]
    all_labels = ["Shadows", "Singles", "Genome"]
    x_positions = [0, 1, 2]

    for i, values in enumerate(all_te_values):
        bottom = 0
        for val, color in zip(values, te_colors):
            ax[1].bar(i, val, width=0.5, bottom=bottom, color=color, edgecolor='black')
            if val < 5:
                ax[1].text(i + 0.3, bottom + val / 2, f'{val:.1f}%', ha='left', va='center', fontsize=18)
            else:
                ax[1].text(i, bottom + val / 2, f'{val:.1f}%', ha='center', va='center', fontsize=18)
            bottom += val

    ax[1].set_ylabel('Percentage')
    ax[1].set_xticks(x_positions)
    ax[1].set_xticklabels(all_labels, fontsize=18)

    # Legend construction
    dummy_rna = Line2D([], [], color='none', label=r'$\mathbf{RNA\ Class}$')
    dummy_dna = Line2D([], [], color='none', label=r'$\mathbf{DNA\ Class}$')
    rna_handles = [
        Patch(facecolor=te_colors[0], edgecolor='black', label='LTR'),
        Patch(facecolor=te_colors[1], edgecolor='black', label='LINE'),
        Patch(facecolor=te_colors[2], edgecolor='black', label='SINE')
    ]
    dna_handles = [
        Patch(facecolor=te_colors[3], edgecolor='black', label='RC'),
        Patch(facecolor=te_colors[4], edgecolor='black', label='Other')
    ]
    handles = [dummy_rna] + rna_handles + [dummy_dna] + dna_handles
    ax[1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1.15, 1), frameon=True)

    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.tight_layout()

    plt.savefig("/Users/jillianness/Desktop/Figures_shadowbirth/v2_TE_class_and_family_breakdown.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
    plt.show()
analyze_enhancers(
    '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/FINAL_TE_shadow_ouput3.bed',
    '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/singels/filtered50bp_single_TEoutput.bed',
'/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/final_merged_cleaned_TE.bed')
##