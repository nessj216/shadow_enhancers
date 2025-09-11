import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Use consistent style and font settings for all figures
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({'font.size': 18, 'font.family': 'DejaVu Sans'})

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

    new_categories = ["Shadows", "All enhancers"]
    x_positions = [0, 1]

    # Plot 1: RNA/DNA/Both breakdown
    file1_values = percentages1[:3]
    file2_values = percentages2[:3]
    colors = ["#9ecae1", "#fdae6b", "#7f7f7f"]
    labels = ['RNA Class+', 'DNA Class+', 'Both+']

    for idx, (file_values, xpos) in enumerate(zip([file1_values, file2_values], x_positions)):
        cumulative = 0
        for val, color in zip(file_values[::-1], colors[::-1]):
            ax[0].bar(xpos, val, width=0.5, bottom=cumulative, color=color, edgecolor='black')
            if val < 5:
                ax[0].text(xpos + 0.3, cumulative + val / 2, f'{val:.1f}%', ha='left', va='center')
            else:
                ax[0].text(xpos, cumulative + val / 2, f'{val:.1f}%', ha='center', va='center')
            cumulative += val

    ax[0].set_title('TE Class of TE+ Enhancers')
    ax[0].set_ylabel('Percentage')
    ax[0].set_xticks(x_positions)
    ax[0].set_xticklabels(new_categories)
    ax[0].legend(labels=labels, frameon=True)

    # Plot 2: TE family breakdown
    te_labels = ['LTR', 'LINE', 'SINE', 'DNA/RC', 'DNA']
    file1_te_values = percentages1[3:]
    file2_te_values = percentages2[3:]
    te_colors = ['#c6dbef', '#6baed6', '#1f78b4', '#fdd49e', '#f16913']

    for idx, (file_values, xpos) in enumerate(zip([file1_te_values, file2_te_values], x_positions)):
        cumulative = 0
        for val, color in zip(file_values, te_colors):
            ax[1].bar(xpos, val, width=0.5, bottom=cumulative, color=color, edgecolor='black')
            if val < 5:
                ax[1].text(xpos + 0.3, cumulative + val / 2, f'{val:.1f}%', ha='left', va='center')
            else:
                ax[1].text(xpos, cumulative + val / 2, f'{val:.1f}%', ha='center', va='center')
            cumulative += val

    ax[1].set_title('TE Family of TE+ Enhancers')
    ax[1].set_ylabel('Percentage')
    ax[1].set_xticks(x_positions)
    ax[1].set_xticklabels(new_categories)

    # Legend construction
    dummy_rna = Line2D([], [], color='none', label=r'$\mathbf{RNA\ Class}$')
    dummy_dna = Line2D([], [], color='none', label=r'$\mathbf{DNA\ Class}$')
    rna_handles = [
        Patch(facecolor='#c6dbef', edgecolor='black', label='LTR'),
        Patch(facecolor='#6baed6', edgecolor='black', label='LINE'),
        Patch(facecolor='#1f78b4', edgecolor='black', label='SINE')
    ]
    dna_handles = [
        Patch(facecolor='#fdd49e', edgecolor='black', label='RC'),
        Patch(facecolor='#f16913', edgecolor='black', label='Other')
    ]
    handles = [dummy_rna] + rna_handles + [dummy_dna] + dna_handles

    ax[0].legend(loc='upper left', bbox_to_anchor=(1.13, 1), frameon=True)
    ax[1].legend(handles=handles, loc='upper left', bbox_to_anchor=(1.15, 1), frameon=True)

    plt.tight_layout()
    #plt.savefig("/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/TEfamily/TE_class_and_family_breakdown.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
    plt.show()
analyze_enhancers(
    '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/FINAL_TE_shadow_ouput3.bed',
    '/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/singels/filtered50bp_single_TEoutput.bed'
)