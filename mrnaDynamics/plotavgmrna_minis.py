import scipy.io
import numpy as np
import matplotlib.pyplot as plt
'''
# Load MAT file
mat_path = '/Users/jillianness/Downloads/(No subject)-5/gt_minis_fluorescence_data.mat'
mat_data = scipy.io.loadmat(mat_path, struct_as_record=False, squeeze_me=True)

# Extract data
avgprod_all_ap = mat_data['AvgProdAllAP']
egg_length = mat_data['EggLength'].squeeze()

# Define construct names and colors
construct_list = ['gtSE', 'gtsquish','gtsquishier','gt_SEmini_endog']

colors = {
    'gtSE': 'yellowgreen',
    'gtsquish': 'plum',
    'gtsquishier': 'slategrey',
    'gt_SEmini_endog': 'green'
}

# Create plot
plt.figure(figsize=(10, 6))

for i, construct in enumerate(construct_list):
    data = avgprod_all_ap[i]
    avg = data.AvgProd.squeeze()
    ci = data.All95Conf.squeeze()

    # Filter out NaN values
    valid = ~np.isnan(avg) & ~np.isnan(ci) & ~np.isnan(egg_length)
    if np.sum(valid) == 0:
        print(f"No valid data for {construct}")
        continue

    x = egg_length[valid]
    y = avg[valid]
    yerr = ci[valid]

    plt.plot(x, y, label=construct, color=colors[construct])
    plt.fill_between(x, y - yerr, y + yerr, alpha=0.3, color=colors[construct])

plt.xlabel('% egg length', fontsize=15)
plt.ylabel('Integrated fluorescence (AU)', fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.tight_layout()
plt.savefig('avg_mrna_traces_all_xrange.png', dpi=300)
plt.show()


import scipy.io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the .mat file (update the path if needed)
mat_path = '/Users/jillianness/Downloads/gtminisViolinPlot_CV_Data.mat'
cv_data = scipy.io.loadmat(mat_path, struct_as_record=False, squeeze_me=True)
cv_struct = cv_data['ExportStruct']

# Define plot labels and color scheme
colors = {
    'gtSE': 'yellowgreen',
    'gtsquish': 'plum',
    'gtsquishier': 'slategrey',
    'gt_SEmini_endog': 'green'
}

# Extract CV data and label each construct
cv_values = []
cv_labels = []

for construct in colors:
    if hasattr(cv_struct, construct):
        cv_array = np.array(getattr(cv_struct, construct)).squeeze()
        valid_cv = cv_array[~np.isnan(cv_array)]
        cv_values.extend(valid_cv)
        cv_labels.extend([construct] * len(valid_cv))

# Create DataFrame for plotting
df = pd.DataFrame({'Construct': cv_labels, 'CV': cv_values})

# Plot the violin plot
plt.figure(figsize=(8, 6))
sns.violinplot(data=df, x='Construct', y='CV', palette=colors)
plt.ylabel('Temporal CV', fontsize=14)
plt.xlabel('Construct', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Temporal Noise (CV) Across AP Bins 13–18', fontsize=15)
plt.tight_layout()
plt.savefig('temporal_cv_violin_plot.png', dpi=300)
plt.show()'''


import scipy.io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load mRNA trace data
mat_path_mrna = '/Users/jillianness/Downloads/(No subject)-5/gt_minis_fluorescence_data.mat'
#mat_path_mrna='/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/Kr_v3minis_fluorescence_data.mat'
mat_data = scipy.io.loadmat(mat_path_mrna, struct_as_record=False, squeeze_me=True)
avgprod_all_ap = mat_data['AvgProdAllAP']
egg_length = mat_data['EggLength'].squeeze()

# Load CV data
#mat_path_cv = '/Volumes/rkc_wunderlichLab/Imaging Data/Jillian/18to23KrminisViolinPlot_CV_Data.mat'
mat_path_cv  = '/Users/jillianness/Downloads/gtminisViolinPlot_CV_Data.mat'
cv_data = scipy.io.loadmat(mat_path_cv, struct_as_record=False, squeeze_me=True)
cv_struct = cv_data['ExportStruct']

# Define constructs and colors
#construct_list = ['KrSE_0', 'Krsquish','Kr_SEmini_endog', 'Krsquish_mini']
construct_list = ['gtSE', 'gtsquish', 'gtsquishier', 'gt_SEmini_endog']
'''colors = {
    'KrSE_0': 'yellowgreen',
    'Krsquish': 'plum',
    'Krsquish_mini': 'slategrey',
    'Kr_SEmini_endog': 'green'
}'''
colors = {
    'gtSE': 'yellowgreen',
    'gtsquish': 'plum',
    'gtsquishier': 'slategrey',
    'gt_SEmini_endog': 'green'
}
# Extract CV data
cv_values = []
cv_labels = []

for construct in construct_list:
    if hasattr(cv_struct, construct):
        cv_array = np.array(getattr(cv_struct, construct)).squeeze()
        valid_cv = cv_array[~np.isnan(cv_array)]
        cv_values.extend(valid_cv)
        cv_labels.extend([construct] * len(valid_cv))

# Create CV DataFrame
df = pd.DataFrame({'Construct': cv_labels, 'CV': cv_values})

# Plot side-by-side subplots
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

# --- Left subplot: Avg mRNA traces ---
for i, construct in enumerate(construct_list):
    data = avgprod_all_ap[i]
    avg = data.AvgProd.squeeze()
    ci = data.All95Conf.squeeze()

    valid = ~np.isnan(avg) & ~np.isnan(ci) & ~np.isnan(egg_length)
    if np.sum(valid) == 0:
        continue

    x = egg_length[valid]
    y = avg[valid]
    yerr = ci[valid]

    axes[0].plot(x, y, label=construct, color=colors[construct])
    axes[0].fill_between(x, y - yerr, y + yerr, alpha=0.3, color=colors[construct])

axes[0].set_xlabel('% egg length', fontsize=14)
axes[0].set_ylabel('Integrated fluorescence (AU)', fontsize=14)
axes[0].legend()
axes[0].set_title('Average mRNA Traces', fontsize=15)

# --- Right subplot: CV Violin Plot ---
sns.violinplot(data=df, x='Construct', y='CV', palette=colors, ax=axes[1])
axes[1].set_ylabel('Temporal CV', fontsize=14)
axes[1].set_xlabel('')
axes[1].set_title('Temporal Noise (CV)', fontsize=15)
from scipy.stats import mannwhitneyu
from itertools import combinations
from statsmodels.stats.multitest import multipletests

# Perform pairwise Mann–Whitney U tests
pairs = list(combinations(construct_list, 2))
p_values = []

# Collect p-values
for c1, c2 in pairs:
    group1 = df[df['Construct'] == c1]['CV']
    group2 = df[df['Construct'] == c2]['CV']

    if len(group1) > 0 and len(group2) > 0:
        stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
    else:
        p = np.nan
    p_values.append(p)

# Bonferroni correction
_, p_adj, _, _ = multipletests(p_values, method='bonferroni')

# Plot annotations
y_max = df['CV'].max()
height_step = 0.05 * y_max
text_y = y_max + height_step

'''for i, ((c1, c2), p) in enumerate(zip(pairs, p_adj)):
    x1, x2 = construct_list.index(c1), construct_list.index(c2)
    axes[1].plot([x1, x1, x2, x2], [text_y, text_y + height_step, text_y + height_step, text_y], lw=1.5, c='black')
    axes[1].text((x1 + x2) / 2, text_y + height_step + 0.01, f'p={p:.2e}', ha='center', va='bottom', fontsize=11)
    text_y += height_step * 1.5  # Increment for next annotation to avoid overlap'''

plt.tight_layout()
plt.savefig('combined_mrna_cv_plot.png', dpi=300)
plt.show()


import scipy.io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load mRNA trace data

mat_data = scipy.io.loadmat(mat_path_mrna, struct_as_record=False, squeeze_me=True)
avgprod_all_ap = mat_data['AvgProdAllAP']
egg_length = mat_data['EggLength'].squeeze()

# Load CV data

cv_data = scipy.io.loadmat(mat_path_cv, struct_as_record=False, squeeze_me=True)
cv_struct = cv_data['ExportStruct']

# Select only the 3rd and 4th constructs
construct_list = ['gtsquishier', 'gt_SEmini_endog']
colors = {
    'gtsquishier': 'slategrey',
    'gt_SEmini_endog': 'green'
}

# Map construct names to their correct indices in avgprod_all_ap
construct_index_map = {
    'gtsquishier': 2,
    'gt_SEmini_endog': 3
}

# Extract CV data
cv_values = []
cv_labels = []

for construct in construct_list:
    if hasattr(cv_struct, construct):
        cv_array = np.array(getattr(cv_struct, construct)).squeeze()
        valid_cv = cv_array[~np.isnan(cv_array)]
        cv_values.extend(valid_cv)
        cv_labels.extend([construct] * len(valid_cv))

# Create CV DataFrame
df = pd.DataFrame({'Construct': cv_labels, 'CV': cv_values})

# Plot side-by-side subplots
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

# --- Left subplot: Avg mRNA traces ---
for construct in construct_list:
    i = construct_index_map[construct]
    data = avgprod_all_ap[i]
    avg = data.AvgProd.squeeze()
    ci = data.All95Conf.squeeze()

    valid = ~np.isnan(avg) & ~np.isnan(ci) & ~np.isnan(egg_length)
    if np.sum(valid) == 0:
        continue

    x = egg_length[valid]
    y = avg[valid]
    yerr = ci[valid]

    axes[0].plot(x, y, label=construct, color=colors[construct])
    axes[0].fill_between(x, y - yerr, y + yerr, alpha=0.3, color=colors[construct])

axes[0].set_xlabel('% egg length', fontsize=14)
axes[0].set_ylabel('Integrated fluorescence (AU)', fontsize=14)
axes[0].legend()
axes[0].set_title('Average mRNA Traces', fontsize=15)

# --- Right subplot: CV Violin Plot ---
sns.violinplot(data=df, x='Construct', y='CV', palette=colors, ax=axes[1])
axes[1].set_ylabel('Temporal CV', fontsize=14)
axes[1].set_xlabel('')
axes[1].set_title('Temporal Noise (CV)', fontsize=15)

from scipy.stats import mannwhitneyu

# Extract CV values for the two constructs
cv1 = df[df['Construct'] == construct_list[0]]['CV']
cv2 = df[df['Construct'] == construct_list[1]]['CV']

# Perform Mann–Whitney U test
stat, p_value = mannwhitneyu(cv1, cv2, alternative='two-sided')

# Annotate p-value on the violin plot
'''max_cv = max(df['CV'])
y_pos = max_cv + 0.1 * max_cv  # Position above violins
axes[1].text(0.5, y_pos, f'p = {p_value:.2e}', ha='center', fontsize=13)

# Optionally, add a line between violins
axes[1].plot([0, 1], [y_pos * 0.97, y_pos * 0.97], color='black', linewidth=1.5)'''


plt.tight_layout()
plt.savefig('combined_mrna_cv_plot.png', dpi=300)
plt.show()
