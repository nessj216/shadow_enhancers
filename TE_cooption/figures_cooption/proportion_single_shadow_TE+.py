import matplotlib.pyplot as plt
import numpy as np

# ── Global style and font settings ───────────────────────────────────────────
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({
    'font.size': 20,        # default for labels, titles, text
    'xtick.labelsize': 20,
    'ytick.labelsize': 20
})

# ── Data ─────────────────────────────────────────────────────────────────────
labels  = ['Shadows\n(n=1122)', 'Singles\n(n=6673)']
te_pos  = np.array([46, 915])
totals  = np.array([1122, 6673])
te_neg  = totals - te_pos
width   = 0.5

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(4, 6))

# TE‑negative (grey, bottom segment)
ax.bar(labels, te_neg / totals, width,
       label='TE−', color='0.85', edgecolor='black',linewidth=1)

# TE‑positive (green, top segment)
ax.bar(labels, te_pos / totals, width,
       bottom=te_neg / totals,
       label='TE+', color='#74c476', edgecolor='black',linewidth=1)   # choose any green you prefer

# Annotate counts and percentages on TE‑positive segment
for i, (pos, total) in enumerate(zip(te_pos, totals)):
    pct = pos / total * 100
    ax.text(i, te_neg[i] / total + pos / total / 2,
            f'{pos} ({pct:.1f}%)',
            ha='center', va='center', fontsize=18)
from scipy.stats import fisher_exact

# your counts
te_pos   = np.array([46, 915])
totals   = np.array([1122, 6673])
te_neg   = totals - te_pos

# 2×2 contingency table
#        overlaps   non-overlaps
# shadows   46       1122-46
# singles  915       6673-915
table = [
    [te_pos[0], te_neg[0]],
    [te_pos[1], te_neg[1]]
]

odds_ratio, p_value = fisher_exact(table)
print(f"Fisher’s exact test: OR = {odds_ratio:.2f}, p = {p_value:.3e}")




from scipy.stats import chi2_contingency

chi2, chi2_p, _, _ = chi2_contingency(table, correction=False)
print(f"Chi² test: χ² = {chi2:.2f}, p = {chi2_p:.3e}")

# Axes, legend, layout
ax.set_ylabel('Proportion of enhancers',fontsize=20)
# give 10% headroom above 100%
ax.set_ylim(0, 1.10)
ax.set_yticks(np.linspace(0, 1, 6))
ax.set_yticklabels([f'{int(t*100)}' for t in np.linspace(0, 1, 6)])

#ax.set_yticklabels(['0', '100'], fontsize=16)
# Legend with black border
leg = ax.legend(frameon=True, facecolor='white', fontsize=20, loc='best')
leg.get_frame().set_edgecolor('black')
leg.get_frame().set_linewidth(1)

ax.grid(False)
plt.tight_layout()
plt.savefig("/Users/jillianness/Desktop/Figures_shadowbirth/proportion_single_shadow.png", dpi=600,
            bbox_inches='tight', pad_inches=0.1)

plt.show()
