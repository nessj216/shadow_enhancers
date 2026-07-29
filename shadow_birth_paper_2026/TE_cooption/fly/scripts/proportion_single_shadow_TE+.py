import matplotlib.pyplot as plt
import numpy as np
import seaborn
from pathlib import Path
# ── Global style and font settings ───────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial"
})

ROOT = Path(__file__).resolve().parents[1]
PLOT_FIGURES = ROOT / "plot_figures"
PLOT_FIGURES.mkdir(parents=True, exist_ok=True)

# ── Data ─────────────────────────────────────────────────────────────────────
labels  = ['shadows\n(n=1122)', 'singles\n(n=6673)']
te_pos  = np.array([37, 170])
totals  = np.array([1122, 6673])
te_neg  = totals - te_pos
width   = 0.5

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(3, 3))

# TE‑negative (grey, bottom segment)
ax.bar(labels, te_neg / totals, width,
       label='TE−', color='0.85', edgecolor='black',linewidth=1)

# TE‑positive (green, top segment)
ax.bar(labels, te_pos / totals, width,
       bottom=te_neg / totals,
       label='TE+', color='#5ca87c', edgecolor='black',linewidth=1)   # choose any green you prefer

# Annotate counts and percentages on TE‑positive segment
for i, (pos, total) in enumerate(zip(te_pos, totals)):
    pct = pos / total * 100
    # put label slightly above the bar top
    ax.text(i, 1.01, f'{pct:.1f}%', ha='center', va='bottom')
from scipy.stats import fisher_exact



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
ax.set_ylabel('percentage of enhancers', fontsize=12)
# give 10% headroom above 100%
ax.set_ylim(0, 1.24)
ax.set_yticks(np.linspace(0, 1, 6))
ax.set_yticklabels([f'{int(t*100)}' for t in np.linspace(0, 1, 6)])

#ax.set_yticklabels(['0', '100'], fontsize=16)
# Legend with black border
leg = ax.legend(frameon=True, facecolor='white', loc='best')
leg.get_frame().set_edgecolor('black')
leg.get_frame().set_linewidth(1)
# remove box, keep only axes
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# make axes black
ax.spines["left"].set_color("black")
ax.spines["bottom"].set_color("black")
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

# show y tick marks
ax.tick_params(axis="y", which="both",
               left=True, right=False,
               length=4, width=1,
               color="black", labelcolor="black")

# x axis styling too
ax.tick_params(axis="x", which="both",
               bottom=True, top=False,
               length=4, width=1,
               color="black", labelcolor="black")
ax.grid(False)
plt.tight_layout()
plt.savefig(PLOT_FIGURES / "proportion_single_shadow.png", dpi=600,
            bbox_inches='tight', pad_inches=0.1)

plt.show()


import matplotlib.pyplot as plt
import numpy as np
import seaborn
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.proportion import proportions_ztest

# ── Global style and font settings ───────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "font.family": "Arial"
})



# ── Stats ────────────────────────────────────────────────────────────────────
# 2×2 contingency table
#          overlaps   non-overlaps
# shadows      46       1122-46
# singles     915       6673-915
table = [
    [te_pos[0], te_neg[0]],
    [te_pos[1], te_neg[1]]
]

# Fisher's exact test
odds_ratio, fisher_p = fisher_exact(table)
print(f"Fisher’s exact test: OR = {odds_ratio:.4f}, p = {fisher_p:.3e}")

# Chi-square test
chi2, chi2_p, _, _ = chi2_contingency(table, correction=False)
print(f"Chi² test: χ² = {chi2:.4f}, p = {chi2_p:.3e}")

# Two-proportion z-test
z_stat, z_p = proportions_ztest(count=te_pos, nobs=totals, alternative='two-sided')
print(f"Two-proportion z-test: z = {z_stat:.4f}, p = {z_p:.3e}")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(3, 2.5))

# TE-negative (grey, bottom segment)
ax.bar(labels, te_neg / totals, width,
       label='TE−', color='0.85', edgecolor='black', linewidth=1)

# TE-positive (green, top segment)
ax.bar(labels, te_pos / totals, width,
       bottom=te_neg / totals,
       label='TE+', color='#5ca87c', edgecolor='black', linewidth=1)

# Annotate percentages above bars
for i, (pos, total) in enumerate(zip(te_pos, totals)):
    pct = pos / total * 100
    ax.text(i, 1.01, f'{pct:.1f}%', ha='center', va='bottom')



# Axes, legend, layout
ax.set_ylabel('percentage of enhancers', fontsize=12)
ax.set_ylim(0, 1.3)
ax.set_yticks(np.linspace(0, 1, 6))
ax.set_yticklabels([f'{int(t*100)}' for t in np.linspace(0, 1, 6)])

leg = ax.legend(frameon=True, facecolor='white', loc='best')
leg.get_frame().set_edgecolor('black')
leg.get_frame().set_linewidth(1)

ax.grid(False)
# remove box, keep only axes
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# make axes black
ax.spines["left"].set_color("black")
ax.spines["bottom"].set_color("black")
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

# show y tick marks
ax.tick_params(axis="y", which="both",
               left=True, right=False,
               length=4, width=1,
               color="black", labelcolor="black")

# x axis styling too
ax.tick_params(axis="x", which="both",
               bottom=True, top=False,
               length=4, width=1,
               color="black", labelcolor="black")
plt.tight_layout()
plt.savefig(
    PLOT_FIGURES / "proportion_single_shadow.png",
    dpi=600, bbox_inches='tight', pad_inches=0.1
)
plt.show()
