import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "legend.fontsize": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 600,
    "savefig.dpi": 600,
})

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLOT_FIGURES = ROOT / "plot_figures"
PLOT_FIGURES.mkdir(parents=True, exist_ok=True)

# ====== Inputs ======
te_minus_sets = 376   # terminus– (TE–) sets
te_plus_sets  = 36    # terminus+ (TE+) sets
total_sets = te_minus_sets + te_plus_sets

te_plus_breakdown = {
    "1 TE+ enhancer/set": 30,
    "2 TE+ enhancers/set": 5,
    "3 TE+ enhancers/set": 1,
}

# ====== Appearance ======
C_PLUS  = "#5ca87c"   # green (replaces old blue)
C_MINUS = "#E5E7E9"   # light grey

# inner ring: still green, but distinguish slices
C_1 = "#ddeea2"
C_2 = "#bcbc11"
C_3 = "#678905"

TRANSPARENT = (1, 1, 1, 0)
DPI = 600

def nested_donut_terminus_outer_enhancers_inner(ax, startangle=90):
    # matches your original ring geometry
    outer_radius = 1.0
    outer_width  = 0.32
    inner_radius = 0.68
    inner_width  = 0.28

    # ---- Outer ring: terminus – / + ----
    ax.pie(
        [te_plus_sets, te_minus_sets],
        colors=[C_PLUS, C_MINUS],
        startangle=startangle,
        radius=outer_radius,
        wedgeprops=dict(width=outer_width, edgecolor="black", linewidth=.3)
    )

    # ---- Inner ring: enhancers/set within TE+ ----
    v1 = te_plus_breakdown.get("1 TE+ enhancer/set", 0)
    v2 = te_plus_breakdown.get("2 TE+ enhancers/set", 0)
    v3 = te_plus_breakdown.get("3 TE+ enhancers/set", 0)

    # transparent remainder so inner ring aligns to full circle
    remainder = total_sets - (v1 + v2 + v3)

    ax.pie(
        [v1, v2, v3, remainder],
        colors=[C_1, C_2, C_3, TRANSPARENT],
        startangle=startangle,
        radius=inner_radius,
        wedgeprops=dict(width=inner_width, edgecolor="black", linewidth=.3)
    )

    ax.set(aspect="equal")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

# ====== Build figure ======
fig, ax = plt.subplots(1, 1, figsize=(3.2, 3.0), dpi=DPI)
nested_donut_terminus_outer_enhancers_inner(ax)

handles = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_PLUS, markersize=8,
           label=f"terminus+ / TE+ (n={te_plus_sets})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_MINUS, markersize=8,
           label=f"terminus– / TE– (n={te_minus_sets})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_1, markersize=8,
           label=f"1 TE+ enhancer/set (n={te_plus_breakdown['1 TE+ enhancer/set']})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_2, markersize=8,
           label=f"2 TE+ enhancers/set (n={te_plus_breakdown['2 TE+ enhancers/set']})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_3, markersize=8,
           label=f"3 TE+ enhancers/set (n={te_plus_breakdown['3 TE+ enhancers/set']})"),
]

ax.legend(
    handles=handles,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(0.92, 0.5),
    prop={"family": "Arial", "size": 8},
)

plt.tight_layout()
plt.savefig(PLOT_FIGURES / "te_nested_donut_green_cooption.png", dpi=DPI, bbox_inches="tight", transparent=True)
plt.show()






import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "legend.fontsize": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 600,
    "savefig.dpi": 600,
})

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ====== Inputs ======
te_minus_sets = 420-9   # terminus– (TE–) sets
te_plus_sets  = 9    # terminus+ (TE+) sets
total_sets = te_minus_sets + te_plus_sets

te_plus_breakdown = {
    "1 TE+ enhancer/set": 7,
    "2 TE+ enhancers/set": 2,
    "3 TE+ enhancers/set": 0,
}

# ====== Appearance ======
C_PLUS  = "#1a7a6a"   # green (replaces old blue)
C_MINUS = "#E5E7E9"   # light grey

# inner ring: still green, but distinguish slices
C_1 = "#ddeea2"
C_2 = "#bcbc11"
C_3 = "#678905"

TRANSPARENT = (1, 1, 1, 0)
DPI = 600

def nested_donut_terminus_outer_enhancers_inner(ax, startangle=90):
    # matches your original ring geometry
    outer_radius = 1.0
    outer_width  = 0.32
    inner_radius = 0.68
    inner_width  = 0.28

    # ---- Outer ring: terminus – / + ----
    ax.pie(
        [te_plus_sets, te_minus_sets],
        colors=[C_PLUS, C_MINUS],
        startangle=startangle,
        radius=outer_radius,
        wedgeprops=dict(width=outer_width, edgecolor="black", linewidth=.3)
    )

    # ---- Inner ring: enhancers/set within TE+ ----
    v1 = te_plus_breakdown.get("1 TE+ enhancer/set", 0)
    v2 = te_plus_breakdown.get("2 TE+ enhancers/set", 0)
    v3 = te_plus_breakdown.get("3 TE+ enhancers/set", 0)

    # transparent remainder so inner ring aligns to full circle
    remainder = total_sets - (v1 + v2 + v3)

    ax.pie(
        [v1, v2, v3, remainder],
        colors=[C_1, C_2, C_3, TRANSPARENT],
        startangle=startangle,
        radius=inner_radius,
        wedgeprops=dict(width=inner_width, edgecolor="black",linewidth=.3)
    )

    ax.set(aspect="equal")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

# ====== Build figure ======
fig, ax = plt.subplots(1, 1, figsize=(3.2, 3.0), dpi=DPI)
nested_donut_terminus_outer_enhancers_inner(ax)

handles = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_PLUS, markersize=8,
           label=f"terminus+ / TE+ (n={te_plus_sets})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_MINUS, markersize=8,
           label=f"terminus– / TE– (n={te_minus_sets})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_1, markersize=8,
           label=f"1 TE+ enhancer/set (n={te_plus_breakdown['1 TE+ enhancer/set']})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_2, markersize=8,
           label=f"2 TE+ enhancers/set (n={te_plus_breakdown['2 TE+ enhancers/set']})"),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=C_3, markersize=8,
           label=f"3 TE+ enhancers/set (n={te_plus_breakdown['3 TE+ enhancers/set']})"),
]

ax.legend(
    handles=handles,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(0.92, 0.5),
    prop={"family": "Arial", "size": 8},
)

plt.tight_layout()
plt.savefig(PLOT_FIGURES / "te_nested_donut_splitting.png", dpi=DPI, bbox_inches="tight", transparent=True)
plt.show()
