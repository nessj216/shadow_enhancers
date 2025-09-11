import matplotlib.pyplot as plt

# Apply seaborn style
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif'
})

# Color palette
pretty_colors = ['#e6e6e6', '#c3aed6', '#fff3b0']

def draw_custom_scaled_circles(ax, outer_val, outer_label,
                                middle_val=None, middle_label=None,
                                inner_val=None, inner_label=None,
                                colors=None,
                                scale_factor=0.5,
                                full_radius=0.4,
                                outer_text_y_offset=.05,
                                middle_text_y_offset=0.02,
                                inner_text_y_offset=-0.05,
                                label_fontsize=14,
                                external_middle_label=True,
                                external_inner_label=True,
                                extra_label_text=None,
                                extra_label_pos=(0.95, 0.00)):
    if colors is None:
        colors = pretty_colors

    # Radii
    R_outer = full_radius * scale_factor
    R_middle = R_outer * 0.5
    R_inner = R_middle * 0.1

    # Y centers
    Y_outer = R_outer
    Y_middle = R_middle
    Y_inner = R_inner

    # Outer circle
    ax.add_patch(plt.Circle((0.5, Y_outer), R_outer, color=colors[0], ec='white', linewidth=2))
    ax.text(0.5, Y_outer+ outer_text_y_offset, f"{outer_label}\n(n={outer_val})",
            ha='center', va='center', fontsize=label_fontsize)

    # Middle circle
    if middle_val is not None:
        percent_middle = f"{(middle_val / outer_val * 100):.1f}%" if outer_val else ""
        ax.add_patch(plt.Circle((0.5, Y_middle), R_middle, color=colors[1], ec='white', linewidth=2))
        if external_middle_label:
            mid_line_x = 0.5 + R_middle + 0.02
            ax.plot([0.5 + R_middle * 0.7, mid_line_x], [Y_middle, Y_middle + middle_text_y_offset], color='black', linewidth=1)
            ax.text(mid_line_x + 0.02, Y_middle + middle_text_y_offset,
                    f"{middle_label}\n{percent_middle}\n(n={middle_val})",
                    va='center', ha='left', fontsize=label_fontsize)
        else:
            ax.text(0.5, Y_middle + middle_text_y_offset, f"{middle_label}\n{percent_middle}\n(n={middle_val})",
                    ha='center', va='center', fontsize=label_fontsize)

    # Inner circle
    if inner_val is not None:
        percent_inner = f"{(inner_val / middle_val * 100):.1f}%" if middle_val else ""
        ax.add_patch(plt.Circle((0.5, Y_inner), R_inner, color=colors[2], ec='white', linewidth=2))
        if external_inner_label:
            in_line_x = 0.5 + R_inner + 0.02
            ax.plot([0.5 + R_inner * 0.7, in_line_x], [Y_inner, Y_inner + inner_text_y_offset], color='black', linewidth=1)
            ax.text(in_line_x + 0.02, Y_inner + inner_text_y_offset,
                    f"{inner_label}\n{percent_inner}\n(n={inner_val})",
                    va='center', ha='left', fontsize=label_fontsize)
        else:
            ax.text(0.5, Y_inner + inner_text_y_offset, f"{inner_label}\n{percent_inner}\n(n={inner_val})",
                    ha='center', va='center', fontsize=label_fontsize)

    # Extra label for 0% case
    if extra_label_text:
        ax.text(extra_label_pos[0], extra_label_pos[1], extra_label_text,
                ha='center', va='center', fontsize=label_fontsize)

    ax.set_xlim(0, 1.25)
    ax.set_ylim(0, full_radius + 0.1)
    ax.set_aspect('equal')
    ax.axis('off')

# Create plots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: Shadow TE with extra label for D. pseudo (0%)
draw_custom_scaled_circles(
    axes[0],
    outer_val=46,
    outer_label="TE+ Shadows",
    middle_val=15,
    middle_label=r"$\it{D.\ sim}$",
    scale_factor=0.5,
    outer_text_y_offset=0.05,
    external_middle_label=True,
    external_inner_label=False,
    extra_label_text=r"$\it{D.\ pseudo}$" + "\n0%",
    extra_label_pos=(0.65, -.03)
)

# Right: Single TE with both external labels
draw_custom_scaled_circles(
    axes[1],
    outer_val=915,
    outer_label="TE+ Singles",
    middle_val=512,
    middle_label=r"$\it{D.\ sim}$",
    inner_val=5,
    inner_label=r"$\it{D.\ pseudo}$",
    scale_factor=1.0,
    middle_text_y_offset=0.05,
    external_middle_label=True,
    external_inner_label=True
)

plt.tight_layout()
fig.subplots_adjust(wspace=0.3)
plt.savefig("/Users/jillianness/Desktop/Figures_shadowbirth/TEcoooption_species.png", dpi=600, bbox_inches='tight', pad_inches=0.1)

plt.show()
