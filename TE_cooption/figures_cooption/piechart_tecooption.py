import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# Use a clean style with default fonts and update font size globally
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({'font.size': 18})


def plot_te_sets_custom_striped():
    # ---- Data for the Pie Chart ----
    te_minus_sets = 376  # TE– sets (approx. 90%)
    te_plus_sets = 44  # TE+ sets (approx. 10%)
    total_sets = te_minus_sets + te_plus_sets
    te_minus_pct = te_minus_sets / total_sets * 100
    te_plus_pct = te_plus_sets / total_sets * 100

    pie_labels = [
        f'TE– sets\n({te_minus_pct:.1f}%,{te_minus_sets} sets)',
        f'TE+ sets\n({te_plus_pct:.1f}%,{te_plus_sets} sets)'
    ]
    # TE– slice in grey; TE+ slice will be overwritten with a stripe pattern.
    pie_colors = ['#d3d3d3', '#66c2a5']

    # ---- Data for the Stacked Bar (TE+ breakdown) ----
    te_plus_breakdown = dict(sorted({
                                        '1 TE+ enhancer/set': 37,
                                        '2 TE+ enhancers/set': 5,
                                        '3 TE+ enhancers/set': 1,
                                    }.items(), key=lambda x: x[1]))

    total_TE_plus_breakdown = sum(te_plus_breakdown.values())
    # Use different shades of green for stacked bar segments
    breakdown_colors = ['#31a354','#74c476','#c7e9c0']

    # ---- Create the Figure and Axes ----
    fig = plt.figure(figsize=(12, 7))

    # Pie chart axes
    ax_pie = fig.add_axes([0.05, 0.1, 0.4, 0.8])
    explode = [0, 0.09]  # Explode TE+ slice for emphasis

    # Draw the pie. (TE+ slice will later be overlaid with a stripe pattern.)
    wedges, texts, autotexts = ax_pie.pie(
        [te_minus_sets, te_plus_sets],
        labels=pie_labels,
        colors=pie_colors,
        startangle=20,
        autopct='',
        explode=explode,

    )
    #ax_pie.set_title('TE Sets Composition', fontsize=18)

    # ---- Create Stripe Pattern for TE+ Slice ----
    # Get the TE+ wedge (second wedge)
    te_plus_wedge = wedges[1]
    # Remove its facecolor so the stripe image shows through
    te_plus_wedge.set_facecolor('none')

    # Generate a stripe pattern image: diagonal stripes with colors '#a1d99b' and '#31a354'
    height, width = 200, 200  # Resolution of the pattern image
    stripe_width = 10  # Adjust stripe thickness (in pixels)

    # Create indices and a mask for alternating stripes
    i = np.arange(height).reshape(-1, 1)
    j = np.arange(width).reshape(1, -1)
    mask = ((i + j) // stripe_width) % 2 == 0
    pattern = np.empty((height, width, 3))
    pattern[mask] = mcolors.to_rgb('#a1d99b')
    pattern[~mask] = mcolors.to_rgb('#31a354')

    # Overlay the stripe pattern image over the entire pie area.
    # The extent here (-1 to 1) matches the data coordinates of the pie (with center at (0,0) and radius 1).
    im = ax_pie.imshow(pattern, extent=[-1, 1, -1, 1], origin='lower', interpolation='nearest',
                       zorder=te_plus_wedge.get_zorder() - 1)
    # Clip the stripe image to the shape of the TE+ wedge so that stripes appear only in that slice.
    im.set_clip_path(te_plus_wedge)

    # ---- Create Stacked Bar Chart for TE+ Breakdown ----
    # ---- Create Stacked Bar Chart for TE+ Breakdown ----
    # ---- Create Stacked Bar Chart for TE+ Breakdown ----
    ax_bar = fig.add_axes([0.65, 0.35, 0.15, 0.3])  # Narrower and tighter position

    bottom_val = 0
    bar_x = 0  # x-position of the bar
    bar_width = 0.1  # Narrow bar width

    for (label, value), color in zip(te_plus_breakdown.items(), breakdown_colors):
        bar_height = value / total_TE_plus_breakdown * 100
        ax_bar.bar(bar_x, bar_height, width=bar_width, bottom=bottom_val,
                   color=color, edgecolor='black', label=label)
        y_mid = bottom_val + bar_height / 2
        # Shift text slightly to right of bar if short, otherwise center
        if bar_height < 5:
            ax_bar.text(bar_x + bar_width/2 +.01, y_mid, f'{bar_height:.1f}%',
                        ha='left', va='center', color='black')
        else:
            ax_bar.text(bar_x, y_mid, f'{bar_height:.1f}%',
                        ha='center', va='center', color='black', fontsize=18)
        bottom_val += bar_height

    ax_bar.set_xticks([bar_x])
    ax_bar.set_xticklabels(['TE+ Sets'], )
    ax_bar.set_ylabel('Percent of TE+ Sets')
    ax_bar.set_ylim(0, 100)
    #ax_bar.set_title('Breakdown of TE+ Sets', fontsize=18)
    # Move legend outside bar axes
    ax_bar.legend(frameon=True, loc='upper left', bbox_to_anchor=(1.05, 1.0))
    ax_bar.set_xticks([])  # Remove x-axis ticks
    ax_bar.set_yticks([])  # Remove y-axis ticks
    ax_bar.spines['top'].set_visible(False)  # Hide top border
    ax_bar.spines['right'].set_visible(False)  # Hide right border
    ax_bar.spines['left'].set_visible(False)  # Hide left border
    ax_bar.spines['bottom'].set_visible(False)  # Hide bottom border
    ax_bar.set_xlabel("")  # Remove x-axis label if set
    ax_bar.set_ylabel("")  # Remove y-axis label if set
    ax_bar.grid(False)
    ax_bar.hlines(y=0, xmin=bar_x - bar_width / 2, xmax=bar_x + bar_width / 2, colors='black', linewidth=1.5)

    #plt.tight_layout()
    # Save the plot as a high-resolution image
    plt.savefig("/Users/jillianness/Desktop/Figures_shadowbirth/piebarCooption.png", dpi=600, bbox_inches='tight', pad_inches=0.1)
    plt.show()

# Example usage:
plot_te_sets_custom_striped()
