# Determine the maximum y-value for annotations, slightly decrease the starting point
y_max = max([item for sublist in data_list for item in sublist])
y_text = y_max * 1.01  # Starting point for text annotations, decreased from 1.05 to 1.01

bar_height = y_max * 0.02  # Decreased height of the significance bar

# Loop over all combinations
for (idx1, idx2) in combinations:
    _, p_value = mannwhitneyu(data_list[idx1], data_list[idx2])
    # Apply Bonferroni correction
    corrected_p_value = p_value * num_tests

    # Check for significance and set the appropriate symbol
    if corrected_p_value < 0.001:
        sig_symbol = '***'  # Highly significant
    elif corrected_p_value < 0.01:
        sig_symbol = '**'  # Very significant
    elif corrected_p_value < 0.05:
        sig_symbol = '*'  # Significant
    else:
        sig_symbol = 'n.s.'  # Not significant, can be omitted or adjusted as needed

    # Annotate if corrected p-value is less than 0.05 (common significance level)
    if corrected_p_value < 0.05:
        # Draw a line between the pairs and a star for significance
        plt.plot([idx1, idx2], [y_text, y_text], color='black')
        plt.text((idx1 + idx2) / 2, y_text + bar_height / 2, sig_symbol, ha='center', va='bottom')
        y_text += bar_height * 1.5  # Adjusted spacing to decrease the height increment between annotations
