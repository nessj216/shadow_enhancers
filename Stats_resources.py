
#!/usr/bin/env python3
"""
nonparam_tests.py

A single script running all our nonparametric tests.
"""


'''Perform Wilcoxon signed‐rank test'''

from scipy.stats import wilcoxon

# Example paired data (before vs. after treatment)
before = [5.0, 4.2, 5.8, 6.1, 4.9]
after  = [5.8, 4.9, 6.0, 6.5, 5.3]

# Perform Wilcoxon signed‐rank test
# zero_method='wilcox' handles ties/zeros appropriately
# alternative='two-sided' tests for any change in median
stat, p = wilcoxon(before, after, zero_method='wilcox', alternative='two-sided')

# Print the W statistic and two‐sided p‐value
print(f"W = {stat:.3f}, p = {p:.3f}")




'''Perform two-sided Mann–Whitney U test'''

from scipy.stats import mannwhitneyu

# Sample data
group_A = [5.1, 7.3, 6.8, 5.5, 7.0] #swap out for name of data vector
group_B = [8.2, 7.9, 9.1, 8.5, 8.8] #swap out for name of data vector

# Perform two-sided Mann–Whitney U test
u_stat, p_value = mannwhitneyu(group_A, group_B, alternative='two-sided')

print(f"U statistic: {u_stat:.2f}")
print(f"p-value: {p_value:.3f}")

# Interpretation
alpha = 0.05
if p_value < alpha:
    print("→ Reject H0: distributions differ")
else:
    print("→ Fail to reject H0: no evidence of difference")


'''Perform Kruskal–Wallis test '''

from scipy.stats import kruskal

# Example independent groups (three samples of ordinal/continuous data)
group1 = [5.1, 4.9, 5.7, 5.3, 5.4]
group2 = [6.2, 6.0, 5.8, 6.3, 6.1]
group3 = [7.5, 7.2, 7.8, 7.4, 7.3]

# Perform Kruskal–Wallis H Test on the three groups
stat, p = kruskal(group1, group2, group3)

# Print the H statistic and two‐sided p‐value
print(f"H = {stat:.3f}, p = {p:.3f}")




#plotting

import matplotlib.pyplot as plt

# Sample data for three independent groups
group1 = [5.1, 4.9, 5.7, 5.3, 5.4]
group2 = [6.2, 6.0, 5.8, 6.3, 6.1]
group3 = [7.5, 7.2, 7.8, 7.4, 7.3]

# 2) Violin Plot

plt.figure()
plt.violinplot([group1, group2, group3], showmeans=False, showmedians=True)
plt.xticks([1, 2, 3], ['Group 1', 'Group 2', 'Group 3'])
plt.ylabel('Value')
plt.title('Kruskal–Wallis Data: Violin Plot (Closer Medians)')
plt.show()


# 1) Box Plot
plt.figure()
plt.boxplot([group1, group2, group3])
plt.xticks([1, 2, 3], ['Group 1', 'Group 2', 'Group 3'])
plt.ylabel('Value')
plt.title('Kruskal–Wallis Data: Box Plot')
plt.show()
