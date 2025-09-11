import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Sample Data
data = {
    'Class': ['Class I', 'Class I', 'Class I', 'Class II', 'Class II', 'Class II'],
    'Category': ['Genome', 'Single Enhancer', 'Shadow Enhancer', 'Genome', 'Single Enhancer', 'Shadow Enhancer'],
    'Percentage': [88, 71, 80, 12, 29, 20]
}

df = pd.DataFrame(data)

# Compute fold-change relative to the genome baseline
df_baseline = df[df["Category"] == "Genome"].set_index("Class")["Percentage"]
df["Fold Change"] = df.apply(lambda row: row["Percentage"] / df_baseline[row["Class"]]
                             if row["Category"] != "Genome" else 1, axis=1)

# Filter out the genome baseline from the plot since it will always be 1
df_plot = df[df["Category"] != "Genome"]

# Set up figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

# Plot for Class I (default colors)
sns.barplot(
    data=df_plot[df_plot["Class"] == "Class I"],
    x="Category",
    y="Fold Change",
    ax=axes[0]
)
axes[0].set_title("Class I Retroelements", fontsize=18)
axes[0].set_xlabel("Enhancer Type", fontsize=16)
axes[0].set_ylabel("Fold Change Over Genome Baseline", fontsize=16)
axes[0].axhline(1, color='black', linestyle='--', linewidth=1)  # Baseline at 1

# Plot for Class II (default colors)
sns.barplot(
    data=df_plot[df_plot["Class"] == "Class II"],
    x="Category",
    y="Fold Change",
    ax=axes[1]
)
axes[1].set_title("Class II DNA Elements", fontsize=18)
axes[1].set_xlabel("Enhancer Type", fontsize=18)
axes[1].axhline(1, color='black', linestyle='--', linewidth=1)  # Baseline at 1
axes[0].tick_params(axis='x', labelsize=18)
axes[1].tick_params(axis='x', labelsize=18)
axes[0].tick_params(axis='y', labelsize=16)
axes[1].tick_params(axis='y', labelsize=16)
# Adjust layout
plt.tight_layout()
plt.show()
