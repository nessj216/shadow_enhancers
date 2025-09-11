#+-30 kb TE % environment

# (violinplot_perceTE file is for ALL shadows or single enhancer TE environemtns

#violinplot_cooptionhits_perceTEenv file is for plotting the TE env of single or shadow enhancers that are TE co-option hits


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths for the two files
singles_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_envi_singles/v3_merged_single_env__TE_overlap.bed"
shadows_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_env_shadows/Final_percentsums_merged_shadows.txt"


# Load the single enhancer file
singles_df = pd.read_csv(singles_file_path, sep="\t", header=None)
# Extract the relevant column (column index 4, which is the 5th column)
singles_df = singles_df[[4]].copy()
singles_df.columns = ["Percent_TE_Overlap"]
singles_df["Category"] = "Single Enhancer"

# Load the shadow enhancer file (has a header)
shadows_df = pd.read_csv(shadows_file_path, sep="\t", header=0)
# Extract the relevant column (column index 4, which is the 5th column)
shadows_df = shadows_df.iloc[:, [4]].copy()
shadows_df.columns = ["Percent_TE_Overlap"]
shadows_df["Category"] = "Shadow Enhancer"

# Combine the two datasets
combined_df = pd.concat([singles_df, shadows_df], ignore_index=True)

# Convert to numeric to avoid any plotting errors
combined_df["Percent_TE_Overlap"] = pd.to_numeric(combined_df["Percent_TE_Overlap"], errors="coerce")

# Create the violin plot
plt.figure(figsize=(8, 6))
sns.violinplot(x="Category", y="Percent_TE_Overlap", data=combined_df, inner="quartile")

# Customize the plot
plt.xlabel("Enhancer Type")
plt.ylabel("Percent TE Overlap")
plt.title("Violin Plot of TE Overlap in Single and Shadow Enhancers")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Ensure all values are properly numeric and drop NaN values
singles_df["Percent_TE_Overlap"] = pd.to_numeric(singles_df["Percent_TE_Overlap"], errors="coerce")
shadows_df["Percent_TE_Overlap"] = pd.to_numeric(shadows_df["Percent_TE_Overlap"], errors="coerce")

# Drop NaN values if any remain
singles_df = singles_df.dropna()
shadows_df = shadows_df.dropna()

# Check data types to confirm conversion
print(singles_df.dtypes)
print(shadows_df.dtypes)

# Plot the distributions using histogram (to avoid kernel density estimation issues)
plt.figure(figsize=(8, 6))

sns.histplot(singles_df["Percent_TE_Overlap"], label="Single Enhancer", kde=True, bins=30, alpha=0.5)
sns.histplot(shadows_df["Percent_TE_Overlap"], label="Shadow Enhancer", kde=True, bins=30, alpha=0.5)

# Customize plot
plt.xlabel("Percent TE Overlap")
plt.ylabel("Frequency")
plt.title("Distribution of TE Overlap in Single and Shadow Enhancers")
plt.legend()
plt.tight_layout()

# Show the plot
plt.show()


