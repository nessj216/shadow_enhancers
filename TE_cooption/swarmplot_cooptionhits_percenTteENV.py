'''import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths for the two CSV files
singles_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_envi_singles/singles_TEhits_env/TEpercent_singles_Deduplicated_Data.csv"
shadows_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_env_shadows/TEenv_shadowcooptionhits/Deduplicated_File_shadows.csv"

singles_df = pd.read_csv(singles_file_path, sep=",", header=None)
# Extract the relevant column (column index 9, which is the 10th column)
singles_df = singles_df[[9]].copy()
singles_df.columns = ["Tenth_Column"]
singles_df["Category"] = "Single Enhancer"

# Load the shadow enhancer file
shadows_df = pd.read_csv(shadows_file_path, sep=",", header=0)
# Extract the relevant column (column index 9, which is the 10th column)
shadows_df = shadows_df.iloc[:, [9]].copy()
shadows_df.columns = ["Tenth_Column"]
shadows_df["Category"] = "Shadow Enhancer"

# Combine the two datasets
combined_df = pd.concat([singles_df, shadows_df], ignore_index=True)

# Convert to numeric and drop NaNs
combined_df["Tenth_Column"] = pd.to_numeric(combined_df["Tenth_Column"], errors="coerce")
combined_df.dropna(subset=["Tenth_Column"], inplace=True)

# -- Box Plot Only (with quartiles, median, and mean line) --
plt.figure(figsize=(8, 6))

sns.boxplot(
    x="Category",
    y="Tenth_Column",
    data=combined_df,
    showmeans=True,     # Include a line or marker for the mean
    meanline=True,      # Draw the mean as a line (instead of a point)
    width=0.3
)

# Customize the plot - bigger labels, normal (horizontal) x-axis
plt.xlabel("Enhancer Type", fontsize=14)
plt.ylabel("% TE local landscape", fontsize=14)
plt.title("Box Plot of 10th Column in Single and Shadow Enhancers", fontsize=16)
plt.xticks(rotation=0, fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.show()

# -- (Optional) Histograms for each group --
singles_df["Tenth_Column"] = pd.to_numeric(singles_df["Tenth_Column"], errors="coerce").dropna()
shadows_df["Tenth_Column"] = pd.to_numeric(shadows_df["Tenth_Column"], errors="coerce").dropna()

plt.figure(figsize=(8, 6))
sns.histplot(singles_df["% TE local landscape"], label="Single Enhancer", kde=True, bins=30, alpha=0.5)
sns.histplot(shadows_df["% TE local landscape"], label="Shadow Enhancer", kde=True, bins=30, alpha=0.5)

plt.xlabel("% TE local landscape", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.title("Distribution of 10th Column in Single and Shadow Enhancers", fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.show()
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths for the two TXT files (tab-delimited assumed)
singles_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_envi_singles/v3_merged_single_env__TE_overlap.bed"
shadows_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_env_shadows/Final_percentsums_merged_shadows.txt"


# ---------------------------------------------------------
# 1. Read in SINGLE Enhancer file
#    (Assuming NO header row, we will manually name columns)
# ---------------------------------------------------------
singles_df = pd.read_csv(
    singles_file_path,
    sep="\t",      # Change if your .txt is delimited differently (e.g., sep=" " or sep=",")
    header=None    # No header row in this file
)

# Extract the relevant column (column index 9 => the 10th column)
singles_df = singles_df[[4]].copy()
singles_df.columns = ["Tenth_Column"]
singles_df["Category"] = "Single Enhancer"

# ---------------------------------------------------------
# 2. Read in SHADOW Enhancer file
#    (If the second .txt file actually has a header row,
#     use header=0. Otherwise, use header=None.)
# ---------------------------------------------------------
shadows_df = pd.read_csv(
    shadows_file_path,
    sep="\t",    # Change if needed
    header=None  # or header=0 if the file has a header
)

# Extract the relevant column (column index 9 => the 10th column)
shadows_df = shadows_df.iloc[:, [4]].copy()
shadows_df.columns = ["Tenth_Column"]
shadows_df["Category"] = "Shadow Enhancer"

# ---------------------------------------------------------
# 3. Combine and clean data
# ---------------------------------------------------------
combined_df = pd.concat([singles_df, shadows_df], ignore_index=True)
combined_df["Tenth_Column"] = pd.to_numeric(combined_df["Tenth_Column"], errors="coerce")
combined_df.dropna(subset=["Tenth_Column"], inplace=True)

# ---------------------------------------------------------
# 4. Box Plot
# ---------------------------------------------------------
plt.figure(figsize=(8, 6))

sns.boxplot(
    x="Category",
    y="Tenth_Column",
    data=combined_df,
    showmeans=False,   # Include a line or marker for the mean
    meanline=False,    # Draw the mean as a line (instead of a point)
    width=0.3
)

# Customize fonts and orientation
plt.xlabel("Enhancer Type", fontsize=18)
plt.ylabel("% TE local landscape", fontsize=1)
plt.title("Box Plot of 10th Column in Single and Shadow Enhancers", fontsize=16)
plt.xticks(rotation=0, fontsize=12)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.show()

# ---------------------------------------------------------
# 5. Histograms (optional)
# ---------------------------------------------------------
# Convert each group again and drop NaNs
# (They may already be numeric, but this ensures consistency)
singles_df["Tenth_Column"] = pd.to_numeric(singles_df["Tenth_Column"], errors="coerce").dropna()
shadows_df["Tenth_Column"] = pd.to_numeric(shadows_df["Tenth_Column"], errors="coerce").dropna()

plt.figure(figsize=(8, 6))

sns.histplot(singles_df["Tenth_Column"], label="Single Enhancer", kde=True, bins=30, alpha=0.5)
sns.histplot(shadows_df["Tenth_Column"], label="Shadow Enhancer", kde=True, bins=30, alpha=0.5)

plt.xlabel("% TE local landscape", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.title("Distribution of 10th Column in Single and Shadow Enhancers", fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.tight_layout()
plt.show()
