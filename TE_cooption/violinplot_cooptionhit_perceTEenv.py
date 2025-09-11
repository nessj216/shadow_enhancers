import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#+-30 kb TE % environment
#for plotting the TE env of single or shadow enhancers that are TE co-option hits
# (violinplot_perceTE file is for ALL shadows or single enhancer TE environemtns
# File paths for the two CSV files
singles_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_envi_singles/singles_TEhits_env/TEpercent_singles_Deduplicated_Data.csv"
shadows_file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_env_shadows/TEenv_shadowcooptionhits/Deduplicated_File_shadows.csv"

# Load the single enhancer file
singles_df = pd.read_csv(singles_file_path, sep=",", header=None)
# Extract the relevant column (column index 9, which is the 10th column)
singles_df = singles_df[[9]].copy()
singles_df.columns = ["Tenth_Column"]
singles_df["Category"] = "Single Enhancer"

# Load the shadow enhancer file
shadows_df = pd.read_csv(shadows_file_path, sep=",", header=0)  # Assuming the file has a header
# Extract the relevant column (column index 9, which is the 10th column)
shadows_df = shadows_df.iloc[:, [9]].copy()
shadows_df.columns = ["Tenth_Column"]
shadows_df["Category"] = "Shadow Enhancer"

# Combine the two datasets
combined_df = pd.concat([singles_df, shadows_df], ignore_index=True)

# Convert to numeric to avoid any plotting errors
combined_df["Tenth_Column"] = pd.to_numeric(combined_df["Tenth_Column"], errors="coerce")

# Create the violin plot
plt.figure(figsize=(8, 6))
sns.violinplot(x="Category", y="Tenth_Column", data=combined_df, inner="quartile")

# Customize the plot
plt.xlabel("Enhancer Type")
plt.ylabel("Tenth Column Values")
plt.title("Violin Plot of 10th Column in Single and Shadow Enhancers")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Ensure all values are properly numeric and drop NaN values
singles_df["Tenth_Column"] = pd.to_numeric(singles_df["Tenth_Column"], errors="coerce")
shadows_df["Tenth_Column"] = pd.to_numeric(shadows_df["Tenth_Column"], errors="coerce")

# Drop NaN values if any remain
singles_df = singles_df.dropna()
shadows_df = shadows_df.dropna()

# Check data types to confirm conversion
print(singles_df.dtypes)
print(shadows_df.dtypes)

# Plot the distributions using a histogram
plt.figure(figsize=(8, 6))

sns.histplot(singles_df["Tenth_Column"], label="Single Enhancer", kde=True, bins=30, alpha=0.5)
sns.histplot(shadows_df["Tenth_Column"], label="Shadow Enhancer", kde=True, bins=30, alpha=0.5)

# Customize plot
plt.xlabel("Tenth Column Values")
plt.ylabel("Frequency")
plt.title("Distribution of 10th Column in Single and Shadow Enhancers")
plt.legend()
plt.tight_layout()

# Show the plot
plt.show()

