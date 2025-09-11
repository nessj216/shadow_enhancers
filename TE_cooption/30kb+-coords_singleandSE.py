import pandas as pd

# Load the file
file_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/011925_all_shadowsets_DM6.bed"
df = pd.read_csv(file_path, sep="\t", header=None)

# Load chromosome sizes
chrom_size_file = "/Users/jillianness/Desktop/TOBIAS091922/general/dm6.fasta.fai"
chrom_sizes = pd.read_csv(chrom_size_file, sep="\t", header=None, usecols=[0, 1], names=["chrom", "size"])
chrom_size_dict = chrom_sizes.set_index("chrom")["size"].to_dict()

# Ensure columns are numeric
df[1] = pd.to_numeric(df[1])
df[2] = pd.to_numeric(df[2])

# Modify columns 2 and 3
for i in range(len(df)):
    chrom = df.loc[i, 0]
    chrom_max = chrom_size_dict.get(chrom, float('inf'))

    col_min = min(df.loc[i, 1], df.loc[i, 2])
    col_max = max(df.loc[i, 1], df.loc[i, 2])

    new_min = max(10, col_min - 30000)
    new_max = min(chrom_max, col_max + 30000)

    # Assign back ensuring order is maintained
    if df.loc[i, 1] < df.loc[i, 2]:
        df.loc[i, 1], df.loc[i, 2] = new_min, new_max
    else:
        df.loc[i, 2], df.loc[i, 1] = new_min, new_max

# Save the modified file
output_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_landscape/TE_envi_singles/30kb_plusminus_shadows_enhancer.txt"
df.to_csv(output_path, sep="\t", header=False, index=False)

print(f"Modified file saved as: {output_path}")
