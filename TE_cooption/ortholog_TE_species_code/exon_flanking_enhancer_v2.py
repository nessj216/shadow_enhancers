import pandas as pd
import subprocess

# ── Step 1: Load your enhancer file (with extra columns!) ──
enhancer_df = pd.read_csv(
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/single/sort_single_TEoutput.bed",
    sep="\t",
    header=None
)
# Create a key on the chr/start–end so we can merge later:
enhancer_df["key"] = (
    enhancer_df[0].astype(str)
    + ":"
    + enhancer_df[1].astype(str)
    + "-"
    + enhancer_df[2].astype(str)
)

# ── Step 2: Write out the full enhancer table for bedtools ──
enhancer_bed = "all_enhancers.bed"
enhancer_df.drop(columns="key").to_csv(
    enhancer_bed, sep="\t", index=False, header=False
)
# ── (← instead of subsetting to [0,1,2], we write out _all_ columns) ──

# ── Step 3: Prepare the exon file as before ──
exons = pd.read_csv(
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/old_dont_use_flanking_regions_find_orthologs/redo_correct/dm6_ncbiRefSeq_exons_sort.bed",
    sep="\t",
    header=None
)
exon_bed = "exon_unique.bed"
exons[[0,1,2]].drop_duplicates().to_csv(
    exon_bed, sep="\t", index=False, header=False
)

# ── Step 4: Run bedtools closest upstream (5′) and downstream (3′) ──
file5p = "closest_5p.bed"
file3p = "closest_3p.bed"
subprocess.run(
    f"bedtools closest -a {enhancer_bed} -b {exon_bed} -D a -id -io > {file5p}",
    shell=True
)
subprocess.run(
    f"bedtools closest -a {enhancer_bed} -b {exon_bed} -D a -iu -io > {file3p}",
    shell=True
)

# ── Step 5: Read the bedtools outputs ──
# Figure out how many columns your enhancer file had:
num_enh_cols = enhancer_df.shape[1] - 1   # minus the "key" column
# (so if your file was chr, start, end, name, score → num_enh_cols = 5)

# Load the 5′ results
# ── After running bedtools, load the 5′ results and drop the last (distance) column ──
five = pd.read_csv(file5p, sep="\t", header=None)
five = five.iloc[:, :-1]    # drop the very last column
five.columns = (
    list(range(num_enh_cols))  # your original enhancer columns: 0,1,2,…
    + ["five_chr", "five_start", "five_end"]
)

# ── Same for the 3′ results ──
three = pd.read_csv(file3p, sep="\t", header=None)
three = three.iloc[:, :-1]   # drop the distance column here too
three.columns = (
    list(range(num_enh_cols))
    + ["three_chr", "three_start", "three_end"]
)


# ── Step 6: Re‐create the same key on each ──
five["key"] = (
    five[0].astype(str) + ":" + five[1].astype(str) + "-" + five[2].astype(str)
)
three["key"] = (
    three[0].astype(str) + ":" + three[1].astype(str) + "-" + three[2].astype(str)
)

# ── Step 7: Merge them back _onto_ your original enhancer_df ──
# Drop any dedup calls: we want to keep every row as it was
# (so we do _not_ do five.drop_duplicates(...))

# First, merge the 5′ info (left‐join to keep all enhancers)
merged = pd.merge(
    enhancer_df,
    five[["key", "five_chr", "five_start", "five_end"]],
    on="key",
    how="left"
)

# Then merge the 3′ info
merged = pd.merge(
    merged,
    three[["key", "three_chr", "three_start", "three_end"]],
    on="key",
    how="left"
)

# ── Step 8: Reorder and write out ──
# Say you want: chr,  five_end,  three_start,  [all enhancer cols],  enhancer_start, enhancer_end
output_cols = [
    "five_chr", "five_end", "three_chr", "three_start",
    0, 1, 2    # these are your original chr, start, end
]
# plus any other extra col indices (e.g. 3,4,5 ...)
# and note: .iloc needs integer positions, .loc can use labels if we rename
# For simplicity let’s do .loc with column names:

final = merged.loc[
    :,
    ["five_chr", "five_end", "three_start"]
    + list(range(num_enh_cols))
]

# Write it
final.to_csv(
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/single/single_enhancers_with_flanks_and_extras.bed",
    sep="\t",
    index=False,
    header=False
)

import pandas as pd

# 1) Load your file
#    Make sure the path matches wherever you wrote enhancers_with_flanks_and_extras.bed
df = pd.read_csv(
     "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/single/single_enhancers_with_flanks_and_extras.bed",
    sep="\t",
    header=None
)

# 2) Define which zero‐based columns hold your (chr, start, end) triples
#    - enhancer coords are in 1-based cols 4,5,6 → zero‐based [3,4,5]
#    - your “set” coords are in 1-based cols 8,9,10 → zero‐based [7,8,9]
enh_cols = [3, 4, 5]
set_cols = [7, 8, 9]

# 3) Helper to turn 0→A, 1→B, … 25→Z, 26→AA, etc.
def num_to_letters(n: int) -> str:
    s = ""
    while True:
        n, rem = divmod(n, 26)
        s = chr(rem + 65) + s
        if n == 0:
            break
        n -= 1
    return s

# 4) Factorize the enhancer triples → numeric codes 0,1,2…
enh_tuples = df.iloc[:, enh_cols].apply(tuple, axis=1)
enh_codes, enh_uniques = pd.factorize(enh_tuples)

# Map those codes to A, B, … etc.
alpha_ids = [num_to_letters(code) for code in enh_codes]
df["ENH_ID"] = alpha_ids

# 5) Factorize the “set” triples → numeric codes 0,1,2…
set_tuples = df.iloc[:, set_cols].apply(tuple, axis=1)
set_codes, set_uniques = pd.factorize(set_tuples)

# Shift to 1‐based and attach
df["SET_ID"] = set_codes + 1

# 6) Write out the augmented file
df.to_csv(
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/single/single_enhancers_with_flanks__withIDs.bed",
    sep="\t",
    header=False,
    index=False
)


import pandas as pd
import subprocess

# … [all of your existing code up through writing enhancers_with_flanks_and_extras_withIDs.bed] …

# ── Step 9: Load your ID’d BED and write a subset ──
id_file = (
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/single/single_enhancers_with_flanks__withIDs.bed"
)
subset_file = (
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/single/inputfile_SINGLEenhancers_with_flanks_and_extras_withIDs.bed"
)

# Read in with header=None so columns are 0,1,2,…
df_ids = pd.read_csv(id_file, sep="\t", header=None, dtype=str)

# Define the zero‐based columns to extract
# (1,2,3,11,13,14,13 in 1-based → 0,1,2,10,12,13,12 in 0-based)
cols_to_keep = [0, 1, 2, 10, 12, 13, 12]

# Slice and write
df_ids.iloc[:, cols_to_keep].to_csv(
    subset_file,
    sep="\t",
    header=False,
    index=False
)
