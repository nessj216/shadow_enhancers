
import pandas as pd
import subprocess

# Step 1: Load the BED file (no quotes inside the string path!)
enhancer_bed = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/REdoing_TE merge/final_merged_cleaned_TE.bed", sep="\t", header=None)

# Step 2: Extract columns 1, 2, 3 and save to new BED file
enhancer_bed[[0, 1, 2]].to_csv("cleaned_enhancer_coords.bed", sep="\t", header=False, index=False)


# Step 2: Load and deduplicate exon BED (keep only chr/start/end)
exons = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/old_dont_use_flanking_regions_find_orthologs/redo_correct/dm6_ncbiRefSeq_exons_sort.bed", sep="\t", header=None)
exon_minimal = exons[[0, 1, 2]].drop_duplicates()
exon_bed = "exon_unique.bed"
exon_minimal.to_csv(exon_bed, sep="\t", header=False, index=False)
file3p="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/old_dont_use_flanking_regions_find_orthologs/redo_correct/closest_3p.bed"
file5p="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/old_dont_use_flanking_regions_find_orthologs/redo_correct/closest_5p.bed"
# Step 3: Run bedtools to find closest upstream and downstream exons
subprocess.run(f"bedtools closest -a {enhancer_bed} -b {exon_bed} -D a -id > {file5p}", shell=True)
subprocess.run(f"bedtools closest -a {enhancer_bed} -b {exon_bed} -D a -iu > {file3p}", shell=True)

# Step 4: Load bedtools results
import pandas as pd

# ── Load files (only first 6 columns from each BED) ──
five_p_df = pd.read_csv("../../closest_5p.bed", sep="\t", usecols=[0, 1, 2, 3, 4, 5], header=None)
three_p_df = pd.read_csv("../../closest_3p.bed", sep="\t", usecols=[0, 1, 2, 3, 4, 5], header=None)

# ── Rename columns ──
five_p_df.columns = ["enh_chr", "enh_start", "enh_end", "five_chr", "five_start", "five_end"]
three_p_df.columns = ["enh_chr", "enh_start", "enh_end", "three_chr", "three_start", "three_end"]

# ── Create enhancer key for merge ──
five_p_df["key"] = five_p_df["enh_chr"].astype(str) + ":" + five_p_df["enh_start"].astype(str) + "-" + five_p_df["enh_end"].astype(str)
three_p_df["key"] = three_p_df["enh_chr"].astype(str) + ":" + three_p_df["enh_start"].astype(str) + "-" + three_p_df["enh_end"].astype(str)
# ── Deduplicate enhancer entries by key ──
five_p_df = five_p_df.drop_duplicates(subset="key", keep="first")
three_p_df = three_p_df.drop_duplicates(subset="key", keep="first")

# ── Merge based on enhancer coordinate key ──
merged_df = pd.merge(five_p_df, three_p_df, on="key", suffixes=("_5p", "_3p"))

# ── Select final ordered columns ──
temp_df = merged_df[[
    "enh_chr_5p", "enh_start_5p", "enh_end_5p",
    "five_chr", "five_start", "five_end",
    "three_chr", "three_start", "three_end"
]]


# ── Convert coordinates to integers (remove .0) ──
for col in temp_df.columns[1:]:
    temp_df[col] = pd.to_numeric(temp_df[col], errors="coerce").astype("Int64")

# ── Save output ──
temp_df.to_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/old_dont_use_flanking_regions_find_orthologs/redo_correct/enhancer_merged_by_key_cleaned.bed", sep="\t", index=False, header=False)

# ── Pick out columns 1, 5, 9, 1, 2, 3 (1-based → 0,4,8,0,1,2 in 0-based) ──
# If you prefer position-based:
selected_df = temp_df.iloc[:, [0, 4, 8, 0, 1, 2]]

# …or, if you like explicit names (same result):
# selected_df = final_df[
#     ["enh_chr_5p", "five_start", "three_end",
#      "enh_chr_5p", "enh_start_5p", "enh_end_5p"]
# ]

# ── Save without headers or index ──
out_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/old_dont_use_flanking_regions_find_orthologs/redo_correct/flankingcol123_enhancer_col456.bed"
selected_df.to_csv(out_path, sep="\t", header=False, index=False)
