

#needed to add IDs to inital D mel species file then deletes unncessary columns, then send into liftover, then bedtools intersect ovwvevrlap
#then and only then run the following code ex. of input file: " /Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/take2/try_again/input_ID_enhancer_TE_ids_corrected.txt"

# Reload the original BED file (to ensure all rows are present)
import pandas as pd
import string
df = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/take2/try_again/Dpseudo/Dpseudo_overlap.bed", sep='\t', header=None)

# Step 1: Remove duplicate rows
df = df.drop_duplicates()

# Step 2 (relaxed): Keep rows where both col 4 and col 12 share one of the key TE types
def match_te_category(row):
    te_classes = ['LTR', 'LINE', 'SINE', 'DNA']
    for te in te_classes:
        if pd.notna(row[3]) and pd.notna(row[11]) and te in str(row[3]) and te in str(row[11]):
            return True
    return False

df = df[df.apply(match_te_category, axis=1)]

# Step 3: Assign alphanumeric IDs to sets of columns 8, 9, 10 (indexes 7, 8, 9)
flank_sets = df[[7, 8, 9]].drop_duplicates().reset_index(drop=True)

def generate_alphanumeric_ids(n):
    import string
    letters = string.ascii_uppercase
    ids = []
    i = 1
    while len(ids) < n:
        for letter in letters:
            ids.append(f"{letter}{i}")
            if len(ids) == n:
                break
        i += 1
    return ids

flank_sets['flankSet_ID'] = generate_alphanumeric_ids(len(flank_sets))

# Merge back the flankSet_ID
df = df.merge(flank_sets, on=[7, 8, 9], how='left')

# Step 4: Drop column 5 (index 4)
df = df.drop(columns=[4])

# Save to output TXT file
output_path = "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/take2/try_again/Dpseudo/Dpseudo_overlap_filtered_TEtype_match_with_flankSetID.txt"
df.to_csv(output_path, sep='\t', index=False, header=False)


import pandas as pd

# Load File 1 (no header)
file1 = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/take2/try_again/orig_ID_enhancer_TE_ids_corrected.csv", header=None)

# Load File 2 (BED format with tab delimiter, no header)
file2 = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/take2/try_again/Dpseudo/Dpseudo_overlap_filtered_TEtype_match_with_flankSetID.txt", sep='\t', header=None)

# ---- Part 1: Count unique TE_IDs per enhancer_ID from File 1 ----
# enhancer_ID = column index 11, TE_ID = column index 12
counts_file1 = file1.groupby(file1.columns[11])[file1.columns[12]].nunique().reset_index()
counts_file1.columns = ['enhancer_ID', 'num_TE_IDs_file1']

# ---- Part 2: Count TE_IDs per enhancer_ID from File 2,
# Only count one TE_ID per unique (enhancer_ID, col11) group ----
# TE_ID = col4 (index 4), enhancer_ID = col5 (index 5), grouping col = col11 (index 10)
df2_grouped = file2.groupby([5, 10])[4].first().reset_index()
counts_file2 = df2_grouped.groupby(5)[4].nunique().reset_index()
counts_file2.columns = ['enhancer_ID', 'num_TE_IDs_file2_corrected']

# ---- Merge and Output ----
final_counts = counts_file1.merge(counts_file2, on='enhancer_ID', how='outer')

# Save to CSV
final_counts.to_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/take2/try_again/Dpseudo/Dpseudo_TE_ID_counts_by_enhancer.csv", index=False)