import pandas as pd

import pandas as pd

# ——— Load File 2 ———
file2 = pd.read_csv(
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/shadow/take2/Dsim_overlap_results.bed",
    sep='\t', header=None
)

# ——— NEW: keep only rows where col 4 (idx 3) equals col 11 (idx 10) ———
matching = file2[file2[3] == file2[10]]

# (a) write out the filtered rows as a BED
matching.to_csv(
    "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/shadow/take2/Filteredmatch4and11_Dsim_overlap_results.bed",
    sep='\t', header=False, index=False
)

# Load File 1 (no header)
#file1 = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/take2/try_again/orig_ID_enhancer_TE_ids_corrected.csv", header=None)
file1 = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/shadow/take2/enhancers_with_flanks_and_extras_withIDs.bed",sep='\t',  header=None)

# Load File 2 (BED format with tab delimiter, no header)
file2 = pd.read_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/shadow/take2/Filteredmatch4and11_Dsim_overlap_results.bed", sep='\t', header=None)

# ---- Part 1: Count unique TE_IDs per enhancer_ID from File 1 ----
# enhancer_ID = column index 12, TE_ID = column index 13 ; 0 start
counts_file1 = file1.groupby(file1.columns[12])[file1.columns[13]].nunique().reset_index()
counts_file1.columns = ['enhancer_ID', 'num_TE_IDs_file1']

# ---- Part 2: Count TE_IDs per enhancer_ID from File 2,
# Only count one TE_ID per unique (enhancer_ID, col11) group ----
# old file columns: TE_ID = col4 (index 4), enhancer_ID = col5 (index 5), grouping col = col11 (index 10)
# TE_ID = col5 (index 5), enhancer_ID = col6 (index 6), grouping col = col9 (index 9)
#df2_grouped = file2.groupby([5, 10])[4].first().reset_index()
#counts_file2 = df2_grouped.groupby(5)[4].nunique().reset_index()
df2_grouped = file2.groupby([6, 9])[5].first().reset_index()
counts_file2 = df2_grouped.groupby(6)[5].nunique().reset_index()
counts_file2.columns = ['enhancer_ID', 'num_TE_IDs_file2_corrected']

# ---- Merge and Output ----
final_counts = counts_file1.merge(counts_file2, on='enhancer_ID', how='outer')

# Save to CSV
final_counts.to_csv("/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/final_flankingregions_find_orthologs/shadow/take2/filter_TE_ID_counts_by_enhancer.csv", index=False)
