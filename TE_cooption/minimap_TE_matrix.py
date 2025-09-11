import pandas as pd
from te_length_matrix import build_length_matrix  # save te_length_matrix.py somewhere on your PYTHONPATH

df = build_length_matrix(
    paf_dir="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/minimap/splice_tool/minimap_splice_output__sec",
    glob_pat="*.paf",
    flank_len=500,
    min_mapq=20,
    min_flank_cov=0.85,
    out_tsv="te_length_matrix.tsv",   # or None if you don't want to save
)

# inspect / use
print(df.shape)
df.head()
