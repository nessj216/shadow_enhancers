import pandas as pd

# Load your TE list
shadow_df = pd.read_csv("/Users/jillianness/Downloads/test_FIXED_FINAL_TE_shadow.bed", sep="\t", header=None)
shadow_df = shadow_df[[4, 5, 6]]  # Chromosome, Start, End
shadow_df.columns = ["Chr", "Start", "End"]
shadow_df["Chr"] = shadow_df["Chr"].str.replace("chr", "", regex=False)

# Load DGRP reference TE list
dgrp_df = pd.read_csv("/Users/jillianness/Downloads/test_DGRP_TE.txt", sep="\t")
#"/Users/jillianness/Desktop/SEbirthanalysis1024/TE_cooption/DGRP_line_analysis/Old_DGRP_ref_only_updated_annotated.txt
dgrp_df.columns = dgrp_df.columns.str.strip()  # Clean column names
dgrp_df = dgrp_df[["Ch", "Pos", "Stop"]]
dgrp_df.columns = ["Chr", "Start", "End"]
dgrp_df["Chr"] = dgrp_df["Chr"].str.strip()

# Function to check for a match within ±200 bp
def is_match(te_chr, te_start, te_end, window=200):
    matches = dgrp_df[
        (dgrp_df["Chr"] == te_chr) &
        (dgrp_df["Start"] >= te_start - window) &
        (dgrp_df["End"] <= te_end + window)
    ]
    return not matches.empty

# Apply to your TE list
shadow_df["Present_in_DGRP"] = shadow_df.apply(
    lambda row: is_match(row["Chr"], row["Start"], row["End"], window=200),
    axis=1
)

# Save or view results
shadow_df.to_csv("TE_matches_in_DGRP.csv", index=False)
print(shadow_df)
