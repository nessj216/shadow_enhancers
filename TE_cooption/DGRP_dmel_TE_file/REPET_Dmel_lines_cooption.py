#compare my ISO1 TE co-opted hits against the REPET ~40 d mel lines


import os
import subprocess

# Paths to files and directories (example placeholders, update with your actual paths)
my_bed_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_library_DM6/dm6_Final_combinedTEs/final_masker_outputRM_ONECODEfile.bed"
bed_dir = "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TEsREPET-2/add_chr_referenceCoordinates"  # Folder containing multiple BED files
output_dir = "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/masker_repet_overlap"

import os
import subprocess
import pandas as pd

# CONFIGURATION
a_file = "/Users/jillianness/Desktop/SEbirthanalysis1024/TE_library_DM6/dm6_Final_combinedTEs/final_masker_outputRM_ONECODEfile.bed"      # your -a file
b_dir = "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TEsREPET-2/add_chr_referenceCoordinates" # directory containing -b files
window_size = 1000          # example: change to desired window size
bedtools_path = "bedtools"  # path to bedtools executable if not in PATH

# Count lines in a_file
with open(a_file) as f:
    a_lines = sum(1 for _ in f)

results = []

# Loop through files in directory
for b_file in sorted(os.listdir(b_dir)):
    b_path = os.path.join(b_dir, b_file)

    if not os.path.isfile(b_path) or not b_file.endswith(".bed"):
        continue

    cmd = [
        bedtools_path, "window",
        "-a", a_file,
        "-b", b_path,
        "-w", str(window_size)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error with {b_file}: {result.stderr}")
        out_lines = None
    else:
        out_lines = len(result.stdout.strip().split("\n")) if result.stdout.strip() else 0

    results.append({
        "B File": b_file,
        "Lines in A": a_lines,
        "Lines in Output": out_lines,
        "Ratio": (out_lines / a_lines) if (out_lines is not None and a_lines > 0) else None
    })

# Save results to Excel
df = pd.DataFrame(results)
df.to_excel("bedtools_results.xlsx", index=False)

print("Results saved to bedtools_results.xlsx")

"""# Paths
MY_BED="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TE_shadow_for_bedtools.bed"
BED_DIR="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TEsREPET-2/add_chr_referenceCoordinates"
OUT_DIR="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TEsREPET-output"


# Loop over each bed file in the folder
for FILE in "$BED_DIR"/*.bed; do
    BASENAME=$(basename "$FILE" .bed)

    
    # Run bedtools window
    bedtools window -a "$MY_BED" -b "$FILE" -w 1000 > "$OUT_DIR/${BASENAME}_window1000.bed"
done
"""

'''
#!/usr/bin/env bash
set -u  # (no -e, we want to continue on per-file errors)

MY_BED="/Users/jillianness/Desktop/SEbirthanalysis1024/TE_library_DM6/dm6_Final_combinedTEs/final_masker_outputRM_ONECODEfile.bed"
BED_DIR="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TEsREPET-2/add_chr_referenceCoordinates"
OUT_DIR="/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/REPEATM_repet_overlap"
TMP_DIR="${OUT_DIR}/_tmp_clean"

mkdir -p "$OUT_DIR" "$TMP_DIR"

shopt -s nullglob   # if no *.bed, loop won’t run with a literal pattern

for FILE in "$BED_DIR"/*.bed; do
  BASENAME=$(basename "$FILE" .bed)
  CLEAN="${TMP_DIR}/${BASENAME}.clean.bed"
  OUT="${OUT_DIR}/${BASENAME}_window1000.bed"

  echo ">>> Processing: $BASENAME"

  # Minimal normalization:
  awk '
    BEGIN{OFS="\t"}
    { gsub("\r",""); }                 # strip CR
    /^$/ || /^track/ || /^browser/ || /^#/ {next}
    { gsub(/[ ]+/, "\t"); }            # collapse spaces to tabs
    { gsub(",", "", $2); gsub(",", "", $3); }  # remove commas in start/end
    $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ { print; }
  ' "$FILE" > "$CLEAN"

  # quick check: if clean file is empty, skip
  if ! [ -s "$CLEAN" ]; then
    echo "    [SKIP] No valid rows after cleaning"
    continue
  fi

  # Run bedtools (continue on failure, but report)
  if bedtools window -a "$MY_BED" -b "$CLEAN" -w 1000 > "$OUT" 2> "${OUT}.stderr"; then
    echo "    [OK] ${OUT}"
  else
    echo "    [FAIL] bedtools error. See ${OUT}.stderr"
  fi
done

echo "Done. Outputs: $OUT_DIR  | Cleaned: $TMP_DIR"
'''