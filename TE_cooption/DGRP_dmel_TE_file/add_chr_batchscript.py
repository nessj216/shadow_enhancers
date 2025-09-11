#!/usr/bin/env python3
import os
from pathlib import Path

# ==== EDIT THESE ====
input_dir = "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TEsREPET-2/ReferenceCoordinates"
output_dir = "/Users/jillianness/Desktop/SEbirthanalysis1024/TEs_otherspecies/REPET_Rechetal/TEsREPET-2/add_chr_referenceCoordinates"
# ====================

os.makedirs(output_dir, exist_ok=True)

for bed_file in Path(input_dir).glob("*.bed"):
    out_path = Path(output_dir) / bed_file.name
    with open(bed_file, "r") as infile, open(out_path, "w") as outfile:
        for line in infile:
            if line.strip() == "" or line.startswith(("track", "browser", "#")):
                outfile.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts[0].startswith("chr"):
                parts[0] = "chr" + parts[0]
            outfile.write("\t".join(parts) + "\n")
    print(f"Processed: {bed_file.name} -> {out_path}")

print("Done.")
