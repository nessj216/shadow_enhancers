import os
import pandas as pd
import pybedtools

# Define transcription factor (TF) list
TF_list = ['cad', 'hb', 'bcd', 'kni', 'tll', 'Stat92E', 'gt', 'zelda', 'Kr']

# Define enhancer coordinates
start_ = 21112355
end_ = 21113940
chr_ = 'chr2R'

# Specify the directory containing BED files
directory = "/path/to/your/directory"  # Replace this with your directory path

# Initialize output files
output_file = os.path.join(directory, "krproxcat_filesall.txt")
subfile = os.path.join(directory, "krproxlmotifs.txt")

# Function to process each BED file
def process_bed_file(bed_file, tf_name):
    # Load BED file into pybedtools object
    bed = pybedtools.BedTool(bed_file)

    # Define enhancer region
    enhancer_region = pybedtools.BedTool([[chr_, start_, end_]])

    # Intersect BED file with enhancer region
    intersection = bed.intersect(enhancer_region, wa=True)

    # Process the intersection results
    results = []
    for interval in intersection:
        binding_site = interval.start
        abs_val = abs(int(interval.start) - int(interval.end))
        strand = interval.strand
        tf = tf_name
        if strand == "+":
            val = int(interval.start) - start_
        elif strand == "-":
            val = int(interval.start) - start_
        else:
            val = ""
        results.append([binding_site, abs_val, interval.score, val, tf, strand])

    # Save results to subfile
    df = pd.DataFrame(results, columns=['binding_site', 'abs_val', 'score', 'val', 'tf', 'strand'])
    df.to_csv(subfile, mode='a', sep='\t', index=False, header=False)

    print(f"Processed file: {bed_file} and saved output to {subfile}")

# Loop through all files in the specified directory matching the pattern
for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith("bound.bed"):
            bed_file = os.path.join(root, file)
            # Find matching transcription factor
            tf_name = next((tf for tf in TF_list if tf in file), None)
            if tf_name:
                process_bed_file(bed_file, tf_name)

# Combine results from subfile into the final output file
with open(output_file, 'w') as outfile:
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == os.path.basename(subfile):
                with open(os.path.join(root, file), 'r') as infile:
                    outfile.write(infile.read())

print(f"All results have been combined into {output_file}")
