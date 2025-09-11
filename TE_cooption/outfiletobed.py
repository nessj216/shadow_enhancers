input_file = "/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/mm10.fa.out"  # replace with your filename
output_file = "/Users/jillianness/Desktop/mouse_analysis_031925/TE_analysis/mm10.fa.bed"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        # Skip header and empty lines
        if line.startswith("   SW") or line.strip() == "" or line.startswith("score"):
            continue

        parts = line.strip().split()

        # Make sure the line has enough columns
        if len(parts) >= 11:
            chrom = parts[4]
            start = parts[5]
            end = parts[6]
            repeat_name = parts[9]
            repeat_class = parts[10]
            outfile.write(f"{chrom}\t{start}\t{end}\t{repeat_name}\t{repeat_class}\n")
