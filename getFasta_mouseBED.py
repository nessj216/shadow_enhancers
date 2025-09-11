import os
import pandas as pd
#might need to add an additional index
# Define main paths
#mainpath = "/Users/jillianness/Desktop/mouse_analysis_031925/Fasta_files"
#input_bed = "/Users/jillianness/Desktop/mouse_analysis_031925/enhancers_formatted.bed"
#fasta_file = "/Users/jillianness/Desktop/mouse_analysis_031925/GRCm38.p6.genome.fa"

import os
import pandas as pd

mainpath = "/Users/jillianness/Desktop/mouse_analysis_031925/Mouse_enhancers/Fasta_files/Singles_fasta_files"
df = pd.read_csv('/Users/jillianness/Desktop/mouse_analysis_031925/Mouse_enhancers/singlesmouse_file.csv')
df = df.dropna()

Chromosome = df.iloc[:, 0].values.tolist()
Start = df.iloc[:, 1].values.tolist()
End = df.iloc[:, 2].values.tolist()
Enhancer_Index = df.iloc[:, 4].values.tolist()
Gene = df.iloc[:, 3].values.tolist()

# Create the directory if it doesn't exist
bed_folder = mainpath + "/FliesBED/"
if not os.path.exists(bed_folder):
    os.makedirs(bed_folder)

for i in range(len(Enhancer_Index)):
    start_int = int(Start[i])
    end_int = int(End[i])
    temp = f"{Chromosome[i]}\t{start_int}\t{end_int}\t{Enhancer_Index[i]}\t{Gene[i]}"
    filepath = bed_folder
    filename = filepath + str(Enhancer_Index[i]) + ".txt"

    with open(filename, "w") as file:
        file.write(temp)
        file.close()

# Create output folders per gene
for gene in Gene:
    path = mainpath + "/FliesOutput/" + str(gene) + "/"
    if not os.path.exists(path):
        os.makedirs(path)

# Generate FASTA sequences for each enhancer
for i in range(len(Enhancer_Index)):
    bed = bed_folder + str(Enhancer_Index[i]) + ".txt"
    out = mainpath + "/FliesOutput/" + Gene[i] + "/" + str(Enhancer_Index[i]) + "_output.txt"
    cmd = f'bedtools getfasta -fi /Users/jillianness/Desktop/mouse_analysis_031925/GRCm38.p6.genome.fa -bed {bed} -fo {out}'
    os.system(cmd)