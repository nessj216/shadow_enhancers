#UCi undergrads modified code to covert the xlsx file containing SE coordinates into seperate bedfiles, then call getFasta and put
# the enhancer sequnce in the correct gene folder

import os
import pandas as pd

mainpath = "/Users/jillianness/Desktop/sorting_cannavo_data"
df = pd.read_excel('/Users/jillianness/Desktop/SEbirthanalysis1024/011925_all_shadowsets_DM6.bed', engine='openpyxl')
df = df.dropna()

Enhancer_ID = df["CRM Name"].values.tolist()
Chromosome = df["Chromosome"].values.tolist()
Start = df["Start"].values.tolist()
End = df["End"].values.tolist()
Gene = df["Gene Name"].values.tolist()

# Create the directory if it doesn't exist
bed_folder = mainpath + "/FliesBED/"
if not os.path.exists(bed_folder):
    os.makedirs(bed_folder)

for i in range(len(Enhancer_ID)):
    start_int = int(Start[i])
    end_int = int(End[i])
    temp = str(Chromosome[i]) + "\t" + str(start_int) + "\t" + str(end_int) + "\t" + str(Gene[i])
    filepath = mainpath + "/FliesBED/"
    Enhancer_ID[i] = str(Enhancer_ID[i])
    filename = filepath + Enhancer_ID[i] + ".txt"
    file = open(filename, "w")
    file.write(temp)
    file.close()



for gene in Gene:
    path = mainpath + "/FliesOutput/" + str(gene) + "/"
    if not os.path.exists(path):
        os.makedirs(path)

for i in range(len(Enhancer_ID)):
    bed = mainpath + "/FliesBED/" + str(Enhancer_ID[i]) + ".txt"
    out = mainpath + "/FliesOutput/" + Gene[i] + "/" + str(Enhancer_ID[i]) + "output.txt"
    cmd = 'bedtools getfasta -fi ' + '/Users/jillianness/Desktop/TOBIAS091922/general/dm3.fasta' + " " + '-bed ' + bed + " " + '-fo ' + out
    os.system(cmd)
