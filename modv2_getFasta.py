import os
import pandas as pd

mainpath = "/Users/jillianness/Desktop/sorting_cannavo_data"

# Check the available sheet names first
#xls = pd.ExcelFile('/Users/jillianness/Desktop/sorting_cannavo_data/cannavo_etal_files/processed_crm_final_cleaned_reordered.xlsx', engine='openpyxl')
print("Available sheet names:", xls.sheet_names)

# After verifying the sheet names, load the correct sheet
# Replace 'CorrectSheetName' with the actual sheet name you find
df = pd.read_excel(xls, sheet_name=xls.sheet_names[0])  # Load the first sheet if you're unsure
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
    with open(filename, "w") as file:
        file.write(temp)

for gene in Gene:
    path = mainpath + "/FliesOutput/" + str(gene) + "/"
    if not os.path.exists(path):
        os.makedirs(path)

for i in range(len(Enhancer_ID)):
    bed = mainpath + "/FliesBED/" + str(Enhancer_ID[i]) + ".txt"
    out = mainpath + "/FliesOutput/" + Gene[i] + "/" + str(Enhancer_ID[i]) + "output.txt"
    cmd = 'bedtools getfasta -fi ' + '/Users/jillianness/Desktop/TOBIAS091922/general/dm3.fasta' + " -bed " + bed + " -fo " + out
    os.system(cmd)
