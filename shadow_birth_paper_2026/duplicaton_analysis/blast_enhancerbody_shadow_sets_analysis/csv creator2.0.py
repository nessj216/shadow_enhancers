import os
import shutil
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
pathTillGene = str(ROOT / "input" / "random_control_test")
# legacy external path example removed
list = os.listdir(pathTillGene)
onlyFolder = []
for x in list:
    new = os.path.join(pathTillGene, x)
    if os.path.isdir(new):
        onlyFolder.append(new + '/Comparisons_1')

shutil.rmtree(pathTillGene + '/Comparisons_1', ignore_errors=True)

with open(pathTillGene + '/CollatedComparisons_null.csv', mode='w') as csv_file:
    csvFileWriter = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csvFileWriter.writerow(
        ['Gene Name', 'Comparisons', 'qseqid', 'sseqid', 'pident', 'length','qstart', 'qend', 'sstart','send',  'evalue', 'sseq'])

    for compLoc in onlyFolder:
        onlyComp = os.listdir(compLoc)

        for compFile in onlyComp:
            with open(os.path.join(compLoc, compFile), "r") as f:
                data = f.readlines()
                if len(data) != 0:
                    geneName = compLoc.split("/")[-2]
                    for line in data:
                        line_data = line.strip().split('\t')
                        line_data.insert(0, compFile)
                        line_data.insert(0, geneName)
                        csvFileWriter.writerow(line_data)
