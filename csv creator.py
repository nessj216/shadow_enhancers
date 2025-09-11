import os
import shutil
##from tqdm import tqdm
from os import path
import csv

pathTillGene = "/Users/jillianness/Desktop/all_mouse_seq" #put the main dir where all sequences are held
list = os.listdir(pathTillGene)
print(list)
onlyFolder = []
for x in list:
    new = os.path.join(pathTillGene, x)
    if path.isdir(new):
        onlyFolder.append(new+'/Comparisons') #put the folder where the blast output is
print (onlyFolder)
shutil.rmtree(pathTillGene+'/CollatedComparisons.csv',ignore_errors = True)
with open(pathTillGene+'/CollatedComparisons.csv', mode='w') as csv_file:
    csvFileWriter = csv.writer(csv_file, delimiter=',',quotechar='"', quoting = csv.QUOTE_MINIMAL)
    csvFileWriter.writerow(['Gene Name','Comparisons', 'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'evalue'])
    currentGene = ''
    for compLoc in onlyFolder:
        onlyComp = os.listdir(compLoc)

        for compFile in onlyComp:
            f = open(compLoc +'/'+compFile,"r")
            data = f.read()
            if len(data) != 0:
            #print(data)
                data = data.strip()
                geneName = compLoc.split("/")[-2]
                if geneName != currentGene and currentGene != '':
                    csvFileWriter.writerow(['-','-','-','-','-','-'])
                currentGene = geneName
                #print(geneName)
                splitData = data.split(',')
                splitData.insert(0,compFile)
                splitData.insert(0,geneName)
                #print(splitData)
                csvFileWriter.writerow(splitData)