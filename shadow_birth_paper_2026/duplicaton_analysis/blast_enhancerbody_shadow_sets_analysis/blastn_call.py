import os
import shutil
from tqdm import tqdm
from os import path
from pathlib import Path
'''
this iteraviley looks for files of sequences; if in same folder, they are iteraviely blasted using blastn, pairwise.
 the output is a Comparisons folder in each subdir with the blast hits 
the next python script "csv creator" takes the directory and combines all of the blast hits in 
each of the Comparisons output folders
'''
# legacy external path examples removed
ROOT = Path(__file__).resolve().parent
pathTillGene = str(ROOT / "input" / "random_control_test")
list = os.listdir(pathTillGene)
onlyFolder = []
for x in list:
    new=os.path.join(pathTillGene,x)
    if path.isdir(new):
        onlyFolder.append(new)
print (list)
print(onlyFolder)
for current_GENE in tqdm(sorted(onlyFolder)):
    hold = os.listdir(current_GENE)
    sequences = []
    for x in hold:
       # print("DIR", +current_GENE+'/'+x)
        if path.isfile(os.path.join(pathTillGene,current_GENE,x)) and x!='.DS_Store':
            sequences.append(x)
    #print(current_GENE, sequences)
    shutil.rmtree(os.path.join(pathTillGene,current_GENE, 'Comparisons_1'),ignore_errors = True)
    os.mkdir(os.path.join(pathTillGene,current_GENE, 'Comparisons_1'))
    for i in tqdm(range(len(sequences))):
        for j in range(i+1,len(sequences)):
            text1 = sequences[i]
            text2 = sequences[j]
            #print(text1)
            #print(text2)
            file1 = os.path.join(current_GENE,text1)
            file2 = os.path.join(current_GENE,text2)
            #print(file1)
            #print(file2)
            text1 = text1[:-4] # find function to replace extension instead of using removal of 4 characters
            comparison = os.popen("blastn -evalue '.0001' -word_size '10' -gapopen '5' -gapextend '2' -reward '2' -penalty '-3' -dust 'yes' -query "+file1+' -subject '+file2 +' -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue sseq"') ##-outfmt "10  evalue bitscore score pident"
            #comparison = os.popen(" blastn -evalue '.0005' -word_size '7' -gapopen '8' -gapextend '6' -reward '5' -penalty '-4' -dust 'yes' -query "+file1+' -subject '+file2 +' -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue"')

            output = comparison.read()
            #print(output)
            f = open(os.path.join(pathTillGene,current_GENE, 'Comparisons_1', text1+'_'+text2),'w')
            f.write(output)
            f.close()
