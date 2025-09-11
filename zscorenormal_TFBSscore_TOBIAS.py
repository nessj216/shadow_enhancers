import numpy
import os
#import matplotlib.pyplot as plt
import statistics
#import xlwt

directory = "/Users/jillianness/Desktop/TOBIAS091922/d7"
path_enh=0
#for file1 in os.listdir(directory):
    #filename: str = os.fsdecode(file)
   # print(file1)
enhancer = {}
TFscores={}
tempor=[]
tempor1=[]
TFstats={}
norm_scores=[]
norm_dict={}
avg_norm_dict={}
sum_norm_scores={}

TFs=['gt','cad','bcd','kni', 'zelda', 'Kr', 'Stat','hb','tll']
for path, subdir, files in os.walk(directory):

    for file1 in files:
        tempor=[]
        tempor1=[]
        filename = str(os.fsdecode(file1)) #just filename
        filename1 = os.path.join(path, filename) #filename and path

        if os.path.getsize(filename1) == 0:
            continue

        if 'eve2' in filename:

            with open(filename1, "r") as t:
                for line in t:
                    tempor1.append(float(line.rstrip('\n'))) #clear temp list every outer loop iteration
                enhancer[path]=tempor1 #key:path to each file, value is contents of file
        #print(enhancer)
        if 'all.bed.txt'in filename:
            with open(filename1, "r") as t:
                for line in t:
                    tempor.append(float(line.rstrip('\n'))) #clear temp list every outer loop iteration
                TFscores[path]=tempor #key:path to each file, value is contents of file
'''print(TFscores)
for key,value in TFscores.items():
    print("key: %s" % (key))
    print("Val: %s" % (len(value)))'''


for key,value in TFscores.items():

    avg=statistics.mean(value)
    sd=statistics.stdev(value)
    TFstats[key]=(avg, sd)
    #print(key,avg, sd)
    #print(key, '%s' % float('%1g' % avg), '%s' % float('%1g' % sd))
#print(TFstats)

for key in enhancer:
    norm_scores = []
    for i in enhancer[key]:
        #print(TFstats[key][0],TFstats[key][1])

        norm_scores.append((i-TFstats[key][0])/TFstats[key][1])

    norm_dict[key]=norm_scores # norm_dict is the dictionary containing all domain normalized TFBS footprint scores

print(norm_dict)
row=0
for key, value in norm_dict.items():
    sum_norm_scores[key]=sum(value)
print(sum_norm_scores)
for TF in TFs:
    avg_norm_dict[TF] = []
for key,value in sum_norm_scores.items():
    for TF in TFs:

        if TF in key:

            #print(TF,key)
            avg_norm_dict[TF].append(value) #this will give be sum of TF scores per TF in the enhancer


out = [{k: sum(v)/len(v)} for k,v in avg_norm_dict.items() if v!=[]]
print(out)