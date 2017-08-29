#!/usr/bin/env python
import os
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd


def hamming_distance(str1,str2):
    cost = 0
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            cost += 1
    return cost


def find_same(read,string):
    cost = 100
    pos = -1
    for i in range(4,10):
        h = hamming_distance(read[i:i+len(string)],string)
        if h < cost:
            cost = h
            pos = i
    if cost < 4:
        return pos
    else:
        return -1


dir = "/home/labs/amit/diklag/BaseCalls/"
cell_barcodes = "/home/labs/amit/diklag/BaseCalls/5_CBC PCR1-F_well_IDs_V1.csv"
B_primer = "AGCGACCTCGGGTGGGAAC"
A_primer = "GTACACGGCAGGGTCAGG"
samples=pd.read_csv("/home/labs/amit/diklag/BaseCalls/samples.csv")

list_dir=os.listdir(dir)
list_dir=[a for a in list_dir if "fastq.gz" in a]

index1 = "_L001_I1_001.fastq.gz"
index2 = "_L001_I2_001.fastq.gz"
read1 = "_L001_R1_001.fastq.gz"

global_df = pd.DataFrame(columns=["file","sample","loci","index1","index2","read1"])

global_df_A = pd.DataFrame(columns=["sample","index1","index2","read1"])
global_df_B = pd.DataFrame(columns=["sample","index1","index2","read1"])
samples_indices = defaultdict(dict)

samples = samples.assign(loci = 'A' if samples['Sample_ID'].str.find('tcrA') != -1 else 'B' if samples['Sample_ID'].str.find('tcrB') != -1  else 'NA')
#d=samples[samples['Sample_ID'].str.find('tcrA') != -1]
#samples_indices["tcrA"] = dict([(row["index2"],row["Sample_ID"]) for i,row in d.iterrows()])

for sample in samples['Sample_ID'].values:
    list_sample = [a for a in list_dir if sample+"_" in a]
    counter=0
    for path in list_sample:
        if index1 in path:
            counter+=1
            i1=path
        elif index2 in path:
            counter += 1
            i2=path
        elif read1 in path:
            counter += 1
            r1=path
    if counter==3:
        I1_list=subprocess.getoutput("""gunzip -c %s| awk 'NR%s==2'""" % (os.path.join(dir,i1),"%4")).split("\n")
        I2_list = subprocess.getoutput("""gunzip -c %s| awk 'NR%s==2'""" % (os.path.join(dir, i2), "%4")).split("\n")
        R1_list = subprocess.getoutput("""gunzip -c %s| awk 'NR%s==2'""" % (os.path.join(dir, r1), "%4")).split("\n")
        if samples[samples["Sample_id"] == sample]["loci"].values[0] != "":
            df = pd.DataFrame(data={"file":sample,"sample": sample, "loci":samples[samples["Sample_id"] == sample]["loci"].values[0],"index1": I1_list, "index2": I2_list, "read1": R1_list},
                              columns=["file","sample","loci", "index1", "index2", "read1"])
            global_df.append(df)
        else:
            df = pd.DataFrame(data={"file": "Undetermined", "index1": I1_list, "index2": I2_list, "read1": R1_list},
                              columns=["file", "index1", "index2", "read1"])
            for other_sample in samples['Sample_ID'].values:
                if other_sample == sample:
                    continue
                d = df[]


new_colA = global_df_A['read1'].apply(find_same,string=A_primer)
new_colB = global_df_B['read1'].apply(find_same,string=B_primer)

global_df_A = global_df_A.assign(primer_position = new_colA)
global_df_B = global_df_B.assign(primer_position = new_colB)

new_colA = [row['read1'][0:row['primer_position']] if row['primer_position'] != -1 else 'NA' for i,row in global_df_A.iterrows()]
new_colB = [row['read1'][0:row['primer_position']] if row['primer_position'] != -1 else 'NA' for i,row in global_df_B.iterrows()]

global_df_A = global_df_A.assign(N = new_colA)
global_df_B = global_df_B.assign(N= new_colB)

new_colA = [row['read1'][row['primer_position']+len(A_primer):] if row['primer_position'] != -1 else 'NA' for i,row in global_df_A.iterrows()]
new_colB = [row['read1'][row['primer_position']+len(B_primer):] if row['primer_position'] != -1 else 'NA' for i,row in global_df_B.iterrows()]

global_df_A = global_df_A.assign(read = new_colA)
global_df_B = global_df_B.assign(read= new_colB)

global_df_A.to_csv(dir+"global_A.csv")
global_df_B.to_csv(dir+"global_B.csv")

