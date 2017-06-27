#!/usr/bin/env python

import os,sys
import subprocess
from collections import defaultdict
import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2

from gelSeqLib import align_func, io_func

parser = argparse.ArgumentParser(add_help=True,description="find clones of cells with the same VDJ sequence",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('table', metavar="<table>",
                         help='table in csv format')
parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                         help='directory for output')

args = parser.parse_args(sys.argv[1:])
table = args.table
output_dir = args.output_dir

df = pd.read_csv(table)
res_df = pd.DataFrame()
fasta_file = output_dir + "/cdr3.fasta"
with open(fasta_file,'w+') as fas:
    for i,row in df.iterrows():
        if type(row["CDR3 first"]) is not str :
            continue
        fas.write(">"+ "_".join(row["Patient"].split(" ")) + ":" + row["well_id"] + ":" + row["V first"] + "\n")
        fas.write(row["CDR3 first"]+"\n")

alignment_file = align_func.clustalw_align(fasta_file,sys.stdout)

groups = df.groupby("Patient")
'''
for patient, patient_group in groups:
    cdr3_groups = patient_group.groupby(["CDR3 first","V first"])
    d = pd.concat([,cdr3_group])
    res_df = res_df.append([{"Patient":patient, "clone": clone,}])
'''
