#!/usr/bin/env python


import numpy as np
import pandas as pd
import matplotlib
import networkx
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from itertools import combinations
import pickle
import numpy
from Bio import pairwise2
from Bio import AlignIO
# from gelSeqLib import align_func, io_func


'''


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

'''

aln = AlignIO.read('/home/labs/amit/diklag/output/cdr3.aln','clustal')

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
with open('/home/labs/amit/diklag/output/dm.pkl','wb') as f:
    pickle.dump(dm,f,protocol=0)

l = list(combinations(range(len(dm.names)),2))
distmat = np.repeat(np.inf, len(l))

for index in range(len(l)):
    distmat[index] = dm.matrix[l[index][1]][l[index][0]]


with open('/home/labs/amit/diklag/output/distmap.pkl','wb') as f:
    pickle.dump(distmat,f,protocol=0)
