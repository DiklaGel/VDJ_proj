#!/usr/bin/env python
import pickle
import sys
import numpy as np
from itertools import combinations
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from matplotlib import pyplot as plt
from matplotlib import axes
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
import pandas as pd
import subprocess
sys.path.insert(0, '/home/labs/amit/diklag/PycharmProjects/VDJ_proj/')
from gelSeqLib import align_func, io_func
from Bio import AlignIO
from scipy.cluster.hierarchy import fcluster


def cluster_by_cdr3(results_table, output_dir):
    df1 = pd.read_csv(results_table)
    res_df = pd.DataFrame()
    fasta_file = output_dir + "/cdr3.fasta"
    with open(fasta_file, 'w+') as fas:
        for i, row in df1.iterrows():
            if type(row["CDR3 first"]) is not str:
                continue
            fas.write(">" + "_".join(row["Patient"].split(" ")) + ":" + row["well_id"] + ":" + row["V first"] + "\n")
            fas.write(row["CDR3 first"] + "\n")

    alignment_file = align_func.clustalw_align(fasta_file, sys.stdout)

    aln = AlignIO.read(alignment_file, 'clustal')

    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    with open(output_dir+'/dm.pkl', 'wb') as f:
        pickle.dump(dm, f, protocol=0)

    l = list(combinations(range(len(dm.names)), 2))
    distmat = np.repeat(np.inf, len(l))

    for index in range(len(l)):
        distmat[index] = dm.matrix[l[index][1]][l[index][0]]

    with open(output_dir+'/distmat.pkl', 'wb') as f:
        pickle.dump(distmat, f, protocol=0)

    Z = linkage(distmat, method='average')

    max_d = 0.05
    clusters = fcluster(Z, max_d, criterion='distance')

    patient_col = [x.split("_W")[0] for x in dm.names]
    well_col = ['W' + x.split("_W")[1].split("_")[0] for x in dm.names]
    df2 = pd.DataFrame(data={"cluster": clusters, "Patient": patient_col, "well_id": well_col})

    table = pd.merge(df2, df1, on=["Patient", "well_id"], how="inner")

    table = table[['cluster', 'Patient', 'Amp Batch', 'well_id', 'cell_name', '#reads', '#umi distribution', "V first",
                   "V first counts",
                   "V second", "V second counts",
                   "D first", "D first counts",
                   "D second", "D second counts",
                   "J first", "J first counts",
                   "J second", "J second counts",
                   "CDR3 first", "CDR3 first counts",
                   "CDR3 second", "CDR3 second counts"]]

    table.to_csv(output_dir+'/full_results.csv')



cluster_by_cdr3("/home/labs/amit/diklag/output/130817/result.csv","/home/labs/amit/diklag/output/130817")






