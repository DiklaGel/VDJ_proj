
import sys, os,csv, operator, tempfile, subprocess, numpy, Bio
import pandas as pd
import numpy as np
from itertools import groupby, count
from collections import Counter
from Bio.Blast import NCBIXML
from Bio import pairwise2
from collections import defaultdict
import warnings
import copy
import io

# running multiple alignment via clustalOmega
def clustalo_align(fasta,log_fd):
    in_file = fasta
    out_file = ".".join(fasta.split(".")[0:len(fasta.split(".")) - 1]) + ".aln"
    log_fd.write(subprocess.getoutput("/apps/RH7U2/gnu/clustal-omega/1.2.4/bin/clustalo -i " + in_file + " -o " + out_file + " --maxseqlen=95 -v"))
    return out_file

# running multiple alignment via clustaW2
# This program is slower and older than clustalo, but it contains the option to customize the penalties for gaps
def clustalw_align(fasta,log_fd):
    in_file = fasta
    out_file = ".".join(fasta.split(".")[0:len(fasta.split(".")) - 1]) + ".aln"
    log_fd.write(subprocess.getoutput("/apps/RH7U2/gnu/clustalw/2.1/bin/clustalw2 -infile=" + in_file + " -GAPOPEN=15 -GAPEXT=15"))
    return out_file

# generate a consensus sequence by major rule
# @format_alignment = clustal (in clustalw case) /fasta (in clustalo case)
def make_consensus(alignment_file, format_alignment):
    from Bio import AlignIO
    from Bio.Align.AlignInfo import SummaryInfo
    format = format_alignment
    alignment = AlignIO.read(alignment_file, format)
    summary = SummaryInfo(alignment)
    mat = summary.pos_specific_score_matrix()
    consensus = ""
    for row in mat:
        key = max(row.items(), key=operator.itemgetter(1))[0]
        consensus = "".join([consensus,key])
    consensus = consensus.strip("-")
    return consensus

# simple hamming distance
def hamming_distance(str1,str2):
    cost = 0
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            cost += 1
    return cost

from math import log
def gap_function_read(x, y):  # x is gap position in seq, y is gap length
    if x < 4 or x > 70: # gaps in the start/end of the sequence are ok
        return 0
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Gap open penalty
        return -2
    return - (2 + y/4.0 + log(y)/2.0)

def gap_function_consensus(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    return -100 # new gaps in the consensus sequence are not welcome
