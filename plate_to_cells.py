
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
from . import io, align_func

################################################################
#  --------------------------to do-----------------------------#
# need to think of a better threshold function
# this function determines what is the minimal number for
################################################################

def threshold(groups):
    hist, not_important = np.histogram(groups, density=True, bins=groups[len(groups) - 1])
    sum = 0
    i = -1
    while sum + hist[i+1] < 0.99 and i < 500:
        sum += hist[i+1]
        i += 1
    return i



def reads_to_fasta(cell_barcodes, dir_name, fastq1, fastq2, cell_name):
    io.makeOutputDir(dir_name)
    # finding the lines (list of numbers) in the fastq2 file that start with those 15-mers
    lines = io.find_lines(cell_barcodes, fastq2)
    # extracting the reads (list of strings) in the fastq2 file that are written in the lines we found above
    reads2 = io.get_full_reads(lines, fastq1)
    fastq_path = os.path.join(dir_name, cell_name + ".fastq")
    io.reads_to_fastq(reads2, fastq_path)
    # converting fastq to fasta for further using in multiple alignment
    fasta_path = io.fastq_to_fasta(fastq_path)
    return fasta_path

def cell_consensus(most_significant_barcode, dir_name, fastq1, fastq2, cell_name,log_fd):
    # extracting the most significant 15-mers of the i cell
    cell_barcodes = [
        "".join([most_significant_barcode.iloc[i]["cell_barcode"], most_significant_barcode.iloc[i]["umi_barcode"]]) for
        i in range(0, len(most_significant_barcode))]
    # finding the lines (list of numbers) in the fastq2 file that start with those 15-mers
    lines = io.find_lines(cell_barcodes, fastq2)
    # extracting the reads (list of strings) in the fastq2 file that are written in the lines we found above
    reads2 = io.get_full_reads(lines, fastq1)
    fastq_path = os.path.join(dir_name, cell_name + ".fastq")
    io.reads_to_fastq(reads2, fastq_path)
    # converting fastq to fasta for further using in multiple alignment
    fasta_path = io.fastq_to_fasta(fastq_path)
    alignment_file = align_func.clustalw_align(fasta_path,log_fd)
    # generating a consensus sequence for the reads with the most abundant 15-mer barcode in cell
    consensus = align_func.make_consensus(alignment_file, "clustal")
    print("consensus for cell " + cell_name, flush=True)
    print(consensus,flush=True)
    return consensus

def filter_real_reads(consensus, other_significant_barcode,dir_name, fastq1, fastq2, cell_name,log_fd):
    # extracting the most significant 15-mers of the i cell
    cell_barcodes = [
        "".join([other_significant_barcode.iloc[i]["cell_barcode"], other_significant_barcode.iloc[i]["umi_barcode"]]) for
        i in range(0, len(other_significant_barcode))]
    # finding the lines (list of numbers) in the fastq2 file that start with those 15-mers
    lines = io.find_lines(cell_barcodes, fastq2)
    # extracting the reads (list of strings) in the fastq2 file that are written in the lines we found above
    min = 10000
    max = -10000
    sum = 0
    count = 0
    for line in lines:
        read = io.get_reads([line], fastq1)[0]
        #alignment = pairwise2.align.globalmc(consensus,read,1,0,gap_function_consensus,gap_function_read, gap_char = '_', score_only = True)
        # score = max([alignment[i][2] for i in range(0,len(alignment))])
        score = pairwise2.align.globalmc(consensus,read,1,0,align_func.gap_function_consensus,align_func.gap_function_read, gap_char = '-', score_only = True)
        sum += score
        count += 1
        if score < min:
            min = score
        if score > max:
            max = score
    avg = sum/count
    print(cell_name + " average score against consesnsus = " + str(avg),flush=True)
    print(cell_name + " min score against consesnsus = " + str(min),flush=True)
    print(cell_name + " max score against consesnsus = " + str(max),flush=True)


def filter_abundant_barcodes(fastq2):
    # reading all 15-mers in the current plate (current fastq2 - fastq file of read2), sorting them by frequency
    column = subprocess.getoutput("""gunzip -c %s | awk 'NR%s==2' | cut -b 1-15 | sort | uniq -c | sort -n """ % (fastq2, "%4")).split("\n")
    # splitting the list of "count 15-mer" to 3-dim list [count][cell barcode][umi barcode]
    columns = [(int(column[i].strip().split(" ")[0]), column[i].strip().split(" ")[1][0:7],
                column[i].strip().split(" ")[1][7:15]) for i in range(0, len(column))]
    # creating a data frame from the 3-dim list above
    plate_seqs = pd.DataFrame(columns, columns=["num", "cell_barcode", "umi_barcode" ])
    # creating a np.array for making an histogram of the counts
    x = np.array([int(column[i].strip().split(" ")[0]) for i in range(0, len(column))])
    # find the threshold number
    t = threshold(x)

    # save only the 15-mers that appear more than t (threshold) times
    high_confidence_barcodes = plate_seqs[plate_seqs["num"] > t]
    if len(high_confidence_barcodes) == 0:
        return
    high_confidence_barcodes = high_confidence_barcodes.sort_values(by= "num", ascending=False)

    return high_confidence_barcodes


def split_by_cells(high_confidence_barcodes,wells_cells_file,output_dir,fastq1,fastq2):
    cells_to_path = dict()
    map_cell_to_barcode = pd.read_csv(wells_cells_file, delimiter='\t',
                                      usecols=['well_coordinates', 'Cell_barcode', 'Amp_batch_ID', 'plate_ID'])
    # grouping to cells by barcodes
    checked_cells, cell_barcode_mapping = group_to_cells(high_confidence_barcodes, map_cell_to_barcode)

    # generate fasta file for each cell
    for cell_barcode in checked_cells.keys():
        # extracting the rows of the current cell
        cell_barcodes_rows = pd.DataFrame()
        for barcode in cell_barcode_mapping[cell_barcode]:
            cell_barcodes_rows = pd.concat([cell_barcodes_rows,high_confidence_barcodes[high_confidence_barcodes["cell_barcode"] == barcode]])
        cell_name = map_cell_to_barcode[map_cell_to_barcode['Cell_barcode'] == cell_barcode]['well_coordinates'].tolist()[0]

        cell_dir = os.path.join(output_dir,cell_name)
        io.makeOutputDir(cell_dir)

        cell_barcodes_rows = cell_barcodes_rows.sort_values(by= "num", ascending=False)
        cell_barcodes = ["".join([cell_barcodes_rows.iloc[i]["cell_barcode"], cell_barcodes_rows.iloc[i]["umi_barcode"]]) for i in range(0,len(cell_barcodes_rows))]
        fasta_path = reads_to_fasta(cell_barcodes,cell_dir+"/"+'filtered_reads',fastq1,fastq2,cell_name)
        cells_to_path[cell_name] = fasta_path
    return cells_to_path


def group_to_cells(high_confidence_barcodes, map_cell_to_barcode):
    checked_cells = defaultdict(list)
    cell_barcode_mapping = defaultdict(list)
    for index,row_high in high_confidence_barcodes.iterrows():
        cell_barcode = row_high["cell_barcode"]
        umi_barcode = row_high["umi_barcode"]
        if cell_barcode in checked_cells.keys(): # if we already pass this barcode
            continue
        if len(map_cell_to_barcode['Cell_barcode'].str.contains(cell_barcode)[map_cell_to_barcode['Cell_barcode'].str.contains(cell_barcode)
                                                                             == True]) > 0:
            checked_cells[cell_barcode].append(umi_barcode)
            cell_barcode_mapping[cell_barcode].append(cell_barcode)
        else: # if the cell barcode is not a "real" barcode, probably there exsit some "real" barcode with hamming distance of 1
            same_umi_rows = high_confidence_barcodes[high_confidence_barcodes["umi_barcode"] == umi_barcode]
            for umi in list(set(high_confidence_barcodes["umi_barcode"].tolist())):
                #iterating over all the umi barcodes is the data frame
                if align_func.hamming_distance(umi,umi_barcode) == 1:
                    # this umi is very similar to our current umi_barcode
                    same_umi_rows = same_umi_rows.append(high_confidence_barcodes[high_confidence_barcodes["umi_barcode"] == umi])
            for index,row_umi in same_umi_rows.iterrows():
                # iterating over the cell barcodes in the data frame we constructed
                other_barcode = row_umi["cell_barcode"]
                if other_barcode in map_cell_to_barcode['Cell_barcode'].tolist() and align_func.hamming_distance(other_barcode,
                                                                                                      cell_barcode) == 1:
                    cell_barcode_mapping[other_barcode].append(cell_barcode)
                    if other_barcode not in checked_cells.keys():
                        checked_cells[other_barcode].append(umi_barcode)
                    break
    return checked_cells, cell_barcode_mapping

def cell_to_tcr(dir_name):
    output_dir = os.path.join(dir_name,"IgBLAST_output")
    io.makeOutputDir(output_dir)
    for file in os.listdir(dir_name):
        if "fasta" in file:
            run_igblast(os.path.join(dir_name,file),output_dir)


def run_igblast(fasta_file,output_dir):
    print("##Running IgBLAST##")
    igblast = "/apps/RH7U2/gnu/ncbi-igblast/1.7.0/bin/igblastn"

    igblast_index_location = "/home/labs/amit/diklag/tracer/resources/Hsap/igblast_dbs"
    imgt_seq_location = "/home/labs/amit/diklag/tracer/resources/Hsap/raw_seqs"

    file_name = os.path.basename(fasta_file).split(".")[0]
    igblast_out = "{output_dir}/{cell_name}.IgBLASTOut".format(
        output_dir=output_dir, cell_name=file_name)

    ig_seqtype = 'TCR'

    databases = {}
    for segment in ['V', 'D', 'J']:
        databases[segment] = "{}/{}_{}.fa".format(igblast_index_location, ig_seqtype,
                                                  segment)

    # Lines below suppress Igblast warning about not having an auxliary file.
    # Taken from
    # http://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
    DEVNULL = open(os.devnull, 'wb')

    command = [igblast, '-germline_db_V', databases['V'],
               '-germline_db_D', databases['D'],
               '-germline_db_J', databases['J'], '-domain_system',
               'imgt', '-organism', 'human',
               '-ig_seqtype', ig_seqtype, '-show_translation',
               '-num_alignments_V', '5',
               '-num_alignments_D', '5', '-num_alignments_J', '5',
               '-outfmt', '7', '-query', fasta_file]


    with open(igblast_out, 'w') as out:
        # print(" ").join(pipes.quote(s) for s in command)
        subprocess.check_call(command, stdout=out, stderr=DEVNULL)

    DEVNULL.close()


