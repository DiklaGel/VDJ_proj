import argparse
import os
import pickle
import subprocess
import sys
import warnings

import pandas as pd
sys.path.insert(0, '/home/labs/amit/diklag/PycharmProjects/VDJ_proj/python_lsf_wrapper/')
from LSF import LSF, wait_for_jobs

from configparser import ConfigParser

from gelSeqLib import plate_to_cells,VDJ_func, io_func
from gelSeqLib import base_dir



class Task:
    base_parser = argparse.ArgumentParser(add_help=False,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_parser.add_argument('--ncores', '-p', metavar="<CORES>",
                             help='number of processor cores to use', type=int,
                             default=1)
    base_parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>",
                             help='config file to use',
                             default='~/.gelseqrc')

    config = None

    def run(self):
        pass

    def get_binary(self, name):
        tool_key = name.lower() + '_path'
        user_path = None
        if self.config.has_option('tool_locations', tool_key):
            user_path = self.resolve_relative_path(
                self.config.get('tool_locations', tool_key))
        return io_func.check_binary(name, user_path)

    def read_config(self, config_file):
        # Read config file
        if not config_file:
            config_file = '~/.gelseqrc'
        config_file = os.path.expanduser(config_file)
        if not os.path.isfile(config_file):
            print(
                "Config file not found at ~/.gelseqrc. Using default gelSeq.conf in repo...")
            config_file = os.path.join(base_dir, 'gelSeq.conf')
        VDJ_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)
        return config

    def resolve_relative_path(self, path):
        if not path.startswith("/"):
            base_directory = os.path.abspath(os.path.dirname(__file__))
            full_path = os.path.normpath(
                "/{}/../{}".format(base_directory, path))
        else:
            full_path = path
        return full_path

    def get_index_location(self, name):
        location = os.path.join(base_dir, 'resources', self.species, name)

        return location

    def get_resources_root(self, species):
        resources_dir = os.path.join(base_dir, 'resources')
        resources_root = os.path.join(resources_dir, species)
        return resources_root


    def get_available_species(self):
        resources_dir = os.path.join(base_dir, 'resources')
        species_dirs = next(os.walk(resources_dir))[1]
        return species_dirs


class Plate_Task(Task):

    def __init__(self, **kwargs):
        if not kwargs:
            # get list of all available species in resources

            self.parser = argparse.ArgumentParser(add_help=True,description="Process fastq files from single plate, split the reads by cell basrcodes",
                parents=[self.base_parser],
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

            self.parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")

            self.parser.add_argument('--species', '-s',
                                help='Species to use for reconstruction',
                                choices=self.get_available_species(),
                                default='Hsap')

            self.parser.add_argument('fastq1', metavar="<FASTQ1>",
                                help='first fastq file - read1')
            self.parser.add_argument('fastq2', metavar="<FASTQ2>",
                                help='second fastq file - read2')
            self.parser.add_argument('plate_name', metavar="<PLATE_NAME>",help='name of plate for file labels')
            self.parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<plate_name>')

            self.parser.add_argument('--receptor_name',
                                help="Name of receptor to reconstruct",
                                default='TCR')
            self.parser.add_argument('--loci',
                                help="Space-separated list of loci to reconstruct for receptor",
                                default=['A', 'B'], nargs='+')
            self.parser.add_argument('--full',
                                help="Continue the full process - after splitting to cells, create new job for each cell",
                                action="store_true")


            args = self.parser.parse_args(sys.argv[2:])
            self.plate_name = args.plate_name
            self.fastq1 = args.fastq1
            self.fastq2 = args.fastq2
            self.ncores = str(args.ncores)
            self.species = args.species
            self.resume_with_existing_files = args.resume_with_existing_files
            self.output_dir = args.output_dir
            # self.cell_to_path = dict()
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            self.full = args.full
            config_file = args.config_file


        else:
            self.plate_name = kwargs.get('plate')
            self.fastq1 = kwargs.get('fastq1')
            self.fastq2 = kwargs.get('fastq2')
            self.ncores = kwargs.get('ncores')
            self.species = kwargs.get('species')
            self.resume_with_existing_files = kwargs.get(
                'resume_with_existing_files')
            self.output_dir = kwargs.get('output_dir')
            # self.cell_to_path = dict()
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')
            self.full = kwargs.get('full')
            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)
        # self.locus_names = ["TCRA", "TCRB"]


    def run(self, **kwargs):
        # Set-up output directories
        root_output_dir = os.path.abspath(self.output_dir)
        io_func.makeOutputDir(root_output_dir)
        self.output_dir = os.path.join(root_output_dir , self.plate_name)
        io_func.makeOutputDir(self.output_dir)

        # Perform core functions

        self.split_to_cells()
        mapper = [_CELLrun(fasta.replace('.fasta',''), self.output_dir + "/" + fasta, self.output_dir)
                  for fasta in os.listdir(self.output_dir) if ".fasta" in fasta]
        for job in mapper:
            job.submit_command(cpu_cores=5, memory=300, queue="new-short")
            # count_jobs(mapper[0:mapper.index(job) + 1])
        wait_for_jobs(mapper)

        df1 = pd.read_csv(self.output_dir + "/final_output.csv")
        list_csv = [self.output_dir + "/" + file for file in os.listdir(self.output_dir) if ".csv" in file and "output" not in file and "high" not in file]
        df2 = pd.DataFrame(columns=["cell_name","V first","V first counts",
                                                                          "V second","V second counts",
                                                                          "D first","D first counts",
                                                                          "D second","D second counts",
                                                                          "J first", "J first counts",
                                                                          "J second", "J second counts",
                                                                          "CDR3 first", "CDR3 first counts",
                                                                          "CDR3 second", "CDR3 second counts"])
        for file in list_csv:
            df2 = df2.append(pd.read_csv(file), ignore_index=True)
        #for file in list_csv:
            #os.remove(file)
        df1 = pd.merge(df1,df2,on="cell_name",how="inner",left_index=False,right_index=False)
        df1 = df1[['well_id','cell_name','#reads','#umi distribution',"V first","V first counts",
                                                                          "V second","V second counts",
                                                                          "D first","D first counts",
                                                                          "D second","D second counts",
                                                                          "J first", "J first counts",
                                                                          "J second", "J second counts",
                                                                          "CDR3 first", "CDR3 first counts",
                                                                          "CDR3 second", "CDR3 second counts"]]
        df1.to_csv(self.output_dir + "/final_table.csv")

        for file in os.listdir(self.output_dir):
            if ".fasta" in file or ".err" in file or ".log" in file or ".pkl" in file:
                os.remove(self.output_dir + "/" + file)


    def split_to_cells(self):
        high_confidence_barcodes = plate_to_cells.filter_abundant_barcodes(self.fastq2)
        high_confidence_barcodes.to_csv(self.output_dir+"/high_conf.csv")
        wells_cells_file = self.config.get('project directories','wells_cells_file')
        plate_to_cells.split_by_cells(self.plate_name,high_confidence_barcodes, wells_cells_file, self.output_dir, self.fastq1, self.fastq2)


class Cell_Task(Task):
    def __init__(self, **kwargs):
        if not kwargs:
            self.parser = argparse.ArgumentParser(add_help=True,
                                                  description="Reconstruct TCR sequences from RNAseq reads for a single cell",
                                                  parents=[self.base_parser],
                                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

            self.parser.add_argument('--resume_with_existing_files', '-r',
                                     help='look for existing intermediate files and use those instead of starting from scratch',
                                     action="store_true")

            self.parser.add_argument('--species', '-s',
                                     help='Species to use for reconstruction',
                                     choices=self.get_available_species(),
                                     default='Hsap')

            self.parser.add_argument('fasta', metavar="<FASTA>",
                                     help='fasta file')
            self.parser.add_argument('cell_name', metavar="<CELL_NAME>", help='name of cell for file labels')
            self.parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                     help='directory for output as <output_dir>/<cell_name>')

            self.parser.add_argument('--receptor_name',
                                         help="Name of receptor to reconstruct",
                                         default='TCR')
            self.parser.add_argument('--loci',
                                         help="Space-separated list of loci to reconstruct for receptor",
                                         default=['A', 'B'], nargs='+')

            args = self.parser.parse_args(sys.argv[2:])
            self.cell_name = args.cell_name
            self.fasta = args.fasta
            self.ncores = str(args.ncores)
            self.species = args.species
            self.resume_with_existing_files = args.resume_with_existing_files
            self.output_dir = args.output_dir
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            config_file = args.config_file

        else:
            self.cell_name = kwargs.get('cell')
            self.fasta = kwargs.get('fasta')
            self.ncores = kwargs.get('ncores')
            self.species = kwargs.get('species')
            self.resume_with_existing_files = kwargs.get(
                'resume_with_existing_files')
            self.output_dir = kwargs.get('output_dir')
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')

            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)
        # self.locus_names = ["TCRA", "TCRB"]

    def run(self, **kwargs):
        '''
        if not os.path.exists(self.output_dir+'/reads'):
            print("dir {}/reads not exists ".format(self.output_dir))
            return

        data_dirs = ['IgBLAST_output','filtered_{receptor}_seqs'.format(
                         receptor=self.receptor_name)]
        for d in data_dirs:
            io_func.makeOutputDir("{}/{}".format(self.output_dir, d))
        '''
        self.collapse_uniq_seqs()

        cell = self.ig_blast()

        summary = cell.choose_recombinants()


        with open(
                "{output_dir}/{cell_name}.pkl".format(
                    output_dir=self.output_dir,
                    cell_name=cell.name,
                    receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

        if len(summary[self.receptor_name]) != 0:
            for locus in self.loci:
                summary[self.receptor_name][locus].to_csv("{output_dir}/{cell_name}_{receptor}_{locus}.csv".format(
                            output_dir=self.output_dir,
                            cell_name=cell.name,
                            receptor=self.receptor_name, locus=locus))


    def collapse_uniq_seqs(self):
        new_fasta =".".join(self.fasta.split(".")[0:len(self.fasta.split(".")) -1]) + "_collapsed.fasta"
        temp_uniq = subprocess.getoutput("""cat %s | awk 'NR%s==0' | sort | uniq -c | sort -n""" % (self.fasta,"%2")).split("\n")
        with open(new_fasta,'w') as fasta_file:
            for row in temp_uniq:
                query_name = ">" + row.strip().split(" ")[0] + "\n"
                seq = row.strip().split(" ")[1] + "\n"
                fasta_file.write(query_name)
                fasta_file.write(seq)
        self.fasta = new_fasta


    def ig_blast(self):
        igblastn = self.get_binary('igblastn')

        # Reference data locations
        igblast_index_location = self.get_index_location('igblast_dbs')
        imgt_seq_location = self.get_index_location('raw_seqs')
        aux_file_location = self.get_index_location('gl.aux')
        igblast_seqtype = self.config.get('IgBlast_options', 'igblast_seqtype')

        # IgBlast of assembled contigs

        VDJ_func.run_IgBlast(igblastn, self.fasta, self.receptor_name, self.loci,
                                self.output_dir, self.cell_name, self.species,
                                igblast_index_location,
                                igblast_seqtype,aux_file_location,self.resume_with_existing_files)
        print()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # cell = io.parse_IgBLAST(self.receptor_name, self.loci, self.output_dir, self.cell_name, imgt_seq_location,
            # self.species, self.seq_method, self.invariant_sequences)
            cell = io_func.parse_IgBLAST(self.receptor_name, self.loci,
                                    self.output_dir, self.cell_name,imgt_seq_location,
                                    self.species, 'imgt',
                                    50)
            '''
            if cell.is_empty:
                self.die_with_empty_cell(self.cell_name, self.output_dir,
                                         self.species)
            '''

        return cell


    def print_cell_summary(self, cell, output_file, receptor_name, loci):
        out_file = open(output_file, 'w')
        out_file.write(
            '------------------\n{name}\n------------------\n'.format(
                name=cell.name))

        # summarise the productive/total recombinants
        for l in loci:
            out_file.write(
                '{receptor}_{locus} recombinants: {summary}\n'.format(
                    receptor=receptor_name, locus=l,
                    summary=cell.summarise_productivity(receptor_name, l)))

        out_file.write('\n\n')

        for l in loci:
            out_file.write(
                "#{receptor}_{locus}#\n".format(receptor=receptor_name,
                                                locus=l))
            rs = cell.recombinants[receptor_name][l]
            if rs is None:
                out_file.write(
                    "No {receptor}_{locus} recombinants found\n\n".format(
                        receptor=receptor_name, locus=l))
            else:
                a = cell.pre_filtering_reads[receptor_name][l]
                out_file.write("# records before filtering: %d \n" % a)

                if cell.after_filtering_reads is not None and len(cell.after_filtering_reads) != 0 :
                    b = cell.after_filtering_reads[receptor_name][l]
                    out_file.write("# records after filtering: %d \n" % b)

                for r in rs:
                    out_file.write(r.get_summary())
                    out_file.write("\n\n")

        out_file.close()


class _CELLrun(LSF):
    def __init__(self, name, fasta, output_dir,loci='B',receptor_name='TCR', species= 'Hsap'):
        """\
        name        - cell name
        fasta       - fasta path of the cell's reads
        output_dir  - dir of plate of origin
        species     - cell species (Human=Hsap, Mouse = Mmus)

        """

        # build the alignment comand
        cell_cmd = "python3.5 /home/labs/amit/diklag/PycharmProjects/VDJ_proj/gelSeq.py cell -s "+ species +" --loci=" + loci + " --receptor_name=" + receptor_name + " " + fasta + " " + name + " " + output_dir

        self.cmd = cell_cmd

        # add error and log files for the lsf output
        self.cmd = " -o " + fasta.replace(".fasta", ".log") + " " + self.cmd
        self.cmd = " -e "+ fasta.replace(".fasta", ".err") + " " + self.cmd





