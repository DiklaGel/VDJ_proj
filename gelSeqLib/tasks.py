import argparse
import os
import sys
import warnings
import pickle
import subprocess
sys.path.insert(0, '/home/labs/amit/diklag/python_libs/python_lsf_wrapper/')
from LSF import LSF, wait_for_jobs

from configparser import ConfigParser

from gelSeqLib import plate_to_cells,VDJ_func, io_func, align_func
from gelSeqLib import base_dir, wells_cells_file



class Task:
    base_parser = argparse.ArgumentParser(add_help=False,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_parser.add_argument('--ncores', '-p', metavar="<CORES>",
                             help='number of processor cores to use', type=int,
                             default=1)
    base_parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>",
                             help='config file to use',
                             default='~/.tracerrc')

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
            config_file = '~/.tracerrc'
        config_file = os.path.expanduser(config_file)
        if not os.path.isfile(config_file):
            print(
                "Config file not found at ~/.tracerrc. Using default tracer.conf in repo...")
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
        self.output_dir = root_output_dir + "/" + self.plate_name

        io_func.makeOutputDir(self.output_dir)

        # Perform core functions
        self.split_to_cells()

        #if self.full:



    def split_to_cells(self):
        high_confidence_barcodes = plate_to_cells.filter_abundant_barcodes(self.fastq2)
        high_confidence_barcodes.to_csv(self.output_dir+"/high_conf.csv")
        plate_to_cells.split_by_cells(high_confidence_barcodes, wells_cells_file, self.output_dir, self.fastq1, self.fastq2)


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
        if not os.path.exists(self.output_dir+'/reads'):
            print("dir {}/reads not exists ".format(self.output_dir))
            return


        data_dirs = ['IgBLAST_output',
                     'unfiltered_{receptor}_seqs'.format(
                         receptor=self.receptor_name),
                     'expression_quantification',
                     'filtered_{receptor}_seqs'.format(
                         receptor=self.receptor_name)]
        for d in data_dirs:
            io_func.makeOutputDir("{}/{}".format(self.output_dir, d))

        cell = self.ig_blast()

        self.print_cell_summary(
            cell,
            "{output_dir}/unfiltered_{receptor}_seqs/unfiltered_{receptor}s.txt".format(
                output_dir=self.output_dir,
                receptor=self.receptor_name),
            self.receptor_name, self.loci)

        with open(
                "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}.pkl".format(
                    output_dir=self.output_dir,
                    cell_name=self.cell_name,
                    receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

        summary = cell.choose_recombinants()

        filt_file = "{output_dir}/filtered_{receptor}_seqs/filtered_{receptor}s.txt".format(
                output_dir=self.output_dir,
                receptor=self.receptor_name)
        self.print_cell_summary(
            cell,filt_file,self.receptor_name, self.loci)

        cdr3_consensus, consensus = self.create_cdr3_consensus(cell, filt_file)

        with open(
                "{output_dir}/filtered_{receptor}_seqs/{cell_name}.pkl".format(
                    output_dir=self.output_dir,
                    cell_name=cell.name,
                    receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

        if len(summary[self.receptor_name]) != 0:
            with open(
                    "{output_dir}/filtered_{receptor}_seqs/{cell_name}.txt".format(
                        output_dir=self.output_dir,
                        cell_name=cell.name,
                        receptor=self.receptor_name), 'w') as pf:
                for locus in self.loci:
                    pf.write("locus = " + locus + '\n')
                    pf.write(summary[self.receptor_name][locus] +'\n')
                    pf.write("CDR3:" + cdr3_consensus + '\n')
                    pf.write("CDR3_freq:" + str(consensus) + '\n')

    def create_cdr3_consensus(self, cell, filt_file):
        cdr3_list = subprocess.getoutput("grep ^CDR3: %s" % (filt_file)).split("\n")
        cdr3_list = ["> \n" + x[6:] for x in cdr3_list if "Couldn" not in x]
        cdr3_list.sort()
        print(cdr3_list)
        cdr3_file = "{output_dir}/filtered_{receptor}_seqs/{cell_name}_cdr3.fasta".format(
            output_dir=self.output_dir,
            cell_name=cell.name,
            receptor=self.receptor_name)
        with open(cdr3_file, 'w') as pf:
            for cdr3 in cdr3_list:
                pf.write(cdr3 + '\n')
        junk_file = open("junk_file", 'w')
        al_file = align_func.clustalo_align(cdr3_file, junk_file)
        os.remove("junk_file")
        cdr3_consensus, freq = align_func.make_consensus(al_file, "fasta")
        consensus = [(cdr3_consensus[i],freq[i]) for i in range(0,len(cdr3_consensus))]
        return cdr3_consensus, consensus

    def ig_blast(self):
        igblastn = self.get_binary('igblastn')

        # Reference data locations
        igblast_index_location = self.get_index_location('igblast_dbs')
        imgt_seq_location = self.get_index_location('raw_seqs')

        igblast_seqtype = self.config.get('IgBlast_options', 'igblast_seqtype')

        # IgBlast of assembled contigs
        VDJ_func.run_IgBlast(igblastn, self.fasta, self.receptor_name, self.loci,
                                self.output_dir, self.cell_name, self.species,
                                igblast_index_location,
                                igblast_seqtype,self.resume_with_existing_files)
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
                for r in rs:
                    out_file.write(r.get_summary())
                    out_file.write("\n\n")

        out_file.close()






