import argparse
import os
import sys
from configparser import ConfigParser

from gelSeqLib import plate_to_cells,VDJ_func, io_func
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
                                default='Mmus')

            self.parser.add_argument('fastq1', metavar="<FASTQ1>",
                                help='first fastq file - read1')
            self.parser.add_argument('fastq2', metavar="<FASTQ2>",
                                help='second fastq file - read2')
            self.parser.add_argument('plate_name', metavar="<PLATE_NAME>",help='name of plate for file labels')
            self.parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<plate_name>')

            self.vdj_parser.add_argument('--receptor_name',
                                help="Name of receptor to reconstruct",
                                default='TCR')
            self.vdj_parser.add_argument('--loci',
                                help="Space-separated list of loci to reconstruct for receptor",
                                default=['A', 'B'], nargs='+')

            args = self.parser.parse_args(sys.argv[2:])
            self.plate_name = args.plate_name
            self.fastq1 = args.fastq1
            self.fastq2 = args.fastq2
            self.ncores = str(args.ncores)
            self.species = args.species
            self.resume_with_existing_files = args.resume_with_existing_files
            self.output_dir = args.output_dir
            self.cell_to_path = dict()
            self.receptor_name = args.receptor_name
            self.loci = args.loci
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
            self.cell_to_path = dict()
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')

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


    def split_to_cells(self):
        high_confidence_barcodes = plate_to_cells.filter_abundant_barcodes(self.fastq2)
        self.cell_to_path = plate_to_cells.split_by_cells(high_confidence_barcodes, wells_cells_file, self.output_dir, self.fastq1, self.fastq2)


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
                                     default='Mmus')

            self.parser.add_argument('fasta', metavar="<FASTA>",
                                     help='fasta file')
            self.parser.add_argument('cell_name', metavar="<CELL_NAME>", help='name of cell for file labels')
            self.parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                     help='directory for output as <output_dir>/<cell_name>')

            self.vdj_parser.add_argument('--receptor_name',
                                         help="Name of receptor to reconstruct",
                                         default='TCR')
            self.vdj_parser.add_argument('--loci',
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










