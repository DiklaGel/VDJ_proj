from __future__ import print_function
import matplotlib as mpl

mpl.use('pdf')
import argparse
import sys

from plate import Plate_Task, VDJ_Plate_Task


def launch():
    parser = argparse.ArgumentParser(
        description='gelSeq: reconstruction of TCR sequences from single-cell RNAseq data',
        usage=''' gelSeq <mode> [<args>]

              Modes are :

              - plate: process regular plate fastq file - split them by cell barcodes
              - vdj: assemble TCR sequences from single-cell RNA-sequencing reads

              use gelSeq <mode> -h for specific help
              ''')
    parser.add_argument('mode', metavar="<MODE>", help='gelSeq mode (plate, vdj)',
                        choices=['plate', 'vdj'])
    args = parser.parse_args(sys.argv[1:2])

    task_mapper = {
        'plate': Plate_Task,
        'vdj': VDJ_Plate_Task,
    }

    if args.mode not in task_mapper:
        print('Unrecognised mode')
        parser.print_help()
        exit(1)

    Task = task_mapper[args.mode]
    Task().run()
