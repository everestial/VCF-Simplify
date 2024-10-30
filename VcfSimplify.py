import time
import sys
import argparse
import os

import metadata_parser
import records_parser

from assign_task import sub_argparsers, author_name
from assign_task.perform_operation import vcf_solver

## fixing some import issues here to absolute imports. Can be fixed later.
## append top level directory path to assist in absolute imports
sys.path.append(os.path.dirname(metadata_parser.__file__))
sys.path.append(os.path.dirname(records_parser.__file__))

""" Purpose of the program: mine the data from the VCF files and convert it into Haplotype file.
The output file consists of ReadBackPhased genotypes as PI (block keys) and PG (genotype values). 
Other GT with no phased
state may be extracted as well."""


def main():

    # print author name 
    author_name.print_author_name()

    # define argument variables
    main_parser = argparse.ArgumentParser(
            prog="VCF-Simplify", formatter_class=argparse.RawTextHelpFormatter)
    
    # Create sub_parsers (for SimplifyVCF, BuildVCF and ViewVCF)
    subparsers = main_parser.add_subparsers(help="Choose one of the following method.")

    # Add three subparsers defined in sub_parsers file
    sub_argparsers.view_argparser(subparsers)
    sub_argparsers.simplify_argparser(subparsers)
    sub_argparsers.build_argparser(subparsers)

    # creating an argument variable to handle the loaded argument
    args = main_parser.parse_args()

    """ Step 02: Based on positional arugments and task go to specific function """
    task = sys.argv[1]
    vcf_solver(task, args)
    # assign_task.perform_operation.vcf_solver(task, args)


if __name__ == "__main__":
    main()
