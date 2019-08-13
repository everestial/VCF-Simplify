
import argparse
from . import simplify_arg, build_arg



## *to do: add gt base as numeric vs IUPAC, 
    # control sample names (prefix, suffix, match)?? 
    # add INFO:tag in table output file 
    # write unit test
    # application optimization (for loop optimization)
    # cythonize etc. 


# to format the help message 
class BlankLinesHelpFormatter (argparse.RawTextHelpFormatter):
    def _split_lines(self, text, width):
        return super()._split_lines(text, width) + ['']



def view_argparser(subparsers):

    """Step 01 - Parser A: Sub parser to view and extract metadata information from the VCF file. """
    parser_a = subparsers.add_parser(
        "ViewVCF", help = "  View and extract metadata from the VCF.", 
        formatter_class = BlankLinesHelpFormatter
    )

    # available parameters 
    #** if possible add more in the future 
    metadata_keys = [
        "VCFspec",
        "reference",
        "contig",
        "samples",
        "INFO",
        "FORMAT",
        "FILTER",
        "GATKCommandLine",
        "GVCFBlock",
    ]

    # available output data types 
    data_types = ["table", "json", "dict"]

    # add argument keys within "ViewVCF"
    parser_a.add_argument("-inVCF", required=True, help="Sorted vcf file.")
    parser_a.add_argument(
        "-outFile",
        required = False,
        help = "Name of the output file without file extension.",
    )
    parser_a.add_argument(
        "-outType",
        required = False,
        choices = data_types,
        nargs = "+",
        help = "Space separated list of output data types.\n" 
        "Multiple types can be requested.",
    )
    parser_a.add_argument(
        "-metadata",
        required=False,
        nargs="+",
        help="Space separated list of metadata of interest." + "\n"
        + "Allowed values are: " + "\n"
        + "  " + ", ".join(metadata_keys) + ".\n" 
        + "Multiple choices can be requested.",
        metavar=None,
    )


def simplify_argparser(subparsers):

    """Step 01 - Parser B: Sub parser to extract record information from the VCF file. """
    parser_b = subparsers.add_parser(
        "SimplifyVCF", help="  Simplify VCF -> to {haplotype, table} data.", 
        formatter_class=BlankLinesHelpFormatter
    )

    # add argument keys within SimplifyVCF
    parser_b.add_argument(
        "-toType",
        help="Type of the output file.",
        nargs=1,
        choices=["haplotype", "table"],
        required=True,
    )
    parser_b.add_argument("-inVCF", help="Sorted vcf file.", required=True)
    parser_b.add_argument("-outFile", help="Name of the output file.", required=True)
    parser_b.add_argument(
            "-outHeaderName", 
            help = "Write the VCF raw METADATA HEADER to a separate output file.\n" 
            "Default: no output.", 
            required=False)
    
    
    genotype_options = ["GT:numeric", "GT:iupac", "PG:iupac", "....."]
    parser_b.add_argument(
        "-GTbase", type=str,
        required=False, nargs = "+",
        default=['GT:numeric'],
        help="Write the genotype (GT, PG etc.) field as IUPAC base code.\n" 
        "Default = GT:numeric (i.e write 'GT' bases as numeric)\n"
        "Choices : [%s]. " % (", ".join(genotype_options)) + "\n"
        "Multiple option can be requested for each genotype fields.\n",
        metavar=None
    )


    """Parser A (01) : sub arguments for - From VCF to Haplotype """
    vcf_to_haplotype = parser_b.add_argument_group(
            title = 'Additional arguments for "VCF To -> Haplotype"')
    
    """Part A (02) : sub arguments for - From VCF to Table"""
    vcf_to_table = parser_b.add_argument_group(
            title = 'Additional arguments for "VCF To -> Table"')
    
    # now, pass the parser to respective function that collects arguments .. 
     # .. for each group     
    simplify_arg.parse_vcf_to_haplotype(vcf_to_haplotype)    
    simplify_arg.parse_vcf_to_table(vcf_to_table)
    

def build_argparser(subparsers):
    """Step 01 - Part C: Sub parser to build VCF file from table or haplotype data."""
    parser_c = subparsers.add_parser(
        "BuildVCF", help="  Build VCF <- from {haplotype, table} data.", 
        formatter_class=BlankLinesHelpFormatter
    )

    # add argument keys within BuildVCF
    parser_c.add_argument(
        "-fromType",
        required=True,
        help="Type of the input file the VCF is being prepared from. ", 
        choices = ['haplotype', 'table']
    )
    parser_c.add_argument(
        "-inFile",
        required=True,
        help="Sorted table or haplotype file.\n"
        "Note: \n"
        "  Haplotype file should be in the format created by 'phase-Stitcher' or 'phase-Extender'.\n" 
        "  Table file should be in the format created by 'VCF-Simplify'\n" 
        "Only wide?? table is supported for now.",
    )
        

    parser_c.add_argument("-outVCF", help="Name of the output VCF file.", 
                          required=True)
    
    parser_c.add_argument(
        "-vcfHeader",
        required=True,
        help="A custom VCF header to add to the VCF file.\n"
        "The VCF header should not contain the line with #CHROM .... \n"
        "#CHROM ... line is auto populated while creating the VCF file." ,
    )
    

    """Part B (01) : sub arguments for - From Table To VCF. """
    table_to_vcf = parser_c.add_argument_group(
        'Additional arguments for "Table To VCF"')   
    
    """Part B (02) : sub arguments for - From haplotype To VCF. """
    haplotype_to_vcf = parser_c.add_argument_group(
        'Additional arguments for "Haplotype To VCF"'
    )
    
    # now, pass the parser to respective function that collects arguments .. 
     # .. for each group
     
    build_arg.parse_table_to_vcf(table_to_vcf)
    build_arg.parse_haplotype_to_vcf(haplotype_to_vcf)
    
    
    
    
    