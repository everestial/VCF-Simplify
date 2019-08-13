    
    
def parse_table_to_vcf(table_to_vcf):
    """ Arguments for converting Table file to VCF """

    table_to_vcf.add_argument(
        "-samples", nargs = "+",
        help="Name of the samples -> "
        "space separated name of the samples that needs to be converted.\n"
            "Default = all"
        "to VCF format",
        default=["all"],
        required=False,
    )

    table_to_vcf.add_argument(
        "-formats", nargs = "+",
        help="Name of the FORMAT tags to write -> " "space separated FORMAT tags name.\n"
            "Default = all",
        default=["all"],
        required=False,
    )

    table_to_vcf.add_argument(
        "-infos", nargs = "+",
        help="Name of the INFO tags to write -> " "space separated INFO tags name.\n"
            "Default = all",
        default=["all"],
        required=False,
    )
    
    genotype_options = ["GT:numeric", "GT:iupac", "PG:iupac", "....."]
    table_to_vcf.add_argument(
        "-GTbase", type=str,
        required=False, nargs = "+",
        default=['GT:numeric'],
        help="Suggest if the genotype fields (GT, PG etc.) are in IUPAC base code.\n" 
        "Default = GT:numeric (i.e assumes 'GT' bases are numeric)\n"
        "Choices : [%s]. " % (", ".join(genotype_options)) + "\n"
        "Multiple option can be suggested for each genotype fields.\n",
        metavar=None
    )
    

## **to do: need to explan this later by borrowing method from 
    # haplotype to VCF python file 
def parse_haplotype_to_vcf(haplotype_to_vcf):
    """ Arguments for converting haplotype file to VCF """    
    haplotype_to_vcf.add_argument(
            "-haplotypeFormat", default = 'iupac',
            help = "report which format (numeric vs. iupac) the haplotype file is in.\n" 
            "Default = iupac", 
            metavar = '')

