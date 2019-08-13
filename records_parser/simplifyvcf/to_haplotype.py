import time
import sys
import itertools

from metadata_parser.vcf_metadata_parser import VcfReadMetaData
from records_parser.vcf_records_parser import VcfReadIndividualRecord

## to do ?? Bishwa:
# have a method to replace args.PI, args.PG with any format tags (not just GT and PGs)
# add method to return genotype as numeric vs IUPAC 
# method to write unphased genotypes 



"""Step 03 (A): Function for VCF To Haplotype """

# print('\navailable paths to haplotype level: %s' % sys.path)
print()


def fnc_vcf_to_haplotype(args):

    print("Creating Haplotype file from VCF file")
    start_time01 = time.time()
    pi_tag = args.PI
    pg_tag = args.PG
    
    if args.includeUnphased == '0' or args.includeUnphased == 'no':
        include_unphased = 'no'
    else: include_unphased = 'yes'
    
    #update default argument for VCF to haplotype
    if args.GTbase == ['GT:numeric']: 
        args.GTbase = ['PG:iupac']

    # extract metadata and raw header from the input VCF
    metadata, only_header, record_keys = VcfReadMetaData(args.inVCF).read_metadata()
    sample_ids = [x["name"] for x in metadata["samples"]]

    print("%i samples found" % len(sample_ids))
    print()

    # method to write raw header to user provided filename
    if args.outHeaderName:
        print("Writing the header to a separate output file.")
        with open(args.outHeaderName, "w") as vcfheader:
            # vcfheader.write(vcf_file.raw_header)
            vcfheader.write(only_header)
    else:
        print("Skipping the header.")
    print()
    
    # find the requested output format for the genotypes of interest 
    gtbase = parse_genotypes_format(args.GTbase)    
    for gts in gtbase:
        print("sample genotypes tag '%s' are written as '%s' bases" % (gts[0], gts[1]))

    with open(args.outFile, "w") as write_block, open(args.inVCF) as invcf:

        ## write the columns names for the haplotype file
        column_names = ["CHROM", "POS", "RefAndAlts"]
        for name in sample_ids:
            # to store "phased index"
            column_names.append(name + ":PI")

            # to store the "phased genotype" as IUPAC bases
            column_names.append(name + ":PGal")

        # write the header of the output file
        # print('\t'.join(column_names), file=write_block)
        write_block.write("\t".join(column_names) + "\n")

        # to keep track of which chromosome/contig is being processed
        chr_on_process = ""

        print("parsing records ... ")
        # skips line starting with '#'
        records_gen = itertools.dropwhile(lambda line: line.startswith("#"), invcf)
        for records in records_gen:
            mapped_record = VcfReadIndividualRecord(
                recordKeys=record_keys, inRecord=records, 
                sample_names=sample_ids, gtbase_is = gtbase
            ).read_vcfRecord()

            """Now, this mapped records can be mined for extracting required data."""
            contig = mapped_record["CHROM"]
            pos = mapped_record["POS"]
            ref_ = [mapped_record["REF"]]
            alt_ = mapped_record["ALT"].split(",")
            all_alleles = ref_ + alt_

            # find which chr is in the process
            if chr_on_process != contig:
                print("contig %s is being processed" % str(contig))
                print()
                chr_on_process = contig

            # write data to the output file
            # ('2', 15881018, 'G', ['G', 'A', 'C'], [0.0, 1.0])
            write_block.write("\t".join([contig, pos, ",".join(all_alleles)]))
            
            #print('mapped record')
            #print(mapped_record)

            # this is returned as array. for now it is just easy to mine data from the "all_alleles" variable
            # ** so this may be used in the future
            #gt_bas = [mapped_record[sample]["GT"] for sample in sample_ids]

            # start mining the phased_index (PI) and phased_genotype (PG) values for each sample
            # this outputs a list of PI and PG for all available sample in the VCF file
            if pi_tag == "CHROM":
                pi_values = [contig for _ in range(len(sample_ids))]
            else:
                #pi_values = [mapped_record[sample]["PI"] for sample in sample_ids]
                pi_values = [mapped_record[sample][pi_tag] for sample in sample_ids]

            if pg_tag == "GT":
                pg_values = [mapped_record[sample]["GT"] for sample in sample_ids]
            else:
                pg_values = [mapped_record[sample].get(pg_tag, '.') for sample in sample_ids]

            # incase the "PI" and "PG" tags are missing in the VCF file
            # e.g: sometimes "bcftools" merged variants have records with missing the PI and PG tags.
            if len(pi_values) == 0:
                pi_values = ["." for _ in range(len(sample_ids))]

            if len(pg_values) == 0:
                pg_values = ["." for _ in range(len(sample_ids))]

            """ now, start extracting sample level information """
            ## ?? to do - this mapping process requires some optimization in the future
            ## Bhuwan and I will do it together. Do not worry about it right now.
            for ith, pi_value in enumerate(pi_values):

                pg_allele = compute_pg_allele(
                    pg_values[ith], all_alleles, include_unphased)

                write_block.write("\t" + "\t".join([pi_value, pg_allele]))
            write_block.write("\n")

    print("elapsed time: ", time.time() - start_time01)
    print()


def compute_pg_allele(pg_val, all_alleles, include_unphased):
    if pg_val in (".", "./.", ".|."):
        return "."
    else:
        sep = "/" if "/" in pg_val else "|"
        if (sep == "/") and (include_unphased == "no"):
            return "."
        else:
            return pg_val            
            #a, b = pg_val.split(sep)
            #return all_alleles[int(a)] + sep + all_alleles[int(b)]
        
def parse_genotypes_format(gtbase):
    gt_output = []
    for gts in gtbase:
        gt_output.append(gts.split(':'))        
    return tuple(gt_output)
