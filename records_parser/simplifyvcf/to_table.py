import time
import sys
import itertools
import warnings

from metadata_parser.vcf_metadata_parser import VcfReadMetaData
from records_parser.vcf_records_parser import VcfReadIndividualRecord

"""Step 03 (B) : Function for VCF To Table."""


def fnc_vcf_to_table(args):
    # extract metadata and raw header from the input VCF
    metadata, only_header, record_keys = VcfReadMetaData(args.inVCF).read_metadata()


    ## Alternative way to access the class attributes 
    ## read_metadata should return self instead of metadata, only_header, record_keys for the following to work 
    # vcf_obj_meta = VcfReadMetaData(args.inVCF).read_metadata()
    # metadata = vcf_obj_meta.metadict
    # record_keys = vcf_obj_meta.record_keys
    # only_header = vcf_obj_meta.raw_header


    all_info = [info["ID"] for info in metadata["INFO"]]
    all_samples = [x["name"] for x in metadata["samples"]]
    all_format = [format_["ID"] for format_ in metadata["FORMAT"]]
    all_preheader = record_keys[0:7]

    print("%i samples found" % len(all_samples))
    print()

    # method to write raw header to user provided filename
    if args.outHeaderName: 
        print("Writing the header to a separate output file.")
        with open(args.outHeaderName, "w") as vcfheader:
            vcfheader.write(only_header)
    else:
        print("Skipping the header.")
    print()

    # Step 02: Now, pipe the "input arguments" to a variable
    pre_header = args.preHeader

    if args.mode == "1" or args.mode == "long":
        mode = "long"
        print("Writing table in a long format.")
        print()
    else:
        mode = "wide"
        print("Writing table in a wide format.")
        print()

    # find the requested output format for the genotypes of interest 
    gtbase = parse_genotypes_format(args.GTbase)    
    for gts in gtbase:
        print("sample genotypes tag '%s' are written as '%s' bases" % (gts[0], gts[1]))


    # Step 03: Read vcf file using cyvcf2 and start mining the data
    start_time01 = time.time()
    # with open("simplified_vcf.txt", 'w') as write_block:
    with open(args.outFile, "w") as write_block, open(args.inVCF) as invcf:

        # get field values according to arguments and all possible tags
        my_preheader = process_fields(args.preHeader, all_preheader, argument_flag = '-preHeader')
        my_infos = process_fields(args.infos, all_info, argument_flag = '-infos')
        my_formats = process_fields(args.formats, all_format, argument_flag = '-formats')

        matches = ["prefix:", "suffix:", "match:"]
        if any(elem.startswith(x) for x in  matches for elem in args.samples):
            my_samples = find_sample_name_by_stringmatch(args.samples, all_samples)
        else:
            my_samples = process_fields(args.samples, all_samples, argument_flag = '-samples')
        #print(my_samples)
        
        # ensure that sample names and format tags are in sync
        check_sample_and_format(my_samples, my_formats)
        
        # create a list to store the column header of the output table file
        output_header = []

        # Step 03-A: now, write the appropriate front part of the header
        output_header += my_preheader

        # Step 03-B: add INFO tags of interest to the output header
        # also add INFO: prefix to info tags
        my_infos_column = ['INFO:' + x for x in my_infos] 
        output_header += my_infos_column

        ## Step 03-C: first add tags from the "FORMAT" of interest to the header
        # *Note: the format tags should however be combined with sample name
        # as SAMPLE:FORMAT tag

        ## Step 03-D: Decide between "long" vs. "wide" way of representing FORMAT field values
        # Now, simplify the TAGS from the FORMAT field by assigning it to each SAMPLE names
        # add the tags from the "FORMAT" field as it is. Each sample is represented in different line
        # if len(my_samples) != 0:
        if mode == "long":
            output_header.append("SAMPLE")
            output_header.extend(my_formats)

            # write the final header to the output file
            # print('\t'.join(output_header), file=write_block)
            write_block.write("\t".join(output_header) + "\n")

        elif mode == "wide":
            # add the tags from the "FORMAT" field as suffix to the sample of interest (on output header)
            for name in my_samples:
                for tags in my_formats:
                    output_header.append(name + ":" + tags)

            # write the final header to the output file
            # print('\t'.join(output_header), file=write_block)
            write_block.write("\t".join(output_header) + "\n")

        print()
        print("- Output header of the simplified VCF file:")
        print(output_header)
        print()

        chr_on_process = ""

        """ Step 04: now, start parsing the VCF file using cyVCF2 and add the data for each header fields"""
        print("Reading the input vcf file ... ")
        print()

        # skips line starting with '#'
        records_gen = itertools.dropwhile(lambda line: line.startswith("#"), invcf)
        for records in records_gen:
            line_to_write = []  # create new emtpy variable

            mapped_record = VcfReadIndividualRecord(
                recordKeys=record_keys, inRecord=records, 
                sample_names=all_samples, gtbase_is = gtbase
            ).read_vcfRecord()
            

            """Now, this mapped records can be mined for extracting required data."""
            contig = mapped_record["CHROM"]

            # find which chr is in the process
            if chr_on_process != contig:
                print("Contig %s is being processed ... " % str(contig))
                print()
                chr_on_process = contig

            # Step 04-A : mine the values for 'pre header' of interest
            #line_to_write += [mapped_record.get(prex, '.') for prex in pre_header]
            line_to_write += [mapped_record[prex] for prex in my_preheader]
            #raise KeyError('key does not exist')                


            # Step 04-B: mine the values for the INFO tags of interest
            infos_to_write = process_info(my_infos, mapped_record)
            line_to_write += infos_to_write

            # Step 04-C: compute values for the FORMAT fields of interest for each SAMPLE names of interest
            # so, we need to use both format_fields and sample_names together
            # and pass it to a defined function
            if mode == "wide":
                format_to_write = process_format_wide(
                    mapped_record, my_samples, my_formats
                )

                line_to_write += format_to_write
                write_block.write("\t".join(line_to_write) + "\n")

            elif mode == "long":         
                for sample in my_samples:
                    newline_to_write = []
                    
                    format_to_write = process_format_long(
                        mapped_record, sample, my_formats)
                    newline_to_write = line_to_write + format_to_write

                    write_block.write("\t".join(newline_to_write) + "\n")
                    
            # write_block.write('\n')
    print("elapsed time: ", time.time() - start_time01)
    print("Completed converting the VCF file to table output.")


def process_fields(given_tags, all_fields, argument_flag):
    if given_tags[0] == "all":
        return all_fields
    elif given_tags[0] == "0":
        return []
    else:
        #print(f'given_tags: {given_tags}') 
        #print(f'all_tags: {all_fields}')        
        if all(elem in all_fields for elem in given_tags):
            return given_tags            
        else:
            non_matching_key = [x for x in given_tags if x not in all_fields]
            #print('non matching keys\n', non_matching_key, argument_flag)
            
            ## ?? Bhuwan - this warning needs further rendering 
            warnings.warn(
                    (f"\n    The provided {non_matching_key} for '{argument_flag}' is not present in VCF metadata.\n" 
                    "    Please make sure your vcf file is valid or your input {given_tags} is valid."), stacklevel=4  )
            
            if argument_flag == '-infos' or argument_flag == '-formats':
                print("    The %s for '%s' will be populated with '.' if not present in vcf records." 
                      % (non_matching_key, argument_flag))
                return given_tags
            
            ## only exit when for the following arugments 
            elif argument_flag == '-preHeader' or argument_flag == '-samples':
                print("remove or fix the %s tags from the argument '%s'" 
                      % (non_matching_key, argument_flag))
                sys.exit(0)
            

def process_info(info_fields, mapped_record):    
    """ map and write info """
    
    infos_to_write = [
        mapped_record["INFO"].get(inftags, ".") for inftags in info_fields
    ]
    return infos_to_write

    # now, for each SAMPLE compute and write the FORMAT's field of interest


def process_format_wide(mapped_record, my_samples, my_formats):

    format_to_write = []
    for name in my_samples:
        format_to_write += [
            mapped_record[name].get(fmtags, ".") for fmtags in my_formats
        ]

    # ?? to do (Bishwa) - add method to convert GT tags numeric values to IUPAC bases
    return format_to_write


def process_format_long(mapped_record, sample, my_formats):
    format_to_write = [sample]
    format_to_write += [mapped_record[sample].get(fmtags, ".") for fmtags in my_formats]
    # ?? to do (Bishwa) - add method to convert GT tags numeric values to IUPAC bases

    return format_to_write


def find_sample_name_by_stringmatch(soi, all_samples):
    # for situation when sample name are given in prefix, suffix, or match
    # e.g: -sample MA605 prefix:ms match:A6
    
    selected_samples = []
    
    matches = ["prefix:", "suffix:", "match:"]
    
    for name in soi:
        if name.startswith("prefix:"):
            nameprefix = name.lstrip("prefix:")
            selected_samples += [x for x in all_samples if x.startswith(nameprefix)]
            
        elif name.startswith("suffix:"):
            namesuffix = name.lstrip("suffix:")
            selected_samples += [x for x in all_samples if x.endswith(namesuffix)]
            
        elif name.startswith("match:"):
            namematch = name.lstrip("match:")
            selected_samples += [x for x in all_samples if namematch in x]
            
        else:
            selected_samples += [name]
    
    # to make sure the sample names are unique, but will randomize the order of the names
    selected_samples = set(selected_samples)
    
    # this ensures that order of sample name is returned according to order of file sample names        
    sample_set = [sample for sample in all_samples if sample in selected_samples]
    
    # further check sample names are valid
    sample_set = process_fields(sample_set, all_samples, argument_flag='-samples')    
    
    return sample_set
    

def check_sample_and_format(sample_names, format_tags):
    """ check that format and sample names are raised properly 
        - if no samples are raised no formats can be raised and vice versa. 
    """
    if len(sample_names) == 0 and len(format_tags) != 0:
        print("if no sample names then there should be no format tags")
        print("exiting")
        sys.exit(0)

    if len(sample_names) != 0 and len(format_tags) == 0:
        print("if no format tags then there should be no sample names")
        print("exiting")
        sys.exit(0)
        
        
def parse_genotypes_format(gtbase):
    gt_output = []
    for gts in gtbase:
        gt_output.append(gts.split(':'))
        
    return tuple(gt_output)



