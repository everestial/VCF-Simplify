import time
import sys

"""Step 03 (D): Function for Table To VCF"""


def fnc_table_to_vcf(args):
    print("converting Table file to VCF")

    begin_time = time.time()

    """Assign some input variables. """   
    
    infile = args.inFile
    meta_header = args.vcfHeader
    outfile = args.outVCF
    samples = args.samples
    formats = args.formats
    infos = args.infos
    
    # find the genotype tags that are in iupac bases 
    genotype_is = args.GTbase     
    gt_tag_as_iupac = []
    for gts_tag in genotype_is:
        tag_format = gts_tag.split(':')
        if tag_format[1] == 'iupac':
            gt_tag_as_iupac.append(tag_format[0])


    with open(infile) as tablefile, open(meta_header) as meta_header, open(
        outfile, "w+") as vcf_out:

        """Start reading the haplotype file as generator. This saves memory. """
        for line in tablefile:

            ## find and set the indexes ...
            # ... of pre-fields, INFO, FORMAT and SAMPLE level information

            """ Step 01: The very first line of the file is read;
             - to find the variable name and it's index position in the input file.
             - almost all the variable created downstream are "global variables".
             - SAMPLE level information is automatically identified unless explicitly given.
               The sample names is identified using ":" in the column names, so other names
               should not have ":" at all.
             - FORMAT level tags can also be provided as list, or can be mined automatically
               along with SAMPLE by using ":" matching.
             - All the preHeader tags, ie. CHROM  POS  ID  REF  ALT  QUAL  FILTER are reserved and
               updated by matching the names in text header line.
             """

            # to use the "header" name that have already been taken
            # this will help in finding appropriate "INFO" level tags from the header file
            used_header = []

            if line.startswith("CHROM") or line.startswith("#CHROM"):
                header_line = line.rstrip("\n").split("\t")
                
                
                ################# function 01 ######################
                ## ?? Bhuwan - move this as a preheader to another function and optimize if posible

                if "CHROM" in header_line:
                    contig_idx = header_line.index("CHROM")
                    
                    # update the taken header "labels"
                    used_header += ["CHROM"]
                    
                elif "#CHROM" in header_line:
                    contig_idx = header_line.index("#CHROM")
                                                   
                    # update the taken header "labels"
                    used_header += ["#CHROM"]
                else:
                    print("CHROM field does not exist in the input table file. Update your file")
                    print("Exiting the program")
                    sys.exit(0)

                if "POS" in header_line:
                    pos_idx = header_line.index("POS")
                    used_header += ["POS"]
                else:
                    print("POS field does not exist. Update your file")
                    print("Exiting the program")
                    sys.exit()                

                if "ID" in header_line:
                    id_idx = header_line.index("ID")
                else:
                    id_idx = None
                used_header += ["ID"]

                if "REF" in header_line:
                    ref_idx = header_line.index("REF")
                else:
                    ref_idx == None
                used_header += ["REF"]

                if "ALT" in header_line:
                    alt_idx = header_line.index("ALT")
                else:
                    alt_idx = None
                used_header += ["ALT"]

                if "QUAL" in header_line:
                    qual_idx = header_line.index("QUAL")
                else:
                    qual_idx = None
                used_header += ["QUAL"]

                if "FILTER" in header_line:
                    filter_idx = header_line.index("FILTER")
                else:
                    filter_idx = None
                used_header += ["FILTER"]
                ################## function 01 ends here ######
                
                
                ###############function 02 ##################################
                ## ?? Bhuwan - move this to a function and optimize the process 
                """INFO tags are identified by matching "INFO:" in the column names."""
                infos_in_header = [x for x in header_line if x.startswith("INFO:")]
                all_infos = [x.replace("INFO:", "") for x in infos_in_header]
                
                if len(all_infos) == 0:
                    print("INFO tags are not available.")
                    print("INFO field will be populated with empty '.' value")
                    info_tags = []
                    
                elif infos[0] == 'all':
                    info_tags = all_infos
                    print("Using the following metrics as INFO tags: ")
                    print("    %s" % info_tags)
                else: 
                    info_tags = infos
                    
                    ## find any missing INFO tags or any nonsense tag
                    non_matching_infos = list(set(info_tags) - set(all_infos))
                    non_used_infos = list(set(all_infos) - set(info_tags))
                    
                    if len(non_matching_infos) > 0:
                        print("the following user provided infos are not available in input file")
                        print("  %s" % non_matching_infos)
                        
                    if len(non_used_infos) > 0:
                        print("the following INFO tags won't be put in INFO fields of output VCF")
                        print("  %s" %non_used_infos)
                        
                # also find the position of the info tags on header line
                infos_idx = []
                if len(info_tags) != 0:
                    for inftag in info_tags:
                        infos_idx.append(header_line.index("INFO:" + inftag))
                else:
                    infos_idx = None                      
                ##########################  #######################
                

                ##############function 03 ################################
                ## ?? Bhuwan - move this to a separate function and optimize the process if possible 
                """SAMPLE names and FORMAT tags are identified using ":" delimiter in the column names, 
                    after excluding the INFO fields."""
                possible_samples = [x for x in header_line if ':' in x]
                
                # remove the INFO fields 
                samples_and_formats = list(set(possible_samples) - set(infos_in_header))
                
                # split and set to collect unique sample names and unique format tags 
                # make sure to add process so the order is maintained
                all_samples = list(set([x.split(':')[0] for x in samples_and_formats]))
                all_formats = list(set([x.split(':')[1] for x in samples_and_formats]))
                              
                
                # find the available format tags
                # ?? Bhuwan - write this as a separate function or subfunction
                #### sub function 03 A 
                ### prepare sample names 
                if formats[0]== "all":
                    format_tags = all_formats
                elif len(formats) == 0:
                    print("No format tags available.")
                    format_tags = []  
                else: format_tags =  formats 

                # In the available FORMAT tags, move "GT" field to the beginning.
                if "GT" in format_tags:
                    format_tags.remove("GT")
                    format_tags.insert(0, "GT")
                         
                
                
                ## ?? Bhuwan - write as sub function 03 B 
                ### prepare sample names 
                if samples[0]== "all":
                    sample_names = all_samples
                elif len(samples) == 0:
                    print("No sample available.")
                    sample_names = []  
                else: 
                    sample_names = samples
                    nonsense_sample_names = [x for x in samples if not x in all_samples]
                    if len(nonsense_sample_names) > 0:
                        print("The following sample names %s are not available in table file and not valid." 
                              % nonsense_sample_names)
                        sys.exit(0)                         
                
                used_header += sample_names

                
                

                ### ?? Bhuwan - write as function 04 and optimize 
                """ Now, Read the meta header and add it to the output VCF file. """
                print('\nReading meta header from file "%s" ' % (meta_header.name))

                if meta_header != None:
                    meta_info = meta_header.readlines()
                    # if the meta header has "#CHROM	POS	REF ...." line then delete it
                    if meta_info[-1].startswith("#CHROM\tPOS"):
                        meta_info = "".join(meta_info[:-1]).rstrip("\n")
                    else:
                        meta_info = "".join(meta_info).rstrip("\n")

                else:
                    print("Header with meta information is not provided")
                    print("Exiting the program")
                    sys.exit(0)

                # add meta header to the output VCF file
                meta_info += "\n"
                meta_info += (
                    "\t".join(
                        [
                            "#CHROM",
                            "POS",
                            "ID",
                            "REF",
                            "ALT",
                            "QUAL",
                            "FILTER",
                            "INFO",
                            "FORMAT",
                        ]
                    )
                    + "\t"
                )

                # add SAMPLE fields to output VCF file
                meta_info += "\t".join(sample_names)

                # Finally, write the header part of the output VCF
                vcf_out.write(meta_info + "\n")
                
                ######### function 04 ends here ###########
                

                continue
                # break
            

            """' Now, extract the required data from each of the remaining lines add to output VCF. """
            updated_line = table_to_vcf(
                line,
                contig_idx,
                pos_idx,
                id_idx,
                ref_idx,
                alt_idx,
                qual_idx,
                filter_idx,
                infos_idx,
                info_tags,
                format_tags,
                sample_names,
                gt_tag_as_iupac,
                header_line,
            )
            vcf_out.write(updated_line)
            vcf_out.write("\n")

        print('Elapsed time : "%s".' % (time.time() - begin_time))


"""Function part of Table to VCF """


def table_to_vcf(
    line_in,
    contig_idx,
    pos_idx,
    id_idx,
    ref_idx,
    alt_idx,
    qual_idx,
    filter_idx,
    infos_idx,
    info_tags,
    format_tags,
    sample_names,
    gt_tag_as_iupac,
    header_line,
): 
    

    line = line_in.rstrip("\n").split("\t")
    
    
    if contig_idx is not None:
        chrom = line[contig_idx]
    else:
        chrom = "."
        
    if pos_idx is not None:
        pos = line[pos_idx]
    else:
        pos = "."

    if id_idx is not None:
        ids = line[id_idx]
    else:
        ids = "."
        
    if ref_idx is not None:
        ref = line[ref_idx]
    else:
        ids = "."
        
    if alt_idx is not None:
        alt = line[alt_idx]
    else:
        alt = "."

    if qual_idx is not None:
        qual = line[qual_idx]
    else:
        qual = "."

    if filter_idx is not None:
        filter_ = line[filter_idx]
    else:
        filter_ = "."

    # Update "info tags and value". This is little complex
    if info_tags != []:
        info_ = []
        for ith, itemi in enumerate(info_tags):
            tag_val = "=".join([itemi, line[infos_idx[ith]]])
            info_.append(tag_val)
        info_ = ";".join(info_)
    elif info_tags == []:
        info_ = "."

    # write the tags names of the FORMAT column
    if format_tags != None:
        format_ = ":".join(format_tags)
    else:
        format_ = "."

    # update the output line
    line_out = (
        "\t".join([chrom, pos, ids, ref, alt, qual, filter_, info_, format_]) + "\t"
    )

    # Further update the SAMPLE-to-FORMAT values
    # pass the line to another function
    format_to_sample_vals = update_sample_format(
        line, ref, alt, sample_names, format_tags, header_line, gt_tag_as_iupac
    )
    line_out = line_out + format_to_sample_vals

    return line_out


""" Function part of Table to VCF """


def update_sample_format(
    line, ref, alt, sample_names, format_tags, header_line, gt_tag_as_iupac
):
    
    # The "line" variable is passed into this function.
    # The global variables are "genotype_is", "sample_names" and "format_tags"

    # to store updated line
    format_sample_line = []
    all_alleles = [ref] + alt.split(",")

    for namex in sample_names:
        namex_vals = []
        for tagx in format_tags:
            sample_format_tag = namex + ":" + tagx
            sample_format_idx = header_line.index(sample_format_tag)
            sample_format_val = line[sample_format_idx]

            """ further update the sample:format value if GT in table is as IUPAC base """
            if tagx in gt_tag_as_iupac:
                if (
                    sample_format_val == "."
                    or sample_format_val == "./."
                    or sample_format_val == ".|."
                ):
                    continue

                elif "/" in sample_format_val:
                    sample_format_val = sample_format_val.split("/")

                    sample_format_val = [
                        all_alleles.index(sample_format_val[0]),
                        all_alleles.index(sample_format_val[1]),
                    ]

                    sample_format_val = "/".join(str(xth) for xth in sample_format_val)

                elif "|" in sample_format_val:
                    sample_format_val = sample_format_val.split("|")

                    sample_format_val = [
                        all_alleles.index(sample_format_val[0]),
                        all_alleles.index(sample_format_val[1]),
                    ]

                    sample_format_val = "|".join(str(xth) for xth in sample_format_val)

            namex_vals.append(sample_format_val)

        format_sample_line.append(":".join(namex_vals))

    sample_format_final = "\t".join(format_sample_line)

    return sample_format_final
