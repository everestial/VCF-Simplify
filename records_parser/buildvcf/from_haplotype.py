import sys
import time
import re

from metadata_parser.utils import time_memory_track

"""Step 03 (C): Function for Haplotype To VCF"""
@time_memory_track
def fnc_haplotype_to_vcf(infile, meta_header, outfile, hap_format):

    print("converting Haplotype file to VCF")
    begin_time = time.time()
    
    with open(infile) as hapfile, open(meta_header) as meta_header, open(
        outfile, "w+"
    ) as vcf_out:

        """Start reading the haplotype file as generator. This saves memory. """
        for line in hapfile:
            if line.startswith("CHROM") or line.startswith("#CHROM"):
                header_line = line.rstrip("\n").split("\t")

                if "CHROM" in header_line:
                    contig_idx = header_line.index("CHROM")
                elif "#CHROM" in header_line:
                    contig_idx = header_line.index("#CHROM")
                else:
                    print("CHROM field does not exist. Update your file")
                    print("Exiting the program")
                    sys.exit()

                if "POS" in header_line:
                    pos_idx = header_line.index("POS")
                else:
                    print("POS field does not exist. Update your file")
                    print("Exiting the program")
                    sys.exit()

                if "RefAndAlts" in header_line:
                    all_alleles_idx = header_line.index("RefAndAlts")
                else:
                    print(
                        '"RefAndAlts" field not available in input file. Update your file'
                    )
                    print("Exiting the program")
                    sys.exit()

                """Finally find available SAMPLE names and it's FORMAT tags"""
                sample_namex = []
                for itemx in header_line:
                    if ":" in itemx:
                        sample_namex.append(itemx.split(":")[0])
                sample_names_unordered = list(set(sample_namex))
                
                ## to get the proper order in which sample names appear 
                sample_names_in_order = [x.split(':')[0] for x in header_line if ':' in x]
                
                sample_namex = []
                for name_ in sample_names_in_order:
                    if name_ not in sample_namex and name_ in sample_names_unordered:
                        sample_namex.append(name_)
                        
                

                # assign FORMAT tags - keeping it fixed
                format_tagx = ["GT", "PI", "PG", "PGal"]

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
                meta_info += "\t".join(sample_namex)

                # Finally, write the header part of the output VCF
                vcf_out.write(meta_info)
                vcf_out.write("\n")

                continue

            """' Now, extract the required data from each of the remaining lines add to output VCF. """
            updated_line = haplotype_to_vcf(
                line, header_line, all_alleles_idx, sample_namex, format_tagx, hap_format
            )
            vcf_out.write(updated_line)
            vcf_out.write("\n")

        print('Elapsed time : "%s".' % (time.time() - begin_time))


""" Part B: Function part of Haplotype To VCF """


def haplotype_to_vcf(line, header_line, all_alleles_idx, sample_namex, format_tagx, hap_format):

    line = line.rstrip("\n").split("\t")

    contig = line[0]
    pos = line[1]
    id_ = "."

    all_alleles = line[all_alleles_idx].split(",")
    ref = all_alleles[0]
    alt = ",".join(all_alleles[1:])
    qual = "."
    filter_ = "."
    info_ = "."
    format_ = ":".join(format_tagx)

    line_out = "\t".join([contig, pos, id_, ref, alt, qual, filter_, info_, format_])

    # Now, update SAMPLE:FORMAT values
    # Haplotype file exclusively will have PI and PG_al

    format_sample_values = []
    for namex in sample_namex:
        sample_PG_al = namex + ":PGal"
        sample_PG_al_idx = header_line.index(sample_PG_al)
        sample_PG_al_value = line[sample_PG_al_idx]

        sample_PI_idx = header_line.index(namex + ":PI")
        sample_PI_value = line[sample_PI_idx]

        # split the phased genotype bases (either iupac or numeric) 
        # and also keep the splitter/separator 
        pg_allele_split = re.split(r'([/|])', sample_PG_al_value)  
        
        # to store the updated alleles in their respective format 
        pg_bases_iupac = []; pg_bases_numeric = []
        
        if hap_format == 'numeric':
            # if haplotype in numeric format map the numeric alleles to all avaiable alleles             
            pg_bases_numeric = ''.join(pg_allele_split)            
            pg_bases_iupac = []
            for item in pg_allele_split:
                try:
                    pg_bases_iupac.append(all_alleles[int(item)])
                except ValueError:
                    pg_bases_iupac.append(item)
            pg_bases_iupac = ''.join(pg_bases_iupac)
                    
                    
        elif hap_format == 'iupac':
            # if haplotype is iupac map the iupac alleles to all alleles to find it's index position 
            pg_bases_iupac = ''.join(pg_allele_split)
            pg_bases_numeric = []
            
            for item in pg_allele_split:
                try:
                    pg_bases_numeric.append(str(all_alleles.index(item)))
                except ValueError:
                    pg_bases_numeric.append(str(item))
            pg_bases_numeric = ''.join(pg_bases_numeric)

        format_sample_values.append(
            ":".join(
                [pg_bases_numeric, sample_PI_value, pg_bases_numeric, pg_bases_iupac]
            )
        )

    line_out += "\t" + "\t".join(format_sample_values)

    return line_out
