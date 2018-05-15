#!/home/bin/python3

import argparse
from cyvcf2 import VCF
import numpy as np
import time
import sys
import collections


print()
print('Checking required modules')
print()

''' Purpose of the program: mine the data from the VCF files and convert it into Haplotype file.
The output file consists of ReadBackPhased genotypes as PI (block keys) and PG (genotype values). Other GT with no phased
state may be extracted as well.'''


def main():

    # define argument variables
    main_parser = argparse.ArgumentParser(prog="VCF-Simplify")

    # Create sub_parsers (one for SimplifyVCF, other for BuildVCF)
    subparsers = main_parser.add_subparsers(help='Choose one of the following method.')


    '''Part A: Sub parser for SimplifyVCF - to create simplified table from VCF file. '''
    parser_a = subparsers.add_parser('SimplifyVCF',
                                     help='Simplify VCF : to a haplotype or a table file.')

    # upperlevel parser within SimplifyVCF
    parser_a.add_argument("-toType", help="Type of the output file. Option: haplotype, table", required=True)
    parser_a.add_argument("-inVCF", help="sorted vcf file", required=True)
    parser_a.add_argument("-out", help="name of the output file", required=True)
    parser_a.add_argument("-keepHeader", default='no',
                        help="Write the HEADER data to a separate output file."
                             "Options: 'yes' or 'no' ")

    '''Part A (01) : From VCF to Haplotype '''
    vcf_to_haplotype = parser_a.add_argument_group('Flags for "VCF To Haplotype"')

    vcf_to_haplotype.add_argument("-PG", required=False, default='PG',
                        help="FORMAT tag containing the phased genotype of the SAMPLE. "
                             "Only applicable on 'haplotype file output'. ")
    vcf_to_haplotype.add_argument("-PI", required=False, default='PI',
                        help="FORMAT tag representing the unique index of RBphased haplotype block in the SAMPLE. "
                             "Only applicable on 'haplotype file output'. "
                             "Note: 'CHROM' can also be used as PI if VCF is phased chromosome-wide. ")
    vcf_to_haplotype.add_argument("-unphased", default='no', required=False,
                        help="include unphased variants in the output. "
                             "Aavailable options: yes, no")


    '''Part A (02) : From VCF to Table'''
    vcf_to_table = parser_a.add_argument_group('Flags for "VCF To Table"')

    vcf_to_table.add_argument("-samples",
                        help="SAMPLE of interest; write as comma separated names, "
                             "for e.g: 'sampleA,sampleB' or 'all'.", default='all')
    vcf_to_table.add_argument("-preHeader",
                        help="Comma separated pre-header fields before the 'INFO' field in the input VCF file. "
                             "Write as comma separated fields, for e.g: 'CHR,POS,ID' or 'all'. "
                             "Default: 'all'. ",
                        default='all')
    vcf_to_table.add_argument("-infos",
                        help="INFO tags that are of interest; write as comma separated tags; "
                             "for e.g: 'AC,AF,AN' or 'all'.",
                        default='all')
    vcf_to_table.add_argument("-formats",
                        help="FORMAT tags that are of interest; for e.g: 'GT,PG,PI' or 'all'.",
                        default='all')
    vcf_to_table.add_argument("-mode", help="Structure of the output table."
                                       "Options: wide(0), long(1). Default: 0 .",
                        required=False, default=0)
    vcf_to_table.add_argument("-gtbase", help="write the GT field as IUPAC base code."
                                       "Options: no(0), yes(1). Default: 0 .",
                        required=False, default=0)



    '''Part B: Sub parser for "BuildVCF". To create VCF from simple Table like format. '''
    parser_b = subparsers.add_parser('BuildVCF',
                                     help='Create VCF : from a haplotype or a table file. ')

    # upper level parser within BuildVCF

    parser_b.add_argument("-fromType", required=True,
                        help="Type of the input file the VCF is being prepared from. "
                                      "Options: haplotype, table ")
    parser_b.add_argument("-inFile", required=True,
                        help="Sorted table or haplotype file."
                             "This haplotype file should be obtained from phase-Stitcher, "
                             "phase-Extender. The table file should be in the format "
                             "output by 'VCF-Simplify'; only long format table is supported for now.")

    parser_b.add_argument("-outVCF", help="Name of the output VCF file.", required=True)
    parser_b.add_argument("-vcfHeader", required=True,
                        help="A custom VCF header to add to the VCF file. "
                             "The VCF header should not contain the line with #CHROM .... ")


    '''Part B (01) : Only From Table To VCF. '''
    """ Additional argument parser only to use if "-fromType" is "table" """
    table_to_vcf_parser = parser_b.add_argument_group('Flags for "Table To VCF"')

    table_to_vcf_parser.add_argument("-GTbase", help="Representation of the GT base is : numeric, IUPAC ",
                                     required=False)

    table_to_vcf_parser.add_argument("-samples",
                                     help="Name of the samples -> "
                                          "comma separated name of the samples that needs to be converted "
                                          "to VCF format",
                                     default='all', required=False)

    table_to_vcf_parser.add_argument("-formats",
                                     help="Name of the FORMAT tags to write -> "
                                          "comma separated FORMAT tags name.",
                                     default='all', required=False)

    table_to_vcf_parser.add_argument("-infos",
                                     help="Name of the INFO tags to write -> "
                                          "comma separated INFO tags name. ",
                                     default='all', required=False)



    global args  # creating a global argument variable
    args = main_parser.parse_args()



    """ Step 02: Based on positional arugments and task go to specific function """

    try :
        print('Using option "%s"' %sys.argv[1])
    except IndexError:
        print('Provide one of the positional arguments: SimplifyVCF or BuildVCF ')
        print()
        sys.exit()

    if sys.argv[1] == 'SimplifyVCF':
        if args.toType == 'haplotype':
            '''Go to " 03 (A) : VCF to Haplotype or Table" function if following conditions are true. '''
            print('Converting VCF To Haplotype.')
            fnc_vcf_to_haplotype()

        elif args.toType == 'table':
            '''Go to " 03 (B) : VCF to Table" function if following conditions are true. '''
            print('Converting VCF To Table')
            fnc_vcf_to_table()
        else:
            print('toType is not indicated.')


    elif sys.argv[1] == 'BuildVCF':

        if args.fromType == 'haplotype':
            '''Go from " 03 (C) Haplotype to VCF" if following conditions are true. '''
            print('Converting Haplotype file to VCF')
            fnc_haplotype_to_vcf()

        elif args.fromType == 'table':
            '''Go from " 03 (D) Table to VCF" if following conditions are true. '''
            print('Converting Table to VCF')
            fnc_table_to_vcf()

        else:
            print('fromType is not indicated.')

    else:
        print('Provide one of the positional arguments: SimplifyVCF or BuildVCF ')
        print()
        sys.exit()







'''Step 03 (A): Function for VCF To Haplotype '''
def fnc_vcf_to_haplotype():

    print('Creating Haplotype file from VCF file')
    start_time01 = time.time()
    pi_tag = args.PI
    pg_tag = args.PG

    with open(args.out, 'w') as write_block:
        vcf_file = VCF(args.inVCF)
        sample_ids = vcf_file.samples

        print('%i samples found' %len(sample_ids))
        print()


        # mining header
        # add argument to keep or discard header while writing output file
        header = vcf_file.raw_header.rstrip('\n').split('\n')

        if args.keepHeader == 'yes':
            print("Writing the header to a separate output file.")
            with open('vcf_header.txt', 'w') as vcfheader:
                vcfheader.write(vcf_file.raw_header)
        else:
            print("Skipping the header.")
        print()



        output_header = ['CHROM', 'POS', 'all-alleles']
        for name in sample_ids:
            output_header.append(name + ':PI')  # to store "phased index"
            output_header.append(name + ':PG_al')  # to store the "phased genotype" as IUPAC bases

        # write the header of the output file
        #print('\t'.join(output_header), file=write_block)
        write_block.write('\t'.join(output_header) + '\n')


        chr_on_process = ''

        ''' now, start parsing the VCF file using pyVCF '''
        print('reading the input vcf file %s' % str(args.inVCF))
        print()
        for variant in vcf_file:
            contig = str(variant.CHROM)
            pos = str(variant.POS)
            ref_ = variant.REF
            alt_ = variant.ALT
            all_alleles = [ref_] + alt_


            # find which chr is in the process
            if chr_on_process != contig:
                print('contig %s is being processed' %str(contig))
                print()
                chr_on_process = contig


            # write data to the output file
            # ('2', 15881018, 'G', ['G', 'A', 'C'], [0.0, 1.0])
            write_block.write('\t'.join([contig, pos, ','.join(all_alleles)]))


            # this is returned as array. for now it is just easy to mine data from the "all_alleles" variable
            # ** so this may be used in the future
            gt_bas = variant.gt_bases

            # start mining the phased_index (PI) and phased_genotype (PG) values for each sample
            # this outputs a list of PI and PG for all available sample in the VCF file


            ## ** update this part with args.PI, args.PG
            if pi_tag == 'CHROM':
                pi_values = np.repeat(contig, len(sample_ids), axis=0)
            else:
                pi_values = variant.format(pi_tag)

            if pg_tag == 'GT':
                gts = variant.genotypes
                
                # use class (Genotype()) to further update the pg_values. 
                # ** prolly deprecate in future if cyvcf2 updates
                pg_values = [Genotype(li) for li in gts]
            else:
                pg_values = variant.format(pg_tag)

            # incase the "PI" and "PG" tags are missing in the VCF file
            # seems like "bcftools" merged variants are missing some of the PI and PG tags.
            if pi_values is None:
                pi_values = np.repeat('.', len(sample_ids), axis=0)

            if pg_values is None:
                pg_values = np.repeat('.', len(sample_ids), axis=0)


            ''' now, start extracting sample level information '''
            for ith, pi_value in enumerate(pi_values):
                pg_val = str(pg_values[ith])

                if pg_val == '.':
                    pg_allele = '.'

                elif '/' in pg_val:
                    # if unphased variants are not of interest we return them as '.'
                    if args.unphased == 'no':
                        pg_allele = '.'

                    else:
                        if pg_val == './.':
                            pg_allele = '.'

                        else:
                            pg_val = pg_val.split('/')
                            pg_allele = all_alleles[int(pg_val[0])] + '/' + all_alleles[int(pg_val[1])]

                elif '|' in pg_val:
                    if pg_val == '.|.':
                        pg_allele = '.'
                    else:
                        pg_val = pg_val.split('|')
                        pg_allele = all_alleles[int(pg_val[0])] + '|' + all_alleles[int(pg_val[1])]


                ### Deprecated feature.
                # ** for the future: in some instances there are "*" in the GT, which we could replace as "."
                # this may change in the future
                #if '*' in pg_allele:
                    #pg_allele = '.'

                write_block.write('\t' + '\t'.join([pi_value, pg_allele]))
            write_block.write('\n')

        print('elapsed time: ', time.time() - start_time01)
        print()



'''Class to extract the genotype (numeric bases) when "PG" is set as "GT" '''
class Genotype(object):
    __slots__ = ('alleles', 'phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)
    __repr__ = __str__



'''Step 03 (B) : Function for VCF To Table.'''
def fnc_vcf_to_table():
    global gtbase

    # Step 02: Now, pipe the "input arguments" to a variable
    pre_header = args.preHeader
    info_of_interest = args.infos
    sample_of_interest = args.samples
    format_of_interest = args.formats
    keep_header = args.keepHeader

    if args.mode == '1' or args.mode == 'long':
        mode = 'long'
    else:
        mode = 'wide'

    if args.gtbase == '1' or args.gtbase == 'yes':
        gtbase = 'yes'
    else:
        gtbase = 'no'

    # Step 03: Read vcf file using cyvcf2 and start mining the data
    start_time01 = time.time()
    # with open("simplified_vcf.txt", 'w') as write_block:
    with open(args.out, 'w') as write_block:
        vcf_file = VCF(args.inVCF)
        sample_ids = vcf_file.samples
        print(sample_ids)
        print('- %i samples found.' % len(sample_ids))
        print()

        # mining header
        # add argument to keep or discard header while writing output file
        header = vcf_file.raw_header.rstrip('\n').split('\n')

        if keep_header == 'yes':
            print("Writing the header to a separate output file.")
            with open('vcf_header.txt', 'w') as vcfheader:
                vcfheader.write(vcf_file.raw_header)
        else:
            print("Skipping the header.")
        print()

        # Step 03-A: now, write the appropriate front part of the header
        all_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
        if pre_header != 'all':
            output_header = pre_header.split(',')
        else:
            output_header = all_header

        # Step 03-B: mine fields that are of interest - for rear part of the header and for computing values from variants.
        # all data are returned as list from the defined function
        # the output data are: sample names, infos and formats of interest which can be used in downstream analyses
        my_samples, my_infos, my_formats = process_fields_of_interest(
            header, sample_of_interest, info_of_interest, format_of_interest, sample_ids)

        # create empty dictionary - to store Key, Values from FORMAT field of each sample
        format_ids_values = collections.OrderedDict()
        for ks in my_formats:
            format_ids_values[ks] = None

        my_sample_idx = [sample_ids.index(x) for x in my_samples]  # we can move this elsewhere

        ## Step 03-C: first add tags from the "INFO" of interest to the header
        output_header += my_infos

        ## Step 03-D: Decide between "long" vs. "wide" way of representing FORMAT field values
        # Now, simplify the TAGS from the FORMAT field by assigning it to each SAMPLE names
        # add the tags from the "FORMAT" field as it is. Each sample is represented in different line
        if mode == 'long':
            output_header.append('SAMPLE')
            output_header.extend(my_formats)

            # write the final header to the output file
            #print('\t'.join(output_header), file=write_block)
            write_block.write('\t'.join(output_header) + '\n')


        elif mode == 'wide':
            # add the tags from the "FORMAT" field as suffix to the sample of interest (on output header)
            for name in my_samples:
                for tags in my_formats:
                    output_header.append(name + ':' + tags)

            # write the final header to the output file
            #print('\t'.join(output_header), file=write_block)
            write_block.write('\t'.join(output_header) + '\n')

        print()
        print('- Output header of the simplified VCF file:')
        print(output_header)
        print()

        chr_on_process = ''

        ''' Step 04: now, start parsing the VCF file using cyVCF2 and add the data for each header fields'''
        print('Reading the input vcf file ... ')
        print()
        for variant in vcf_file:

            line_to_write = ''  # create new emtpy variable

            contig = str(variant.CHROM)
            # find which chr is in the process
            if chr_on_process != contig:
                print('Contig %s is being processed ... ' % str(contig))
                print()
                chr_on_process = contig

            ##########################  **************************
            # these methods are deprecated for now, but keeping for future use.
            # pos = str(variant.POS)
            # id_ = variant.ID
            # ref_allele = variant.REF
            # alt_alleles = variant.ALT
            # all_alleles = [ref_allele] + alt_alleles
            # alt_freq = variant.INFO.get('AF')

            # pass the "alt_freq" values to a function to compute "all_freq" as string
            # all_freq = compute_allele_freq(alt_freq)
            ##########################  **************************

            # Step 04-A : Write the values for the pre-fields of the header (i.e CHR upto FILTER)
            # If user desires less number of fields that is also achievable
            chr_upto_filter = str(variant).split('\t')[0:7]

            # to store the data when not all pre-headers are of interest
            pre_header_dict = collections.OrderedDict()
            if pre_header == 'all':
                # write_block.write('\t'.join(chr_upto_filter))
                line_to_write += '\t'.join(chr_upto_filter)

            elif pre_header != 'all':
                for idx, heads in enumerate(pre_header.split(',')):
                    pre_header_dict[heads] = chr_upto_filter[all_header.index(heads)]

                # write_block.write('\t'.join(pre_header_dict.values()))
                line_to_write += '\t'.join(pre_header_dict.values())

            # Step 04-B: compute values for the INFO tags of interest
            infos_to_write = process_info(my_infos, variant)
            # write_block.write('\t' + infos_to_write)
            line_to_write += '\t' + infos_to_write

            # Step 04-C: compute values for the FORMAT fields of interest for each SAMPLE names of interest
            # so, we need to use both format_fields and sample_names together
            # and pass it to a defined function
            if mode == 'wide':
                process_format_wide(variant, my_sample_idx, format_ids_values, write_block, line_to_write)

            elif mode == 'long':
                process_format_long(variant, my_sample_idx, format_ids_values, write_block, line_to_write, sample_ids)

            # write_block.write('\n')

        print('Elapsed time: ', time.time() - start_time01)
        print()

    # mine fields-TAGS that are of interest (INFO, FORMAT and SAMPLE)


def process_fields_of_interest(header, sample_of_interest, info_of_interest, format_of_interest, all_samples):
    print('## Extracting the following fields of interest:')
    print()
    info_tags = []
    format_tags = []
    sample_names = []

    # part 01: mine SAMPLE names of interest
    if sample_of_interest == 'all':
        sample_names = all_samples
    elif sample_of_interest != 'all':
        sample_names = sample_of_interest.split(',')

    print('  # SAMPLE of interest are: ')
    print('  ', sample_names)
    print()

    # part 02: mine INFO tags of interest
    if info_of_interest == 'all':
        for line in header:
            if line.startswith('##INFO'):
                info_tags += [line.split(',')[0].split('=')[2]]
    elif info_of_interest != 'all':
        info_tags = info_of_interest.split(',')

    print('  # INFO of interest are: ')
    print('  ', info_tags)
    print()

    # part 03: mine FORMAT fields of interest
    if format_of_interest == 'all':
        for line in header:
            if line.startswith('##FORMAT'):
                format_tags += [line.split(',')[0].split('=')[2]]
    elif format_of_interest != 'all':
        format_tags = format_of_interest.split(',')

    print('  # FORMAT of interest are: ')
    print('  ', format_tags)
    print()

    return sample_names, info_tags, format_tags

    # function to compute values from INFO's field of interest


def process_info(info_fields, variant):
    infos_to_write = []
    for field in info_fields:
        field_value = variant.INFO.get(field)

        if isinstance(field_value, float):
            field_value = round(field_value, 3)

        if isinstance(field_value, tuple):
            field_value = tuple([round(x, 3) for x in field_value])

        infos_to_write += [field_value]

    infos_to_write = ('\t'.join(str(x) for x in infos_to_write)).replace('(', '').replace(')', '').replace(' ', '')
    # **improve above code in the future

    return infos_to_write

    # now, for each SAMPLE compute and write the FORMAT's field of interest


def process_format_wide(variant, my_sample_idx, format_ids_values, write_block, line_to_write):
    # write the data before the "FORMAT" field begins
    write_block.write(line_to_write)

    # all available FORMAT tags
    format_ids_in_variant = str(variant).split('\t')[8].split(':')

    # all available SAMPLE data
    sample_data_in_variant = str(variant).rstrip('\n').split('\t')[9::]

    # So, now we mine values for requested FORMAT tags for requested SAMPLE
    for nth in my_sample_idx:
        # all the tag values in "nth" sample
        nth_sample = sample_data_in_variant[nth].split(':')

        # create a ZIP of (tags,values) and store as list of tuples
        format_ids_vals = list(zip(format_ids_in_variant, nth_sample))

        for ks, vs in format_ids_vals:
            if ks in format_ids_values.keys():
                format_ids_values[ks] = vs  # condition: if ks in format_ids_values.keys()

        # Now, write the values to each FORMAT field, for each sample
        for ks, vs in format_ids_values.items():
            if gtbase == 'yes':
                if ks == 'GT':
                    vs = variant.gt_bases[nth]
            write_block.write('\t' + str(vs))

    write_block.write('\n')


def process_format_long(variant, my_sample_idx, format_ids_values, write_block, line_to_write, sample_ids):
    # all available FORMAT tags
    format_ids_in_variant = str(variant).split('\t')[8].split(':')

    # all available SAMPLE data
    sample_data_in_variant = str(variant).rstrip('\n').split('\t')[9::]

    # So, now we mine values for requested FORMAT tags for requested SAMPLE
    for nth in my_sample_idx:
        write_block.write(line_to_write + '\t' + sample_ids[nth])

        # all the tag values in "nth" sample
        nth_sample = sample_data_in_variant[nth].split(':')

        # create a ZIP of (tags,values) and store as list of tuples
        format_ids_vals = list(zip(format_ids_in_variant, nth_sample))

        for ki, vi in format_ids_vals:
            if ki in format_ids_values.keys():
                format_ids_values[ki] = vi  # condition: if ks in format_ids_values.keys()

        # Now, write the values to each FORMAT field, for each sample
        for ks, vs in format_ids_values.items():
            if gtbase == 'yes':
                if ks == 'GT':
                    vs = variant.gt_bases[nth]
            write_block.write('\t' + str(vs))

        write_block.write('\n')



'''Step 03 (C): Function for Haplotype To VCF'''
def fnc_haplotype_to_vcf():
    print('converting Haplotype file to VCF')
    begin_time = time.time()

    '''Assign some input variables. '''
    infile = args.inFile
    meta_header = args.vcfHeader
    outfile = args.outVCF

    with open(infile) as hapfile, \
            open(meta_header) as meta_header, \
            open(outfile, 'w+') as vcf_out:

        '''Start reading the haplotype file as generator. This saves memory. '''
        for line in hapfile:
            if line.startswith('CHROM') \
                    or line.startswith('#CHROM'):
                header_line = line.rstrip('\n').split('\t')

                if 'CHROM' in header_line:
                    contig_idx = header_line.index('CHROM')
                elif '#CHROM' in header_line:
                    contig_idx = header_line.index('#CHROM')
                else:
                    print('CHROM field does not exit. Update your file')
                    break

                if 'POS' in header_line:
                    pos_idx = header_line.index('POS')
                else:
                    print('POS field does not exit. Update your file')
                    break

                if 'all-alleles' in header_line:
                    all_alleles_idx = header_line.index('all-alleles')
                else:
                    print('"all-alleles" field not available in input file. Update your file')
                    break

                '''Finally find available SAMPLE names and it's FORMAT tags'''
                sample_namex = []
                for itemx in header_line:
                    if ':' in itemx:
                        sample_namex.append(itemx.split(':')[0])
                sample_namex = list(set(sample_namex))

                # assign FORMAT tags - keeping it fixed
                format_tagx = ['GT', 'PI', 'PG', 'PG_al']

                ''' Now, Read the meta header and add it to the output VCF file. '''
                print('\nReading meta header from file "%s" ' % (meta_header.name))
                if meta_header != None:
                    meta_info = meta_header.read().rstrip('\n')
                    meta_info += '\n'
                else:
                    print('Header with meta information is not provided')
                    sys.exit()

                # add meta header to the output VCF file
                meta_info += '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                                        'QUAL', 'FILTER', 'INFO', 'FORMAT']) + '\t'

                # add SAMPLE fields to output VCF file
                meta_info += '\t'.join(sample_namex)

                # Finally, write the header part of the output VCF
                vcf_out.write(meta_info)
                vcf_out.write('\n')

                continue

            '''' Now, extract the required data from each of the remaining lines add to output VCF. '''
            updated_line = haplotype_to_vcf(
                line, header_line, all_alleles_idx, sample_namex, format_tagx)
            vcf_out.write(updated_line)
            vcf_out.write('\n')

        print('Elapsed time : "%s".' % (time.time() - begin_time))



''' Part B: Function part of Haplotype To VCF '''

def haplotype_to_vcf(line, header_line, all_alleles_idx, sample_namex, format_tagx):

    line = line.rstrip('\n').split('\t')

    contig = line[0]
    pos = line[1]
    id = '.'

    all_alleles = line[all_alleles_idx].split(',')
    ref = all_alleles[0]
    alt = ','.join(all_alleles[1:])
    qual = '.'
    filter_ = '.'
    info_ = '.'
    format_ = ':'.join(format_tagx)

    line_out = '\t'.join([contig, pos, id, ref, alt, qual, filter_, info_, format_])

    # Now, update SAMPLE:FORMAT values
    # Haplotype file exclusively will have PI and PG_al

    format_sample_values = []
    for namex in sample_namex:
        sample_PG_al = namex + ':PG_al'
        sample_PG_al_idx = header_line.index(sample_PG_al)
        sample_PG_al_value = line[sample_PG_al_idx]

        sample_PI_idx = header_line.index(namex + ':PI')
        sample_PI_value = line[sample_PI_idx]

        # to store the values for GT and PG tags
        sample_GT_value = '.'
        sample_PG_value = '.'

        if sample_PG_al_value == '.' \
                or sample_PG_al_value == './.' \
                or sample_PG_al_value == '.|.':
            sample_GT_value = '.'

        elif '/' in sample_PG_al_value:
            sample_GT_value = sample_PG_al_value.split('/')
            sample_GT_value = [all_alleles.index(sample_GT_value[0]),
                               all_alleles.index(sample_GT_value[1])]
            sample_GT_value = '/'.join(str(sth) for sth in sample_GT_value)
            sample_PG_value = sample_GT_value

        elif '|' in sample_PG_al_value:
            sample_GT_value = sample_PG_al_value.split('|')
            sample_GT_value = [all_alleles.index(sample_GT_value[0]),
                               all_alleles.index(sample_GT_value[1])]
            sample_GT_value = '|'.join(str(sth) for sth in sample_GT_value)
            sample_PG_value = sample_GT_value

        format_sample_values.append(
            ':'.join([sample_GT_value, sample_PI_value, sample_PG_value, sample_PG_al_value]))

    line_out += '\t' + '\t'.join(format_sample_values)

    return line_out



##########################################

'''Step 03 (D): Function for Table To VCF'''
def fnc_table_to_vcf():
    print('converting Table file to VCF')

    ## declare globals
    global genotype_is
    global begin_time
    global contig_idx
    global pos_idx
    global id_idx
    global ref_idx
    global alt_idx
    global qual_idx
    global filter_idx

    # INFO, FORMAT and SAMPLE don't have index but tags
    global info_tags
    global infos_idx
    global format_tags
    global reduced_format_tags
    global sample_names
    global header_line

    begin_time = time.time()

    '''Assign some input variables. '''
    genotype_is = args.GTbase
    infile = args.inFile
    meta_header = args.vcfHeader
    outfile = args.outVCF
    samples = args.samples
    formats = args.formats
    infos = args.infos


    with open(infile) as hapfile, \
        open(meta_header) as meta_header,\
            open(outfile, 'w+') as vcf_out:

        '''Start reading the haplotype file as generator. This saves memory. '''
        for line in hapfile:

            ## find and set the indexes ...
            # ... of pre-fields, INFO, FORMAT and SAMPLE level information

            ''' Step 01: The very first line of the file is read;
             - to find the variable name and it's index position in the input file.
             - almost all the variable created downstream are "global variables".
             - SAMPLE level information is automatically identified unless explicitly given.
               The sample names is identified using ":" in the column names, so other names
               should not have ":" at all.
             - FORMAT level tags can also be provided as list, or can be mined automatically
               along with SAMPLE by using ":" matching.
             - All the preHeader tags, ie. CHROM  POS  ID  REF  ALT  QUAL  FILTER are reserved and
               updated by matching the names in text header line.
             '''

            # to use the "header" name that have already been taken
              # this will help in finding appropriate "INFO" level tags from the header file
            used_header = []

            if line.startswith('CHROM') \
                    or line.startswith('#CHROM'):
                header_line = line.rstrip('\n').split('\t')

                if 'CHROM' in header_line:
                    contig_idx = header_line.index('CHROM')
                elif '#CHROM' in header_line:
                    contig_idx = header_line.index('#CHROM')
                else:
                    print('CHROM field does not exit. Update your file')
                    break

                # update the taken header "labels"
                used_header += ['CHROM', '#CHROM']

                if 'POS' in header_line :
                    pos_idx = header_line .index('POS')
                else:
                    print('POS field does not exit. Update your file')
                    break

                # update
                used_header += ['POS']

                if 'ID' in header_line :
                    id_idx = header_line .index('ID')
                else:
                    id_idx = None
                used_header += ['ID']

                if 'REF' in header_line :
                    ref_idx = header_line .index('REF')
                else:
                    ref_idx == None
                used_header += ['REF']

                if 'ALT' in header_line :
                    alt_idx = header_line .index('ALT')
                else:
                    alt_idx = None
                used_header += ['ALT']

                if 'QUAL' in header_line :
                    qual_idx = header_line .index('QUAL')
                else:
                    qual_idx = None
                used_header += ['QUAL']

                if 'FILTER' in header_line :
                    filter_idx = header_line .index('FILTER')
                else:
                    filter_idx = None
                used_header += ['FILTER']


                '''SAMPLE names and FORMAT tags are identified using ":" delimiter in the column names. '''
                if samples != 'all':
                    sample_names = samples.split(',')
                elif samples == 'all':
                    sample_names = []
                    for itemy in header_line:
                        if ':' in itemy:
                            sample_names.append(itemy.split(':')[0])
                    sample_names = list(set(sample_names))

                used_header += [x for x in header_line if ':' in x]

                # find the available format tags
                if formats != 'all':
                    format_tags = formats.split(',')

                elif formats == 'all':
                    format_tags = []
                    for names in sample_names:
                        for itemx in header_line :
                            if itemx.startswith(names):
                                format_tags.append(itemx.split(':')[1])

                    format_tags = list(set(format_tags))


                # In the available FORMAT tags, move "GT" field to the beginning.
                if 'GT' in format_tags:
                    format_tags.remove('GT')
                    format_tags.insert(0, 'GT')



                ''' Finally, update the tag names of the "INFO" field '''
                #** Note: Any column names in the header line that is not taken so far is ..
                  # .. considered a "INFO" field.
                remaining_cols = [itx for itx in header_line if itx not in set(used_header)]

                if infos != 'all':
                    info_tags = infos.split(',')
                elif infos == 'all' and len(remaining_cols) > 0:
                    info_tags = remaining_cols
                else:
                    info_tags = None
                    print('INFO tags are not available.')

                # also find the position of the info tags on header line
                infos_idx = []
                if info_tags != None:
                    for inftag in info_tags:
                        infos_idx.append(header_line.index(inftag))
                else:
                    infos_idx = None



                ''' Now, Read the meta header and add it to the output VCF file. '''
                print('\nReading meta header from file "%s" ' %(meta_header.name))

                if meta_header != None:
                    meta_info = meta_header.read()
                else:
                    print('Header with meta information is not provided')
                    break

                # add meta header to the output VCF file
                meta_info += '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                                         'QUAL', 'FILTER', 'INFO', 'FORMAT']) + '\t'

                # add SAMPLE fields to output VCF file
                meta_info += '\t'.join(sample_names)


                # Finally, write the header part of the output VCF
                vcf_out.write(meta_info + '\n')

                continue
                #break


            '''' Now, extract the required data from each of the remaining lines add to output VCF. '''
            updated_line = table_to_vcf(line)
            vcf_out.write(updated_line)
            vcf_out.write('\n')

        print('Elapsed time : "%s".' %(time.time()-begin_time))



'''Function part of Table to VCF '''
def table_to_vcf(line_in):

    line = line_in.rstrip('\n').split('\t')
    chrom = line[contig_idx]
    pos = line[pos_idx]

    if id_idx is not None:
        ids = line[id_idx]
    else:
        ids = '.'

    ref = line[ref_idx]
    alt = line[alt_idx]

    if qual_idx is not None:
        qual = line[qual_idx]
    else: qual = '.'

    if filter_idx is not None:
        filter = line[filter_idx]
    else: filter = '.'

    # Update "info tags and value". This is little complex
    if info_tags !=None:
        info_ = []
        for ith, itemi in enumerate(info_tags):
            tag_val = '='.join([itemi, line[infos_idx[ith]]])
            info_.append(tag_val)
        info_ = ';'.join(info_)
    elif info_tags is None:
        info_ = '.'

    # write the tags names of the FORMAT column
    if format_tags != None:
        format_ = ':'.join(format_tags)
    else:format_ = '.'


    # update the output line
    line_out = '\t'.join([chrom, pos, ids, ref, alt, qual, filter, info_, format_]) + '\t'

    # Further update the SAMPLE-to-FORMAT values
    # pass the line to another function
    format_to_sample_vals = update_sample_format(line, ref, alt)
    line_out = line_out + format_to_sample_vals

    return line_out

''' Function part of Table to VCF '''
def update_sample_format(line, ref, alt):

    # The "line" variable is passed into this function.
    # The global variables are "genotype_is", "sample_names" and "format_tags"

    # to store updated line
    format_sample_line = []

    all_alleles = [ref] + alt.split(',')

    for namex in sample_names:
        namex_vals = []
        for tagx in format_tags:
            sample_format_tag = namex + ':' + tagx
            sample_format_idx = header_line.index(sample_format_tag)
            sample_format_val = line[sample_format_idx]

            ''' further update the sample:format value if GT in table is as IUPAC base '''
            if tagx == 'GT' and genotype_is == 'IUPAC':
                if sample_format_val == '.' or \
                        sample_format_val == './.' or \
                        sample_format_val == '.|.':
                    continue

                elif '/' in sample_format_val:
                    sample_format_val = sample_format_val.split('/')

                    sample_format_val = [all_alleles.index(sample_format_val[0]),
                                         all_alleles.index(sample_format_val[1])]

                    sample_format_val = '/'.join(str(xth) for xth in sample_format_val)

                elif '|' in sample_format_val:
                    sample_format_val = sample_format_val.split('|')

                    sample_format_val = [all_alleles.index(sample_format_val[0]),
                                         all_alleles.index(sample_format_val[1])]

                    sample_format_val = '|'.join(str(xth) for xth in sample_format_val)

            namex_vals.append(sample_format_val)

        format_sample_line.append(':'.join(namex_vals))

    sample_format_final = '\t'.join(format_sample_line)

    return sample_format_final




####################################









### Deprecated function (for VCF to Haplotype file). Probably useful in future.
''' function to compute allele frequencies at each position of the VCF file'''
def compute_allele_freq(alt_freq):
    # alt_freq is returned as either float or tuple of floats; so we extract data in more universal format
    if isinstance(alt_freq, tuple):
        alt_freq = [round(x, 3) for x in alt_freq]
        ref_freq = round(1 - sum(alt_freq), 3)
        all_freq = [ref_freq] + alt_freq

    elif isinstance(alt_freq, float):
        alt_freq = round(alt_freq, 3)
        ref_freq = round(1 - alt_freq, 3)
        all_freq = [ref_freq, alt_freq]


if __name__ == '__main__':
    main()

    