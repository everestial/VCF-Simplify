#!/home/bin/python3

import argparse

from cyvcf2 import VCF
from cyvcf2 import VCFReader

import numpy as np

import collections
import time




print()
print('- Checking required modules...')
print()

print('''- Purpose of the program:
    Simplify the data representation of a VCF file.
    The output file consists of data from INFO and FORMAT field representing SAMPLE of interest.''')
print()


def main():

    # Step 01: define argument variables
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="Sorted VCF file as input", required=True)
    parser.add_argument("--out", help="Name of the output file that contains simplified VCF as table.",
                        required=True)
    parser.add_argument("--samples",
                        help="SAMPLE of interest; write as comma separated names, "
                             "for e.g: 'sampleA,sampleB' or 'all'.", default='all')
    parser.add_argument("--pre_header",
                        help="Comma separated pre-header fields before the 'INFO' field in the input VCF file. "
                             "Write as comma separated fields, for e.g: 'CHR,POS,ID' or 'all'. "
                             "Default: 'all'. ",
                        default='all')
    parser.add_argument("--infos",
                        help="INFO tags that are of interest; write as comma separated tags; "
                             "for e.g: 'AC,AF,AN' or 'all'.",
                        default='all')
    parser.add_argument("--formats",
                        help="FORMAT tags that are of interest; for e.g: 'GT,PG,PI' or 'all'.",
                        default='all')
    parser.add_argument("--keep_header", default='no',
                        help="Keep the HEADER data in the output file."
                             "Options: 'yes' or 'no' ")

    parser.add_argument("--mode", help="Structure of the output table."
                                       "Options: wide(0), long(1). Default: 0 .",
                        required=False, default=0)
    parser.add_argument("--gtbase", help="write the GT field as IUPAC base code."
                                       "Options: no(0), yes(1). Default: 0 .",
                        required=False, default=0)




    global args  # creating a global argument variable
    args = parser.parse_args()

    global gtbase


    # ********************  only activate during non-interactive mode
    # Step 02: Set the parameters that are of interest in VCF file
    #pre_header = 'CHR,POS'
    #pre_header = 'all'

    #info_of_interest = 'AC,AN'
    #info_of_interest = 'all'

    #sample_of_interest = 'MA611,ms02g'
    #sample_of_interest = 'all'

    #format_of_interest = 'GT,PG,PL'
    #format_of_interest = 'all'

    # keep_header = 'yes'  # if 'yes' then add the header to the output file, but set default at 'no'
    # keep_header = 'no'
    # ************************


    # Step 02: Now, pipe the "input arguments" to a variable
    pre_header = args.pre_header
    info_of_interest = args.infos
    sample_of_interest = args.samples
    format_of_interest = args.formats
    keep_header = args.keep_header

    if args.mode == '1' or args.mode == 'long':
        mode = 'long'
    else: mode = 'wide'

    if args.gtbase == '1' or args.gtbase == 'yes':
        gtbase = 'yes'
    else: gtbase = 'no'


    # Step 03: Read vcf file using cyvcf2 and start mining the data
    start_time01 = time.time()
    #with open("simplified_vcf.txt", 'w') as write_block:
    with open(args.out, 'w') as write_block:
        #vcf_file = VCF('input_test.vcf')
        vcf_file = VCF(args.vcf)
        sample_ids = vcf_file.samples
        #print(sample_ids)
        print('- %i samples found.' %len(sample_ids))
        print()

        # mining header
        # add argument to keep or discard header while writing output file
        header = vcf_file.raw_header.split('\n')

        if keep_header == 'yes':
            write_block.write(vcf_file.raw_header)
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
            header, sample_of_interest, info_of_interest , format_of_interest, sample_ids)


        # create empty dictionary - to store Key, Values from FORMAT field of each sample
        format_ids_values = collections.OrderedDict()
        for ks in my_formats:
            format_ids_values[ks] = None

        my_sample_idx = [sample_ids.index(x) for x in my_samples]  # we can move this elsewhere


        ## Step 03-C: add tags from the "INFO" of interest to the header
        output_header += my_infos


        ## Step 03-D: Decide between "long" vs. "wide" way of representing FORMAT field values
        # Now, simplify the TAGS from the FORMAT field by assigning it to each SAMPLE names
        # add the tags from the "FORMAT" field as it is. Each sample is represented in different line
        if mode == 'long':
            output_header.append('SAMPLE')
            output_header.extend(my_formats)

            # write the final header to the output file
            print('\t'.join(output_header), file=write_block)


        elif mode == 'wide':
            # add the tags from the "FORMAT" field as suffix to the sample of interest (on output header)
            for name in my_samples:
                for tags in my_formats:
                    output_header.append(name + '_' + tags)

            # write the final header to the output file
            print('\t'.join(output_header), file=write_block)

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
            #pos = str(variant.POS)
            #id_ = variant.ID
            #ref_allele = variant.REF
            #alt_alleles = variant.ALT
            #all_alleles = [ref_allele] + alt_alleles
            #alt_freq = variant.INFO.get('AF')

            # pass the "alt_freq" values to a function to compute "all_freq" as string
            #all_freq = compute_allele_freq(alt_freq)
            ##########################  **************************


            # Step 04-A : Write the values for the pre-fields of the header (i.e CHR upto FILTER)
            # If user desires less number of fields that is also achievable
            chr_to_filter = str(variant).split('\t')[0:7]

            # to store the data when not all pre-headers are of interest
            pre_header_dict = collections.OrderedDict()
            if pre_header == 'all':
                #write_block.write('\t'.join(chr_to_filter))
                line_to_write += '\t'.join(chr_to_filter)

            elif pre_header != 'all':
                for idx, heads in enumerate(pre_header.split(',')):
                    pre_header_dict[heads] = chr_to_filter[all_header.index(heads)]

                #write_block.write('\t'.join(pre_header_dict.values()))
                line_to_write += '\t'.join(pre_header_dict.values())


            # Step 04-B: compute values for the INFO tags of interest
            infos_to_write = process_info(my_infos, variant)
            #write_block.write('\t' + infos_to_write)
            line_to_write += '\t' + infos_to_write


            # Step 04-C: compute values for the FORMAT fields of interest for each SAMPLE names of interest
            # so, we need to use both format_fields and sample_names together
            # and pass it to a defined function
            if mode == 'wide':
                process_format_wide(variant, my_sample_idx, format_ids_values, write_block, line_to_write)

            elif mode == 'long':
                process_format_long(variant, my_sample_idx, format_ids_values, write_block, line_to_write, sample_ids)


            #write_block.write('\n')

        print('Elapsed time: ', time.time() - start_time01)
        print()


# mine fields-TAGS that are of interest (INFO, FORMAT and SAMPLE)
def process_fields_of_interest(header, sample_of_interest, info_of_interest, format_of_interest, all_samples):

    print('- Extracting the following fields of interest:')
    print()
    info_tags = []
    format_tags = []
    sample_names = []

    # part 01: mine SAMPLE names of interest
    if sample_of_interest == 'all':
        sample_names = all_samples
    elif sample_of_interest != 'all':
        sample_names = sample_of_interest.split(',')

    print('\tSAMPLE of interest are: ')
    print('\t', sample_names)
    print()

    # part 02: mine INFO tags of interest
    if info_of_interest == 'all':
        for line in header:
            if line.startswith('##INFO'):
                info_tags += [line.split(',')[0].split('=')[2]]
    elif info_of_interest  != 'all':
        info_tags = info_of_interest.split(',')

    print('\tINFO of interest are: ')
    print('\t', info_tags)
    print()

    # part 03: mine FORMAT fields of interest
    if format_of_interest == 'all':
        for line in header:
            if line.startswith('##FORMAT'):
                format_tags += [line.split(',')[0].split('=')[2]]
    elif format_of_interest != 'all':
        format_tags = format_of_interest.split(',')

    print('\tFORMAT of interest are: ')
    print('\t', format_tags)
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
                format_ids_values[ks] = vs    # condition: if ks in format_ids_values.keys()

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
                format_ids_values[ki] = vi    # condition: if ks in format_ids_values.keys()

        # Now, write the values to each FORMAT field, for each sample
        for ks, vs in format_ids_values.items():
            if gtbase == 'yes':
                if ks == 'GT':
                    vs = variant.gt_bases[nth]
            write_block.write('\t' + str(vs))

        write_block.write('\n')




if __name__ == '__main__':
    main()









