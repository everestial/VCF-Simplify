#!/home/bin/python3

import argparse

from cyvcf2 import VCF
from cyvcf2 import VCFReader

import collections
import time

print()
print('Checking required modules')
print()

print('''Purpose of the program: simplify the data representation of the VCF file.
The output file consists of data (from INFO and FORMAT field) representing SAMPLE of interest.
Data are simplified for each requested sample.''')
print()


def main():

    # Step 01: define argument variables
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="sorted vcf file as input", required=True)
    parser.add_argument("--out", help="name of the output file that contains simplified VCF as table", required=True)
    parser.add_argument("--samples", default='all',
                        help="SAMPLE that are of interest; write as comma separated names\n"
                             "for e.g: 'sampleA,sampleB' or 'all' ")
    parser.add_argument("--pre_header", default='all',
                        help="comma separated pre-header fields; write as comma separated fields; "
                             "for e.g: 'CHR,POS,ID' or 'all' ")
    parser.add_argument("--infos", default='all',
                        help="INFO fields that are of interest; write as comma separated tags; "
                             "for e.g: 'AC,AF,AN' or 'all' ")
    parser.add_argument("--formats", default='all',
                        help="FORMAT fields that are of interest; for e.g: 'GT,PG,PI' or 'all' ")
    parser.add_argument("--keep_header", default='no',
                        help="keep the HEADER data in the output file.\n"
                             "options: 'yes' or 'no' ")


    global args  # creating a global argument variable
    args = parser.parse_args()


    # ********************  only activate during non-interactive mode
    # Step 02: Set the parameters that are of interest in VCF file
    #pre_header = 'CHR,POS'
    #pre_header = 'all'

    #info_of_interest = 'AC,AN'
    #info_of_interest = 'all'

    #sample_of_interest = 'MA611,ms02g'
    #sample_of_interest = 'all'

    #format_of_iterest = 'GT,PG,PL'
    #format_of_iterest = 'all'

    # keep_header = 'yes'  # if 'yes' then add the header to the output file, but set default at 'no'
    # keep_header = 'no'
    # ************************


    # Step 02: Now, pipe the "input arguments" to a variable
    pre_header = args.pre_header
    info_of_interest = args.infos
    sample_of_interest = args.samples
    format_of_iterest = args.formats
    keep_header = args.keep_header


    # Step 03: Read vcf file using cyvcf2 and start mining the data
    start_time01 = time.time()
    #with open("simplified_vcf.txt", 'w') as write_block:
    with open(args.out, 'w') as write_block:
        #vcf_file = VCF('F1.phased_updated_variants.Final02.vcf')
        vcf_file = VCF(args.vcf)
        sample_ids = vcf_file.samples
        print(sample_ids)
        print('%i samples found' %len(sample_ids))
        print()

        # mining header
        # add argument to keep or discard header while writing output file
        header = vcf_file.raw_header.split('\n')

        if keep_header == 'yes':
            write_block.write(vcf_file.raw_header)
        print()

        # Step 03-A: now, write the appropriate front part of the header
        if pre_header != 'all':
            output_header = pre_header.split(',')

        else:
            output_header = ['CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']



        # Step 03-B: mine fields that are of interest - for rear part of the header and for computing values from variants.
        # all data are returned as list from the defined function
        # the output data are: sample names, infos and formats of interest which can be used in downstream analyses
        my_samples, my_infos, my_formats = process_fields_of_interest(
            header, sample_of_interest, info_of_interest , format_of_iterest, sample_ids)


        # create empty dictionary - to store Key, Values from FORMAT field for each sample
        format_ids_values = collections.OrderedDict()
        for ks in my_formats:
            format_ids_values[ks] = None

        my_sample_idx = [sample_ids.index(x) for x in my_samples]  # we can move this elsewhere

        # Step 03-C: add "INFO" of interest to the header
        output_header += my_infos

        # Step 03-D: Now, simplify the TAGS from the FORMAT field by assigning it to each SAMPLE names
        # using data from "my_samples", "my_formats" extend the "header" of the output file.
        for name in my_samples:
            for tags in my_formats:
                output_header.append(name + '_' + tags)

        # write the final header to the output file
        print('\t'.join(output_header), file=write_block)

        print('output header of the simplified VCF file')
        print(output_header)
        print()


        chr_on_process = ''

        ''' Step 04: now, start parsing the VCF file using cyVCF2 and add the data for each header fields'''
        print('reading the input vcf file ')
        print()
        for variant in vcf_file:
            contig = str(variant.CHROM)
            # find which chr is in the process
            if chr_on_process != contig:
                print('contig %s is being processed' % str(contig))
                print()
                chr_on_process = contig

            # **************************
            # these methods are deprecated for now, but keeping for future use.
            #pos = str(variant.POS)
            #id_ = variant.ID
            #ref_allele = variant.REF
            #alt_alleles = variant.ALT
            #all_alleles = [ref_allele] + alt_alleles
            #alt_freq = variant.INFO.get('AF')

            # pass the "alt_freq" values to a function to compute "all_freq" as string
            #all_freq = compute_allele_freq(alt_freq)
            # ***********************************


            # Step 04-A : Write the values for the pre-fields of the header (i.e CHR upto FILTER)
            # If user desires less number of fields that is also achievable
            chr_to_filter = str(variant).split('\t')[0:7]

            pre_header_dict = collections.OrderedDict()  # to store the data when not all pre-headers are of interest
            if pre_header == 'all':
                write_block.write('\t'.join(chr_to_filter))

            elif pre_header != 'all':
                for idx, heads in enumerate(pre_header):
                    pre_header_dict[heads] = chr_to_filter[idx]

                write_block.write('\t'.join(pre_header_dict.values()))


            # Step 04-B: compute values for the INFO tags of interest
            infos_to_write = process_info(my_infos, variant)
            write_block.write('\t' + infos_to_write)


            # Step 04-C: compute values for the FORMAT fields of interest for each SAMPLE names of interest
            # so, we need to use both format_fields and sample_names
            # and pass it to a defined function
            process_format(variant, my_sample_idx, format_ids_values, write_block)


            write_block.write('\n')

        print('elapsed time: ', time.time() - start_time01)
        print()


# mine fields-TAGS that are of interest (INFO, FORMAT and SAMPLE)
def process_fields_of_interest(header, sample_names, info_tags, format_tags, all_samples):
    print('extracting the following fields of interest are :')
    print()
    infos_of_interest = []
    formats_of_interest = []
    samples_of_interest = []

    # part 01: mine SAMPLE names of interest
    if sample_names == 'all':
        samples_of_interest = all_samples
    elif sample_names != 'all':
        samples_of_interest = sample_names.split(',')

    print('SAMPLE of interest: ')
    print(samples_of_interest)
    print()

    # part 02: mine INFO tags of interest
    if info_tags == 'all':
        for line in header:
            if line.startswith('##INFO'):
                infos_of_interest += [line.split(',')[0].split('=')[2]]
    elif info_tags != 'all':
        infos_of_interest = info_tags.split(',')

    print('INFO of interest: ')
    print(infos_of_interest)
    print()

    # part 03: mine FORMAT fields of interest
    if format_tags == 'all':
        for line in header:
            if line.startswith('##FORMAT'):
                formats_of_interest += [line.split(',')[0].split('=')[2]]
    elif format_tags != 'all':
        formats_of_interest = format_tags.split(',')

    print('FORMAT of interest: ')
    print(formats_of_interest)
    print()

    return samples_of_interest, infos_of_interest, formats_of_interest


# function to compute the values from INFO's field of interest
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

    return infos_to_write



# now, for each SAMPLE compute the FORMAT's field of interest
def process_format(variant, my_sample_idx, format_ids_values, write_block):

    format_ids_in_variant = str(variant).split('\t')[8].split(':')
    sample_data_in_variant = str(variant).rstrip('\n').split('\t')[9::]


    for nth in my_sample_idx:
        nth_sample = sample_data_in_variant[nth].split(':')

        # update the values for each keys in the dictionary
        format_ids_vals = list(zip(format_ids_in_variant, nth_sample))

        for ks, vs in format_ids_vals:
            if ks in format_ids_values.keys():
                format_ids_values[ks] = vs    # condition: if ks in format_ids_values.keys()


        # Now, write the values to each FORMAT field, for each sample
        for vs in format_ids_values.values():
            write_block.write('\t' + str(vs))


### ************  these two functions are deprecated for now
# this function is deprecated for now - will need some update in the future
def process_add_gt_bases():
    # this is returned as array. for now it is just easy to mine data from the "all_alleles" variable
    # ** so this may be used in the future
    gt_bas = variant.gt_bases

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

    # for the situation when ref and/or alt freq are 'NoneType' and exception arises
    elif alt_freq is None:
        alt_freq = '.'
        ref_freq = '.'
        all_freq = [ref_freq, alt_freq]

    # conver the float instances to string in "freq" data
    all_freq = [str(x) for x in all_freq]
    return all_freq
# ****************************


if __name__ == '__main__':
    main()









