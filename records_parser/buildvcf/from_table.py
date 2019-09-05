import time
import sys

from collections import OrderedDict
from metadata_parser.utils import time_memory_track

"""Step 03 (D): Function for Table To VCF"""
@time_memory_track
def fnc_table_to_vcf(infile, meta_header, outfile, samples, formats, infos, genotype_is):

    print("converting Table file to VCF")

    begin_time = time.time()

    gt_tag_as_iupac = []
    for gts_tag in genotype_is:
        tag_format = gts_tag.split(":")
        if tag_format[1] == "iupac":
            gt_tag_as_iupac.append(tag_format[0])

    with open(infile) as tablefile, open(meta_header) as meta_header, open(
        outfile, "w+"
    ) as vcf_out:

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

                contig_idx, contig_header = check_chrom_in_headerline(header_line)
                used_header.append(contig_header)

                pos_idx = check_pos_headerline(header_line)
                id_idx = header_line.index("ID") if "ID" in header_line else None
                ref_idx = header_line.index("REF") if "REF" in header_line else None
                alt_idx = header_line.index("ALT") if "ALT" in header_line else None
                qual_idx = header_line.index("QUAL") if "QUAL" in header_line else None
                filter_idx = (
                    header_line.index("FILTER") if "FILTER" in header_line else None
                )

                used_header.extend(
                    ["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                )

                """INFO tags are identified by matching "INFO:" in the column names."""
                infos_in_header = [x for x in header_line if x.startswith("INFO:")]
                all_infos = [x.replace("INFO:", "") for x in infos_in_header]

                info_tags = process_fields(
                    given_field=infos, all_fields=all_infos, tag="info"
                )

                # also find the position of the info tags on header line
                infos_idx = []
                if len(info_tags) != 0:
                    for inftag in info_tags:
                        infos_idx.append(header_line.index("INFO:" + inftag))
                else:
                    infos_idx = None

                """SAMPLE names and FORMAT tags are identified using ":" delimiter in the column names, 
                    after excluding the INFO fields."""
                possible_samples = [x for x in header_line if ":" in x]
                # get sample names and unique format tags  by removing info tags
                samples_and_formats = [
                    x for x in possible_samples if x not in infos_in_header
                ]

                # separate samples and formats
                all_samples = [x.split(":")[0] for x in samples_and_formats]
                all_formats = [x.split(":")[1] for x in samples_and_formats]

                # unique sample names and  format tags in order
                all_samples = set_of_list_order(all_samples)
                all_formats = set_of_list_order(all_formats)

                # find the available format tags
                format_tags = process_fields(
                    given_field=formats, all_fields=all_formats, tag="formats"
                )

                # In the available FORMAT tags, move "GT" field to the beginning.
                if "GT" in format_tags:
                    format_tags.remove("GT")
                    format_tags.insert(0, "GT")

                ### prepare sample names
                sample_names = process_fields(
                    given_field=samples, all_fields=all_samples, tag="samples"
                )
                used_header.extend(sample_names)
                print(used_header)

                """ Now, Read the meta header and add it to the output VCF file. """
                print('\nReading meta header from file "%s" ' % (meta_header.name))
                meta_info = get_meta_info(meta_header)

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
                # print(meta_info)

                # Finally, write the header part of the output VCF
                vcf_out.write(meta_info + "\n")

                continue

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

    chrom = line[contig_idx] if contig_idx is not None else "."
    pos = line[pos_idx] if pos_idx is not None else "."
    ids = line[id_idx] if id_idx is not None else "."
    ref = line[ref_idx] if ref_idx is not None else "."
    alt = line[alt_idx] if alt_idx is not None else "."
    qual = line[qual_idx] if qual_idx is not None else "."
    filter_ = line[filter_idx] if filter_idx is not None else "."
    format_ = ":".join(format_tags) if format_tags is not None else "."

    # Update "info tags and value". This is little complex
    if info_tags != []:
        info_ = []
        for ith, itemi in enumerate(info_tags):
            tag_val = "=".join([itemi, line[infos_idx[ith]]])
            info_.append(tag_val)
        info_ = ";".join(info_)
    elif info_tags == []:
        info_ = "."

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
                if sample_format_val in (".", "./.", ".|."):
                    continue
                else:
                    sep = "/" if "/" in sample_format_val else "|"
                    sample_format_val = sample_format_val.split(sep)
                    sample_format_val = [
                        all_alleles.index(sample_format_val[0]),
                        all_alleles.index(sample_format_val[1]),
                    ]

                    sample_format_val = sep.join(str(xth) for xth in sample_format_val)

            namex_vals.append(sample_format_val)

        format_sample_line.append(":".join(namex_vals))

    sample_format_final = "\t".join(format_sample_line)

    return sample_format_final


def check_chrom_in_headerline(header_line):
    if "#CHROM" in header_line:
        return header_line.index("#CHROM"), "#CHROM"
    elif "CHROM" in header_line:
        return header_line.index("CHROM"), "CHROM"
    else:
        print("CHROM field does not exist in the input table file. Update your file")
        print("Exiting the program")
        sys.exit(0)


def check_pos_headerline(header_line):
    if "POS" in header_line:
        return header_line.index("POS")
    else:
        print("POS field does not exist. Update your file")
        print("Exiting the program")
        sys.exit()


def process_fields(given_field, all_fields, tag):
    if given_field[0] == "all":
        return all_fields
    elif len(given_field) == 0:
        print(f"No {tag} available.")
        if tag == "info":
            print("INFO field will be populated with empty '.' value")
        return []
    else:
        if tag == "formats":
            return given_field
        else:
            nonsense_fields = [x for x in given_field if not x in all_fields]
            not_used_fields = [x for x in all_fields if not x in given_field]
            if tag == "info":
                if len(not_used_fields):
                    print(
                        "the following INFO tags won't be put in INFO fields of output VCF"
                    )

            if len(nonsense_fields):
                print(
                    f"The following {nonsense_fields} {tag} are not available in table file and not valid."
                )
                sys.exit(0)
            else:
                return given_field


def get_meta_info(meta_header):
    if meta_header:
        meta_info = meta_header.readlines()
        # if the meta header has "#CHROM	POS	REF ...." line then delete it
        if meta_info[-1].startswith("#CHROM\tPOS"):
            return "".join(meta_info[:-1]).rstrip("\n")
        else:
            return "".join(meta_info).rstrip("\n")

    else:
        print("Header with meta information is not provided")
        print("Exiting the program")
        sys.exit(0)


def set_of_list_order(input_list):
    outtemp = OrderedDict()
    for item in input_list:
        outtemp[item] = None
    return list(outtemp)
