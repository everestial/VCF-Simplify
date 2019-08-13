import time
import sys

from metadata_parser.vcf_metadata_parser import VcfReadMetaData
from metadata_parser.vcf_metadata_writer import Write_Simplified_MetaData
from metadata_parser.utils import vcf_records_as_table


from records_parser.simplifyvcf.to_haplotype import fnc_vcf_to_haplotype
from records_parser.simplifyvcf.to_table import fnc_vcf_to_table
from records_parser.buildvcf.from_table import fnc_table_to_vcf
from records_parser.buildvcf.from_haplotype import fnc_haplotype_to_vcf


def vcf_solver(task, args):
    
    print("Using the following arguments: ")
    print(args)
    print()    
    
    try:
        print('Using option "%s"' % task)
    except IndexError:
        print(
            "Provide only one of the following positional arguments: SimplifyVCF, BuildVCF or ViewVCF"
        )
        print("Exiting the program")
        print()
        sys.exit()

    ## starting VCF metadata parsing
    if task == "ViewVCF":
        start_time = time.time()
        print("  Extracting metadata from the VCF header ...")
        # extract all the available metadata
        metadata, raw_header, record_keys = VcfReadMetaData(args.inVCF).read_metadata()
        # but we will only use metadata dictionary

        if args.outFile:
            outFile = args.outFile

            if args.outType:
                data_types_requested = args.outType
                print("  Requested output data format %s" % (data_types_requested))
                # print(type(data_types_requested))

                # available output data types
                data_types_available = ["dict", "json", "table"]
                # print("\ndata types available\n", data_types_available)

                # match the requested data types with the available ones ..
                # and create truth table to write the output data in requested format
                as_dict, as_json, as_table = list(
                    map(lambda each: each in data_types_requested, data_types_available)
                )

                # pass the extracted metadata to the Writer and write it to a file
                print("  Writing metadata")
                Write_Simplified_MetaData(
                    metadata, outFile, as_dict, as_json, as_table
                ).write_or_print()

        if args.metadata:
            metadata_types_requested = args.metadata
            # Display the requested metadata information in table format
            print()
            print("Displaying requested metadata %s" % (metadata_types_requested))
            print()

            # convert metadata dict to table format and print in terminal
            vcf_records_as_table(metadata, metadata_types_requested)

        # calculate time
        end_time = time.time()
        run_time = end_time - start_time
        print(f"ViewVCF run time = {run_time : .4f} seconds.")

    ## Starting SimiplyVCF
    elif task == "SimplifyVCF":
        start_time = time.time()
        print("  Simplifying the VCF records ...")

        if args.toType[0] == "haplotype":
            # print(args)
            fnc_vcf_to_haplotype(args)

        elif args.toType[0] == "table":
            fnc_vcf_to_table(args)

        else:
            print("Incorrect file type. Enter valid type: table or haplotype")

    ## Starting BuildVCF
    elif task == "BuildVCF":
        start_time = time.time()
        print("  Building VCF ...")
        if args.fromType == "table":
            fnc_table_to_vcf(args)

        elif args.fromType == "haplotype":
            fnc_haplotype_to_vcf(args)

        else:
            print("Incorrect file type. Enter valid type: table or haplotype")

    else:
        print(
            "Provide one of the positional arguments: ViewVCF, SimplifyVCF, BuildVCF "
        )

        print("Exiting the program")
        print()
        sys.exit()
