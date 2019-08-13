"""Collection of several utilities used while parsing VCF."""


def vcf_records_as_table(metadata, metadata_types_requested):
    """Print the requested metadata to the console. """
    for keys in metadata_types_requested:
        for (ks_parent, vs_parent) in metadata.items():  # parent keys and parent values
            if ks_parent == keys:
                print(r"##" + ks_parent)
                vs_keys_to_write = [kys for kys, vys in vs_parent[0].items()]
                print("#" + "\t".join(vs_keys_to_write))

                for nested_vs in vs_parent:
                    vs_vs_to_write = []
                    for _, vs_child in nested_vs.items():
                        vs_vs_to_write.append(str(vs_child))
                    print("\t".join(vs_vs_to_write))
                print()

    # print("completed")
