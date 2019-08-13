import shlex
import collections
import re


## ?? rename this to '_MetadataParser"
class _MetadataParser:
    def __init__(self, line, tag):
        self.lines = line
        self.tag = tag

    @staticmethod
    def split_to_dict(string):
        splitter = shlex.shlex(string, posix=True)
        splitter.whitespace_split = True
        splitter.whitespace = ","
        tags_dict = dict(pair.split("=", 1) for pair in splitter)
        return tags_dict

    def parse_gvcf_block(self):
        """extract the GVCF blocks"""
        # e.g: ##GVCFBlock55-56=minGQ=55(inclusive),maxGQ=56(exclusive)
        gvcf_block = re.search("##GVCFBlock(.*?)=", self.lines, 0).group(1)

        # update tags_dict
        to_replace = self.tag + gvcf_block + "="

        string = self.lines.rstrip("\n").rstrip(">").replace(to_replace, "")
        tags_dict = self.split_to_dict(string)
        tags_dict["Block"] = gvcf_block
        return tags_dict

    def parse_contigs(self):
        """find chromosomes/contig names in the input VCF file"""
        # e.g line: ##contig=<ID=scaffold_1118,length=1005>

        string = self.lines.rstrip("\n").rstrip(">").replace(self.tag, "")
        tags_dict = self.split_to_dict(string)
        return tags_dict

    def parse_reference_genome(self):
        """find the reference genome name used to generated the input VCF/GVCF"""
        # return reference genome name and file path
        # e.g: ##reference=file:///media/02_Alignment_To_Ref/GVCF_Calls_onGenes-RunningTest01/02-JointGenotyping_MySpF1/lyrata_genome.fa
        ref_genome = self.lines.rstrip("\n").replace(self.tag, "")
        return {"reference": ref_genome}

    def parse_gatk_commands(self):
        """find the GATK commands used to generate the input VCF"""
        ## e.g: ##GATKCommandLine.HaplotypeCaller=<ID=HaplotypeCaller....
        gatk_cmd_middle_string = re.search(
            "##GATKCommandLine(.*)=<ID=", self.lines
        ).group(1)
        to_replace = "##GATKCommandLine" + gatk_cmd_middle_string + "=<"
        string = self.lines.rstrip("\n").rstrip(">").replace(to_replace, "")
        tags_dict = self.split_to_dict(string)
        return tags_dict

    @staticmethod  ## unused right now. Use in the future.
    def write_header(self):
        """write the exact metadata header to a separate file"""
        return None

    def parse_sample_names(self):
        """
        Parse the sample names
        - sample names begins at 10th position of the #CHROM header line.
        """
        string = self.lines.rstrip("\n").split("\t")
        sample_names = string[9::]  # extract all the sample names
        sample_pos = [
            string.index(x) + 1 for x in sample_names
        ]  # find sample positions

        # create a list of dictionary of sample {names:position}
        pos_dict = [
            {"name": x, "position": y} for x, y in zip(sample_names, sample_pos)
        ]

        return pos_dict

    def parse_format_info_filter(self):
        """
        Parse the FORMAT, INFO, FILTER from the VCF file. 
        These tags contains same string structure so a single function suffices for all.
        """
        string = self.lines.rstrip("\n").rstrip(">").replace(self.tag, "")
        tags_dict = self.split_to_dict(string)
        return tags_dict

    def parse_vcf_spec(self):
        """identify the VCF format version"""
        string = self.lines.rstrip("\n").replace(r"##", "")
        tags_dict = self.split_to_dict(string)
        return tags_dict

    def parse_record_keys(self):
        """ Parse the VCF header that starts with '#CHROM' to extract required field.
        e.g: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	\
                ms01e	ms02g	ms03g	ms04h	MA611	MA605	MA622 """
        record_keys = self.lines.rstrip("\n").lstrip(r"#").split("\t")
        return record_keys


class VcfReadMetaData:
    """Extract metadata from the header of the VCF file"""

    def __init__(self, inVCF):
        super().__init__()
        if not inVCF:
            raise FileNotFoundError
        self.inVCF = inVCF
        self.is_gvcf = False
        self.metadict = collections.OrderedDict()
        self.vcf_spec = []
        self.infos_ = []
        self.formats_ = []
        self.filters_ = []
        self.contig = []
        self.samples = []
        self.reference_genome = []
        self.gatk_commands = []
        self.gvcf_blocks = []
        self.raw_header = ""
        self.record_keys = []

    def read_metadata(self):
        with open(self.inVCF) as infile:
            for lines in infile:
                """ This can go into the method called capture header"""
                if lines.startswith("#"):
                    self.raw_header += lines

                    ## find VCF version and if the VCF is gVCF
                    if lines.startswith(r"##fileformat"):
                        _vcf_spec = _MetadataParser(
                            lines, tag=r"##fileformat"
                        ).parse_vcf_spec()
                        self.vcf_spec.append(_vcf_spec)

                    elif lines.startswith(r"##GVCFBlock"):
                        self.is_gvcf = True
                        self.vcf_spec[0]["GVCF"] = self.is_gvcf  # returns GVCF is TRUE

                        # parse the GVCF blocks data
                        _gvcf_block = _MetadataParser(
                            lines, tag=r"##GVCFBlock"
                        ).parse_gvcf_block()
                        self.gvcf_blocks.append(_gvcf_block)

                    ## extract existing FORMAT, INFO and FILTER tags
                    elif lines.startswith(r"##FORMAT"):
                        self.formats_.append(
                            _MetadataParser(
                                lines, tag=r"##FORMAT=<"
                            ).parse_format_info_filter()
                        )
                    elif lines.startswith(r"##INFO"):
                        self.infos_.append(
                            _MetadataParser(
                                lines, tag=r"##INFO=<"
                            ).parse_format_info_filter()
                        )
                    elif lines.startswith(r"##FILTER"):
                        self.filters_.append(
                            _MetadataParser(
                                lines, tag=r"##FILTER=<"
                            ).parse_format_info_filter()
                        )

                    ## extract contig/chromosomes included in the VCF file
                    elif lines.startswith(r"##contig"):
                        self.contig.append(
                            _MetadataParser(lines, tag=r"##contig=<").parse_contigs()
                        )

                    ## extract the name/path of the reference genome
                    elif lines.startswith(r"##reference"):
                        self.reference_genome.append(
                            _MetadataParser(
                                lines, tag=r"##reference="
                            ).parse_reference_genome()
                        )
                        ##reference=

                    ## extract the GATK commands used so far in the preparation of the VCF file
                    elif lines.startswith(r"##GATKCommandLine"):
                        self.gatk_commands.append(
                            _MetadataParser(
                                lines, tag=r"##GATKCommandLine"
                            ).parse_gatk_commands()
                        )

                    ## extract sample names
                    elif lines.startswith("#CHROM"):
                        # self.parse_sample_names(lines)
                        self.samples = _MetadataParser(
                            lines, tag=r"#CHROM"
                        ).parse_sample_names()

                        self.record_keys = _MetadataParser(
                            lines, tag=""
                        ).parse_record_keys()

                        break  # the for loop

                    ## ignore now - also start mining the Variant information from each line of the VCF
                    # else:

        ## store the records mmeta-data in dictionary
        self.metadict["VCFspec"] = self.vcf_spec
        self.metadict["FORMAT"] = self.formats_
        self.metadict["INFO"] = self.infos_
        self.metadict["FILTER"] = self.filters_
        self.metadict["contig"] = self.contig
        self.metadict["reference"] = self.reference_genome
        self.metadict["GATKCommandLine"] = self.gatk_commands
        self.metadict["GVCFBlock"] = self.gvcf_blocks
        self.metadict["samples"] = self.samples

        return self.metadict, self.raw_header, self.record_keys
