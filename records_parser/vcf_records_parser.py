import shlex
import collections
from collections import OrderedDict
import sys
import re
import warnings


class _RecordParser:
    def __init__(self, line):
        self.lines = line

    def parse_record_keys(self):
        """ Parse the VCF header that starts with '#CHROM' to extract required field.
        e.g: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	\
                ms01e	ms02g	ms03g	ms04h	MA611	MA605	MA622 """
        # record_keys = self.lines.rstrip('\n').split('\t')
        # return record_keys
        pass

    def parse_record_values(self):
        """extract vcf record values matching each record keys"""
        record_values = self.lines.rstrip("\n").split("\t")
        return record_values

    @staticmethod
    def map_record_keys_to_values(record_keys, record_values):
        record_dict = OrderedDict(zip(record_keys, record_values))
        return record_dict


    @staticmethod
    def map_info_tags_to_values(string):
        # splitter = shlex.shlex(string, posix=True)
        # splitter.whitespace_split = True
        # splitter.whitespace = ";"
        # info_dict = OrderedDict(pair.split("=", 1) for pair in splitter)

        try:
            mapped_info = dict(s.split("=", 1) for s in string.split(";"))
        except ValueError:
            # warnings.warn(
            #     f"This info tag doesnot have '=' sign in info : {string}.\
            #     Such keys will be populated with '.' values"
            # )
            mapped_info = {}
            for s in string.split(";"):
                if "=" in s:
                    k, v = s.split("=")
                    mapped_info[k] = v
                else:
                    mapped_info[s] = "."
        return mapped_info

    @staticmethod    
    def convert_genotypes(ref_and_alt, numeric_genotype):
        iupac_genotype = []  # default value 
        
        numeric_genotype_split = re.split(r'([/|])', numeric_genotype)
        
        for its in numeric_genotype_split:
            try:
                if type(int(its)) is int:
                    iupac_genotype.append(ref_and_alt[int(its)])
            except ValueError:
                iupac_genotype.append(its)
                    
        iupac_genotype = ''.join(iupac_genotype)       
        
        return iupac_genotype
    


    @classmethod
    def map_format_tags_to_sample_values(cls, record_dict, format_tags, sample_names, gtbase_is):
        """map format tags to sample values in VCF file"""

        len_format_tags = len(format_tags)
        split_format_tags = format_tags.split(":")
        ref_and_alts = [record_dict['REF']] + record_dict['ALT'].split(',')

        for name in sample_names:
            sample_value_string = record_dict[name]
            sample_values = sample_value_string.split(":")

            # number of items in each sample should equal to number of format tags
            if len(sample_values) == 1:
                sample_values = sample_values + ["."] * (len_format_tags - 1)
            elif len(sample_values) == 0:
                sample_values = sample_values + ["."] * (len_format_tags)

            # map the format tags to the sample values
            mapped_format_sample = dict(zip(split_format_tags, sample_values))
            
            # update the genotype format for requested genotypes tags 
            #print(gtbase_is)
            for gts in gtbase_is: 
                gt_tag = gts[0]
                gt_output = gts[1]

                if  gt_output == 'iupac':
                    #numeric_gt = mapped_format_sample['GT']
                    numeric_gt = mapped_format_sample[gt_tag]
                    iupac_gt = cls.convert_genotypes(ref_and_alts, numeric_gt)
                    
                    # update the genotype values 
                    mapped_format_sample[gt_tag] = iupac_gt
                    #mapped_format_sample['GT'] = iupac_gt

            # update the sample record (tags, values) in the record_dict
            record_dict[name] = mapped_format_sample
            #print('record dict')
            #print(record_dict)
            #sys.exit(0)

        # after the for loop is complete all the the sample record should be updated
        # then return the record dict
        return record_dict
      


class VcfReadIndividualRecord:
    """Extract metadata from the header of the VCF file"""

    def __init__(self, recordKeys, inRecord, sample_names, gtbase_is):
        super().__init__()
        if not isinstance(inRecord, (str, bytes)):
            # raise string or byte not found error when no string is input
            raise TypeError("Records should be of type string or bytes")
        # self.inVCF = inVCF
        self.inRecord = inRecord
        self.recordKeys = recordKeys
        self.sample_names = sample_names
        self.gtbase_is = gtbase_is

        ## ?? - Bhuwan : not sure if we should keep/use these self variable below this line
        ## ?? - or may be this will be useful in the future, so we can return all the ..
        # .. mapped record as object which can be called using mapped_record.CHROM, mapped_record.REF
        self.record_vals = []
        self.record_dict = collections.OrderedDict()
        self.chrom_ = {}
        self.pos = {}
        self.id_ = {}
        self.ref = {}
        self.alt = {}
        self.qual = {}
        self.filter_flag = {}
        self.info_values = {}
        self.format_values = {}
        self.sample_values = {}

    def read_vcfRecord(self):
        inRecord = self.inRecord
        sample_names = self.sample_names

        if inRecord.startswith(r"#"):
            raise ValueError("Records should not start with #")

        else:
            # extract values for the vcf record
            self.record_vals = _RecordParser(inRecord).parse_record_values()

            """ map the record keys with the record values """
            self.record_dict = _RecordParser.map_record_keys_to_values(
                self.recordKeys, self.record_vals
            )

            """ map info tags to values and update the record_dict """
            info_tag_string = self.record_dict["INFO"]
            mapped_info_tags = _RecordParser.map_info_tags_to_values(info_tag_string)
            self.record_dict["INFO"] = mapped_info_tags

            """ map the FORMAT tags to the SAMPLE values of each sample """
            format_tags = self.record_dict["FORMAT"]
            self.record_dict = _RecordParser.map_format_tags_to_sample_values( 
                self.record_dict, format_tags, sample_names, self.gtbase_is
            )

            return self.record_dict
