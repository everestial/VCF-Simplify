import pytest
import filecmp


from records_parser.buildvcf.from_haplotype import fnc_haplotype_to_vcf
from records_parser.buildvcf.from_table import fnc_table_to_vcf


input_hap = 'tests/hapfromvcf'
meta_header = 'tests/vcf_header.txt'
actual_outfile1 = 'tests/vcf_from_hap.vcf'

def test_hap_to_vcf(tmpdir):
    temp_outfile = tmpdir/ "vcf1.txt"
    fnc_haplotype_to_vcf(infile= input_hap, meta_header = meta_header, outfile=temp_outfile , hap_format= 'iupac')
    assert is_same_file(actual_outfile1, temp_outfile) is True




# Test from table to vcf
intable = 'tests/tablefromvcf.table'
meta_header = 'tests/vcf_header.txt'
actual_outfile2 = 'tests/vcf_from_table.vcf'
samples= ['all']
formats = ['all']
infos = ['all']
genotype_is = ['PG:numeric']
def test_table_to_vcf(tmpdir):
    temp_outfile = tmpdir/ "table1.txt"
    fnc_table_to_vcf(intable, meta_header, temp_outfile, samples, formats, infos, genotype_is)   
    assert is_same_file(actual_outfile2, temp_outfile) is True
    

def is_same_file(file1, file2):
    return filecmp.cmp(file1, file2)
