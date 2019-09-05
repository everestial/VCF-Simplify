import os
import filecmp
import pytest


from records_parser.simplifyvcf.to_haplotype import fnc_vcf_to_haplotype
from records_parser.simplifyvcf.to_table import fnc_vcf_to_table


input_vcf = 'tests/input_test.vcf'
meta_header = 'tests/vcf_header.txt'
actual_outfile1 = 'tests/hapfromvcf'
header_name = 'tests/vcf_header_from_hap.txt'
def test_vcf_to_hap(tmpdir):
    temp_outfile = tmpdir/ "hap1.txt"
    temp_header = tmpdir/ "temp_header.txt"
    fnc_vcf_to_haplotype(input_vcf, temp_outfile, temp_header, pi_tag= 'PI', pg_tag = 'PG', include_unphased= 'yes', gtbase= ['PG:numeric'])
    
    assert is_same_file(actual_outfile1, temp_outfile) and is_same_file(header_name, temp_header) is True


actual_outfile2 = 'tests/tablefromvcf.table'
header_name = 'tests/vcf_header_from_table.txt'
pre_header = ['all']
samples, formats , infos = ['all'], ['all'], ['all']
mode = 'wide'
gtbase = ['PG:numeric']

def test_vcf_to_table(tmpdir):
    temp_outfile = tmpdir/ "table.txt"
    temp_header = tmpdir/ "temp_header2.txt"
    fnc_vcf_to_table(input_vcf, temp_outfile, pre_header, mode, gtbase, temp_header, infos, formats, samples)
    
    assert is_same_file(actual_outfile2, temp_outfile) and is_same_file(header_name, temp_header) is True



def is_same_file(file1, file2):
    return filecmp.cmp(file1, file2)