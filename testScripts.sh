

## All the test below are some working examples with VCFsimplify 

#### A: ViewVCF
# write metadata as table, json, dict
python3 VcfSimplify.py ViewVCF -inVCF exampleInput/input_test.vcf -outFile exampleOutput/testOutput -outType table json dict

# to view only specific metadata of interest (on console)
python3 VcfSimplify.py ViewVCF -inVCF exampleInput/input_test.vcf -outFile exampleOutput/testOutput -metadata VCFspec samples

# to view only specific metadata of interest (write it to a text file)
python3 VcfSimplify.py ViewVCF -inVCF exampleInput/input_test.vcf -outFile exampleOutput/testOutput -metadata VCFspec samples > exampleOutput/VCFspecAndSample.txt


## B: SimplifyVCF

## B-01 : VCF to Table
# wide output, write GT bases as iupac, PG as iupac, write output header, 
python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -infos AF AN BaseQRankSum ClippingRankSum -formats GT PI PG GQ PL -preHeader CHROM POS REF ALT FILTER -mode wide -samples MA605 match:ms -GTbase GT:iupac PG:iupac -outHeaderName exampleOutput/vcf_header02.txt

# long output for the same arguments as above
python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -infos AF AN BaseQRankSum ClippingRankSum -formats GT PI PG GQ PL -preHeader CHROM POS REF ALT FILTER -mode long -samples MA605 match:ms -GTbase GT:iupac PG:iupac -outHeaderName exampleOutput/vcf_header02.txt

# table output without any INFO tags
python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -infos 0 -formats GT PI GQ PL -preHeader CHROM POS REF ALT FILTER -mode wide -samples MA605 ms01e  -GTbase GT:iupac

# no preheader, no infos
# write only sample and it's format tags output
python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -infos 0 -formats GT PI GQ PL -preHeader 0 -mode wide -samples MA605 ms01e

## a test on large VCF file (VCF to table) 
python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/SnpsMY-stitcher.vcf -outFile exampleOutput/simple_table.txt -infos AF AN BaseQRankSum ClippingRankSum -formats GT PI PG GQ PL -preHeader CHROM POS REF ALT FILTER -mode wide -samples match:A62 -GTbase GT:iupac



## B-02 : VCF to Haplotype 

# PI - as haplotype block index
# GT - as genotypes, as numeric, include unphased 
python3 VcfSimplify.py SimplifyVCF -toType haplotype -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_haplotype.txt -outHeaderName exampleOutput/vcf_header02.txt -PI PI -PG GT -GTbase GT:numeric -includeUnphased yes

# PI - as haplotype block index; 
# PG - as genotypes, as iupac, remove unphased 
python3 VcfSimplify.py SimplifyVCF -toType haplotype -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_haplotype.txt -outHeaderName exampleOutput/vcf_header02.txt -PI PI -PG PG -GTbase PG:iupac -includeUnphased no

# if you assume that PG are already phased genome wide you can use 'CHROM' as phased block 
# CHROM/contig - as haplotype block index; PG - as genotypes, as iupac, remove unphased 
python3 VcfSimplify.py SimplifyVCF -toType haplotype -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_haplotype.txt -outHeaderName exampleOutput/vcf_header02.txt -PI CHROM -PG PG -GTbase PG:iupac -includeUnphased no


##### C: BuildVCF 

### C-01 : from table to VCF 
## function 03 : table to VCF (only supported from wide table to VCF) 

# for the following VCF to table output
python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -samples ms01e MA605 -formats GT DP

# we can convert this above table back to VCF 
python3 VcfSimplify.py BuildVCF -fromType table -inFile exampleOutput/simple_table.txt -outVCF exampleOutput/tableToVcf.vcf -vcfHeader exampleOutput/vcf_header02.txt


## say we had created table with GT and PG in iupac bases like:
python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -samples ms01e MA605 -formats GT DP PI PG -GTbase GT:iupac PG:iupac

# we can convert back as 
python3 VcfSimplify.py BuildVCF -fromType table -inFile exampleOutput/simple_table.txt -outVCF exampleOutput/tableToVcf.vcf -vcfHeader exampleOutput/vcf_header02.txt -GTbase GT:iupac PG:iupac



## C-02 : from haplotype to VCF 

# for the following VCF to haplotype output
python3 VcfSimplify.py SimplifyVCF -toType haplotype -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_haplotype.txt -outHeaderName exampleOutput/vcf_header02.txt -GTbase PG:iupac -includeUnphased yes

# we can write back to VCF as 
python3 VcfSimplify.py BuildVCF -fromType haplotype -inFile exampleOutput/simple_haplotype.txt -outVCF exampleOutput/hapToVCF.vcf -vcfHeader exampleOutput/vcf_header.txt -haplotypeFormat iupac


# write to haplotype with no extra flags
# this will by default pick - PI = PI and PG = PG, 
# writes PG bases as numeric 
python3 VcfSimplify.py SimplifyVCF -toType haplotype -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_haplotype.txt -outHeaderName exampleOutput/vcf_header02.txt -GTbase PG:numeric

## convert the haplotype in numeric format back to VCF
python3 VcfSimplify.py BuildVCF -fromType haplotype -inFile exampleOutput/simple_haplotype.txt -outVCF exampleOutput/hapToVCF.vcf -vcfHeader exampleOutput/vcf_header.txt -haplotypeFormat numeric






