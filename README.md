# VCF-Simplify &ensp;&ensp;&ensp; v2.2

A python parser to simplify the vcf file into table like format.
There are several tools available to mainpulate and alter VCF file.
But, a simple and comprehensive tool that can produce a most simple
output required by emperical biologist is still amiss.

- [VCF-Simplify &ensp;&ensp;&ensp; v2.1](#vcf-simplify-enspenspensp-v21)
  - [Tutorial](#tutorial)
    - [Prerequisites](#prerequisites)
    - [Installation and setup](#installation-and-setup)
  - [Usage](#usage)
    - [ViewVCF](#viewvcf)
    - [SimplifyVCF](#simplifyvcf)
    - [BuildVCF](#buildvcf)
  - [Upcoming features:](#upcoming-features)
    - [Citation:](#citation)

This tool performs following three tasks:

1. ViewVCF (Display and extract the metadata from the vcf files.)
2. SimplifyVCF (Conversion of vcf files into table format and haplotype.)
3. BuildVCF (Conversion from table or haplotyoe to vcf files.)

- **Convert VCF to TABLE**\
    This tool takes in sorted vcf file and reports a simplified table output
for `INFO` and `FORMAT` field for each `SAMPLE` of interest. With default
state (minimal code) all the `INFO`, `FORMAT` for all the `SAMPLE` are
simplified. Fields can be further narrowed down using very convenient
and comprehensive scripts.

- **Convert TABLE to VCF**\
    It is also possible to convert the TABLE file into VCF. Controlled,
    workflows are included.

**Exclusively for [phase-Stitcher](https://github.com/everestial/pHASE-Stitcher) and [phase-Extender](https://github.com/everestial/phase-Extender).**\

  - **Convert VCF to Haplotype**\
    It is also possible to convert the TABLE file into VCF. Controlled,
    workflows are included.

  - **Convert Haplotype to VCF**\
    It is also possible to convert the TABLE file into VCF. Controlled,
    workflows are included.

## Tutorial

### Prerequisites

**VCF Simplify** is written in python3, so you need to have python3 installed on your system to run this code locally. If you don't have python installed then, you can install from [here](https://www.python.org/downloads/). For linux; you can get latest python3 by:

` sudo apt-get install python3`


### Installation  and setup

1. Clone this repo.

``` 
git clone https://github.com/everestial/VCF-SimplifyDev
cd VCF-SimplifyDev
```

2. Make virtual env for python and install requirements.

```
python3 -m venv myenv
source myenv/bin/activate   # for linux
myenv\Scripts\activate      # for windows
pip install Cython
```

3. In order to cythonize(so that app runs faster):
   
`python3 setup.py build_ext --inplace`

This makes .so files in your build directory and it will be called when importing modules rather than actual modules.

4. Run help on VCFSimplify by:
   
<pre>
$ python3 VcfSimplify.py -h

    ## VCF Simplify ## : Python application for parsing VCF files.
    Author: Bishwa K. Giri
    Contact: bkgiri@uncg.edu, kirannbishwa01@gmail.com

usage: VCF-Simplify [-h] {ViewVCF,SimplifyVCF,BuildVCF} ...

positional arguments:
  {ViewVCF,SimplifyVCF,BuildVCF}
                        Choose one of the following method.
    ViewVCF               View and extract metadata from the VCF.
    SimplifyVCF           Simplify VCF -> to {haplotype, table} data.
    BuildVCF              Build VCF <- from {haplotype, table} data.

optional arguments:
  -h, --help            show this help message and exit

</pre>

## Usage

### ViewVCF

- Help on view vcf:


<pre>
$ python3 VcfSimplify.py ViewVCF -h

    ## VCF Simplify ## : Python application for parsing VCF files.
    Author: Bishwa K. Giri
    Contact: bkgiri@uncg.edu, kirannbishwa01@gmail.com

usage: VCF-Simplify ViewVCF [-h] -inVCF INVCF [-outFile OUTFILE]
                            [-outType {table,json,dict} [{table,json,dict} ...]]
                            [-metadata METADATA [METADATA ...]]

optional arguments:
  -h, --help            show this help message and exit
                        
  -inVCF INVCF          Sorted vcf file.
                        
  -outFile OUTFILE      Name of the output file without file extension.
                        
  -outType {table,json,dict} [{table,json,dict} ...]
                        Space separated list of output data types.
                        Multiple types can be requested.
                        
  -metadata METADATA [METADATA ...]
                        Space separated list of metadata of interest.
                        Allowed values are: 
                          VCFspec, reference, contig, samples, INFO, FORMAT, FILTER, GATKCommandLine, GVCFBlock.
                        Multiple choices can be requested.
</pre>

- Example A-01 Write metadata from VCF file as table, json, dict
  ```python3 VcfSimplify.py ViewVCF -inVCF exampleInput/input_test.vcf -outFile exampleOutput/testOutput -outType table json dict```
  
  This will create output of vcf headers into three format(table, json and dict) as testOutput.table and so on.

- Example A-02 View only specific metadata of interest (on console)

  ```python3 VcfSimplify.py ViewVCF -inVCF exampleInput/input_test.vcf -outFile exampleOutput/testOutput -metadata VCFspec samples```

  This will print metadata info in the console.

### SimplifyVCF

- Help on SimplifyVCF
<pre>
$ python3 VcfSimplify.py SimplifyVCF -h

    ## VCF Simplify ## : Python application for parsing VCF files.
    Author: Bishwa K. Giri
    Contact: bkgiri@uncg.edu, kirannbishwa01@gmail.com

usage: VCF-Simplify SimplifyVCF [-h] -toType {haplotype,table} -inVCF INVCF
                                -outFile OUTFILE
                                [-outHeaderName OUTHEADERNAME]
                                [-GTbase GTBASE [GTBASE ...]] [-PG PG]
                                [-PI PI] [-includeUnphased {yes,no,0,1}]
                                [-samples SAMPLES [SAMPLES ...]]
                                [-preHeader PREHEADER [PREHEADER ...]]
                                [-infos INFOS [INFOS ...]]
                                [-formats FORMATS [FORMATS ...]]
                                [-mode {wide,long,0,1}]

optional arguments:
  -h, --help            show this help message and exit
                        
  -toType {haplotype,table}
                        Type of the output file.
                        
  -inVCF INVCF          Sorted vcf file.
                        
  -outFile OUTFILE      Name of the output file.
                        
  -outHeaderName OUTHEADERNAME
                        Write the VCF raw METADATA HEADER to a separate output file.
                        Default: no output.
                        
  -GTbase GTBASE [GTBASE ...]
                        Write the genotype (GT, PG etc.) field as IUPAC base code.
                        Default = GT:numeric (i.e write 'GT' bases as numeric)
                        Choices : [GT:numeric, GT:iupac, PG:iupac, .....]. 
                        Multiple option can be requested for each genotype fields.
                        

Additional arguments for "VCF To -> Haplotype":
  -PG PG                FORMAT tag representing the phased genotype.
                        Default: PG 
                        Note: 'GT' can be used if it contains the phased genotype.
                        
  -PI PI                FORMAT tag representing the unique phased haplotype block index.
                        Note: 'CHROM' can also be used as 'PI' if the VCF is phased chromosome or contig wide. 
                        
  -includeUnphased {yes,no,0,1}
                        include unphased variants (genotypes) in the haplotype output.
                        Default: no (0) (i.e do not write unphased variants)
                        

Additional arguments for "VCF To -> Table":
  -samples SAMPLES [SAMPLES ...]
                        SAMPLE of interest:
                          Space separated name of the samples or matching sample names.
                        Matching prefix, suffix or string in the names can be provided too.
                        Choices format:
                         [0, sample A, sample B, prefix:XXXsample, suffix:sampleXXX, match:XXX, all]
                        Multiple choices can be requested.
                        Note: 0 = ignore all the samples; Default = all
                        
  -preHeader PREHEADER [PREHEADER ...]
                        Space separated header fields before the 'INFO' field.
                        Choices:
                         [0, CHROM, POS, ID, REF, ALT, QUAL, FILTER, all]. 
                        Multiple choices can be requested.
                        Note: 0 = ignore all the pre-header-keys; Default = all.
                        
  -infos INFOS [INFOS ...]
                        Space separate INFO tags of interest.
                        Choices :
                         [0, AC, AF, AN, ..., all]. 
                        Multiple choices can be requested.
                        Note: 0 = ignore all the INFO tags; Default = all
                        
  -formats FORMATS [FORMATS ...]
                        Space separate FORMAT tags of interest.
                        Choices : [0, GT, PG, PI, ..., all]. 
                        Multiple choices can be requested.
                        Note: 0 = ignore all the pre-header-keys; Default = all.
                        
  -mode {wide,long,0,1}
                        Structure of the output table. Default = wide (0)

</pre>

- Example B-01 : VCF to Table
  
```python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -infos AF AN BaseQRankSum ClippingRankSum -formats GT PI PG GQ PL -preHeader CHROM POS REF ALT FILTER -mode wide -samples MA605 match:ms -GTbase GT:iupac PG:iupac -outHeaderName exampleOutput/vcf_header02.txt```

 This converts input_test.vcf file into simple_table.txt file using the arguments provided above.

- Example B-02 : VCF to Haplotype
  
```python3 VcfSimplify.py SimplifyVCF -toType haplotype -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_haplotype.txt -outHeaderName exampleOutput/vcf_header02.txt -PI PI -PG GT -GTbase GT:numeric -includeUnphased yes```

This converts input_test.vcf file into simple_haplotype.txt file using the arguments provided above. It also saves the header info intp vcf_header02.txt file.


### BuildVCF

- Help on BuildVCF
  
<pre>
 $ python3 VcfSimplify.py BuildVCF -h   

    ## VCF Simplify ## : Python application for parsing VCF files.
    Author: Bishwa K. Giri
    Contact: bkgiri@uncg.edu, kirannbishwa01@gmail.com

usage: VCF-Simplify BuildVCF [-h] -fromType {haplotype,table} -inFile INFILE
                             -outVCF OUTVCF -vcfHeader VCFHEADER
                             [-samples SAMPLES [SAMPLES ...]]
                             [-formats FORMATS [FORMATS ...]]
                             [-infos INFOS [INFOS ...]]
                             [-GTbase GTBASE [GTBASE ...]] [-haplotypeFormat]

optional arguments:
  -h, --help            show this help message and exit
                        
  -fromType {haplotype,table}
                        Type of the input file the VCF is being prepared from. 
                        
  -inFile INFILE        Sorted table or haplotype file.
                        Note: 
                          Haplotype file should be in the format created by 'phase-Stitcher' or 'phase-Extender'.
                          Table file should be in the format created by 'VCF-Simplify'
                        Only wide?? table is supported for now.
                        
  -outVCF OUTVCF        Name of the output VCF file.
                        
  -vcfHeader VCFHEADER  A custom VCF header to add to the VCF file.
                        The VCF header should not contain the line with #CHROM .... 
                        #CHROM ... line is auto populated while creating the VCF file.
                        

Additional arguments for "Table To VCF":
  -samples SAMPLES [SAMPLES ...]
                        Name of the samples -> space separated name of the samples that needs to be converted.
                        Default = allto VCF format
                        
  -formats FORMATS [FORMATS ...]
                        Name of the FORMAT tags to write -> space separated FORMAT tags name.
                        Default = all
                        
  -infos INFOS [INFOS ...]
                        Name of the INFO tags to write -> space separated INFO tags name.
                        Default = all
                        
  -GTbase GTBASE [GTBASE ...]
                        Suggest if the genotype fields (GT, PG etc.) are in IUPAC base code.
                        Default = GT:numeric (i.e assumes 'GT' bases are numeric)
                        Choices : [GT:numeric, GT:iupac, PG:iupac, .....]. 
                        Multiple option can be suggested for each genotype fields.
                        

Additional arguments for "Haplotype To VCF":
  -haplotypeFormat      report which format (numeric vs. iupac) the haplotype file is in.
                        Default = iupac
</pre>

- Example C01: From table to vcf:
  
  ```python3 VcfSimplify.py BuildVCF -fromType table -inFile exampleOutput/simple_table.txt -outVCF exampleOutput/tableToVcf.vcf -vcfHeader exampleOutput/vcf_header02.txt```

  This converts the input table (simple_table.txt) to vcf file(tableToVcf.vcf) and also saves header into vcf_header02.txt file.

- Example C02: From haplotype to vcf:
  ```python3 VcfSimplify.py BuildVCF -fromType haplotype -inFile exampleOutput/simple_haplotype.txt -outVCF exampleOutput/hapToVCF.vcf -vcfHeader exampleOutput/vcf_header.txt -haplotypeFormat iupac```

  This converts back the haplotype file created using simplify vcf into its original vcf file along with header file.


##  Upcoming features:

  - Ability to add `genotype bases` for fields other than ***"GT"*** .
  - Ability to handle ***symbolic alleles***.
  - Ability to :
    - prepare custom diploid genome.
    - prepare custom GTF, GFF files.
  - Extract gene sequence using ref genome and VCF files for phylogenetic analyses.


### Citation:
  &ensp; ***Giri, B.K, (2018). VCF-simplify: Tool to build and simplify VCF (variant call format) files.***

