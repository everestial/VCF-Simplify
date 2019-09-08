# VCF-Simplify &ensp;&ensp;&ensp; v3.0.1

**A python parser to parse metadata and records information from VCF file.** 

There are several tools available to manipulate and parse VCF file. But, a simple and comprehensive tool that can extract data in a most efficient (minimal code) way is still amiss. VCF simplify is aimed for users from non-programming background, empirical biologist and team who lack bioinformatics support, however can be used by users at any level to extract VCF data in a most efficient way. This application takes the approach of using minimal scripts to generate maximal output.

**This tool performs following three main tasks:**

1. **ViewVCF** : Display and extract metadata information from VCF files.
2. **SimplifyVCF** : Convert VCF files to table or haplotype format.
3. **BuildVCF** : Convert the table or haplotype format file back to VCF.

**Note:** **haplotype** format data is exclusively used in [phase-Stitcher](https://github.com/everestial/pHASE-Stitcher) and [phase-Extender](https://github.com/everestial/phase-Extender) for phasing ReadBackPhased haplotypes.

This application is not designed to merge, compare or compute statistics from VCF files. Tools like [vcflib](https://github.com/vcflib/vcflib) and [bcftools](https://github.com/samtools/bcftools) are the ones better suited for those purpose.

## Table of contents

- [VCF-Simplify &ensp;&ensp;&ensp; v3.0.1](#vcf-simplify-enspenspensp-v301)
  - [Table of contents](#table-of-contents)
  - [Tutorial](#tutorial)
    - [Prerequisites](#prerequisites)
    - [Installation and setup](#installation-and-setup)
  - [Usage](#usage)
    - [ViewVCF](#viewvcf)
    - [SimplifyVCF](#simplifyvcf)
    - [BuildVCF](#buildvcf)
  - [Upcoming features:](#upcoming-features)
    - [Citation:](#citation)

## Tutorial

### Prerequisites

**VCF Simplify** is written in **python3** and only uses standard built-in modules. So, all you need is **`python3`** installed on your system (windows, mac, ubuntu) to run this code locally. If you do not have python installed, you can install it from **[here](https://www.python.org/downloads)**. For linux; you can get latest **python3** by:

```bash
sudo apt-get install python3
```

<br>

### Installation  and setup

1. **Clone this repo**

``` bash
git clone https://github.com/everestial/VCF-Simplify
cd VCF-Simplify

# call the "VCF-Simplify" application
python3 VcfSimplify.py -h
```

2. **Cythonize** (Optional but helpful for faster performance) 

   - **Create and activate a virtual environment**

     ```bash
     # create virtual environment named "myenv"
     python3 -m venv myenv       # rename environment as required
     
     # activate the environment 
     source myenv/bin/activate   # for linux
     myenv\Scripts\activate      # for windows
     ```

   - **Install `Cython` library**

     ```bash
     (myenv)$ pip install Cython
     ```

   - **Cythonize (so that app runs faster)**

     ```bash
     # cythoize using the 'setup.py' file in "VCF Simplify"
     (myenv)$ python3 setup.py build_ext --inplace
     ```

     This makes .so files in your build directory and it will be called when importing modules rather than actual modules.
     

3. **Call `VCFSimplify`** 

```bash
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
```
**Three available functions are shown:**  

- ViewVCF, SimplifyVCF, BuildVCF.
- Each function can be further expanded as shown in the **"Usage"** below.

<br>

## Usage

### ViewVCF

- Help on **ViewVCF**:

```html
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
```

**examples:** 

- Write metadata from VCF file as table, json, dict
  
  ```bash
  python3 VcfSimplify.py ViewVCF -inVCF exampleInput/input_test.vcf -outFile exampleOutput/testOutput -outType table json dict  
  ```
This will write metadata information from vcf headers into three format (table, json and dict) as testOutput.table and so on. 
  
  **Note:** Since no `-metadata` flag is invoked, at default state all the available metadata are extracted.

<br>

- Use `-metadata` flag to extract only specific metadata of interest (on console)

  ```bash
  python3 VcfSimplify.py ViewVCF -inVCF exampleInput/input_test.vcf -metadata VCFspec samples
  ```
  Since no output filename is given, this will print metadata info in the console. 

<br>

### SimplifyVCF

- Help on **SimplifyVCF**
```html
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

```
<br>

**examples:** VCF to Table

- Convert specific metadata into simplified tsv format

```python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt -infos AF AN BaseQRankSum ClippingRankSum -formats GT PI PG GQ PL -preHeader CHROM POS REF ALT FILTER -mode wide -samples MA605 match:ms -GTbase GT:iupac PG:iupac -outHeaderName exampleOutput/vcf_header02.txt```

 This converts [***input_test.vcf***](exampleInput/input_test.vcf) into [***simple_table.txt***](exampleOutput/simple_table.txt) file using the arguments provided above.

 <br>

 - Convert all the VCF into simplified tsv format 

 ```python3 VcfSimplify.py SimplifyVCF -toType table -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_table.txt```

 If no specific argument flag is raised, then all the tags from INFO and FORMAT (for all the samples) are simplified as tsv. 


<br>

**example:** VCF to Haplotype

```python3 VcfSimplify.py SimplifyVCF -toType haplotype -inVCF exampleInput/input_test.vcf -outFile exampleOutput/simple_haplotype.txt -outHeaderName exampleOutput/vcf_header02.txt -PI PI -PG GT -GTbase GT:numeric -includeUnphased yes```

This converts [***input_test.vcf***](exampleInput/input_test.vcf) into [***simple_haplotype.txt***](exampleOutput/simple_haplotype.txt). 

The raised arugments do the following:  

  - `PI` from the FORMAT is used as unique index for the phased block for each sample. 
  - phased genotypes are extracted from `GT`.
  - phased genotypes are represented as numerical bases. 
  - uphased genotypes are also written to haplotype file. 
  - header info is written as [***vcf_header02.txt***](exampleOutput/vcf_header02.txt).

<br>

### BuildVCF

- Help on BuildVCF
  
```
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
```

**example:** Table to VCF

  ```python3 VcfSimplify.py BuildVCF -fromType table -inFile exampleOutput/simple_table.txt -outVCF exampleOutput/tableToVcf.vcf -vcfHeader exampleOutput/vcf_header02.txt```

  This converts the input table [***simple_table.txt***](exampleOutput/simple_table.txt) to vcf file [***tableToVcf.vcf***](exampleOutput/tableToVCF.vcf) and adds [***vcf_header02.txt***](exampleOutput/vcf_header02.txt) as metadata header.

**example:** haplotype to VCF

  ```python3 VcfSimplify.py BuildVCF -fromType haplotype -inFile exampleOutput/simple_haplotype.txt -outVCF exampleOutput/hapToVCF.vcf -vcfHeader exampleOutput/vcf_header.txt -haplotypeFormat iupac```

  This converts the haplotype file [***simple_haplotype.txt***](exampleOutput/simple_haplotype.txt) created using simplify vcf into its original vcf file [***hapToVCF.vcf***](exampleOutput/hapToVCF.vcf) along with header file.

#### Several other useful examples/scripts can be found in `TestScripts.sh` file. 

## Upcoming features:

  - Ability to add `genotype bases` for fields other than ***"GT"*** .
  - Ability to handle ***symbolic alleles***.
  - Ability to :
    - prepare custom diploid genome.
    - prepare custom GTF, GFF files.
  - Extract gene sequence using ref genome and VCF files for phylogenetic analyses.

### Citation:

  &ensp; ***Giri, B.K, (2018). VCF-simplify: Tool to build and simplify VCF (variant call format) files.***
