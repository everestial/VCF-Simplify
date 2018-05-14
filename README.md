# VCF-simplify.v2
A python parser to simplify the vcf file into table like format.
There are several tools available to mainpulate and alter VCF file.
But, a simple and comprehensive tool that can produce a most simple
output required by emperical biologist is still amiss.

  - **Convert VCF to TABLE**\
    This tool takes in sorted vcf file and reports a simplified table output
for `INFO` and `FORMAT` field for each `SAMPLE` of interest. With default
state (minimal code) all the `INFO`, `FORMAT` for all the `SAMPLE` are
simplified. Fields can be further narrowed down using very convenient
and comprehensive scripts.

  - **Convert TABLE to VCF**\
    It is also possible to convert the TABLE file into VCF. Controlled,
    workflows are included.

<br>

**Exclusively for phase-Stitcher and phase-Extender.**
  - **Convert VCF to Haplotype**\
    It is also possible to convert the TABLE file into VCF. Controlled,
    workflows are included.

  - **Convert Haplotype to VCF**\
    It is also possible to convert the TABLE file into VCF. Controlled,
    workflows are included.


## Prerequisites :
Required Python packages and modules
- [python3.x](https://www.python.org/)
- [argparse](https://docs.python.org/3/library/argparse.html)
- [cyvcf2](https://github.com/brentp/cyvcf2/)
- [numpy](http://www.numpy.org/)



## Usage :

### Call for help !

<pre>
$ python3 VCF-Simplify.py -h
  Checking required modules

  usage: VCF-Simplify [-h] {SimplifyVCF,BuildVCF} ...

  positional arguments:
    {SimplifyVCF,BuildVCF}
                          Choose one of the following method.
      SimplifyVCF         Simplify VCF : to a haplotype or a table file.
      BuildVCF            Create VCF : from a haplotype or a table file.

  optional arguments:
    -h, --help            show this help message and exit
</pre>

<br>
<br>

**You can see that there are two options available for conversion:**
  - **SimplifyVCF :** to convert from VCF ----> to TABLE and/or Haplotype file
  - **BuildVCF :** to convert from TABLE and/or Haplotype ---> to VCF
<br>


## SimplifyVCF

<pre>
$ python3 VCF-Simplify.py SimplifyVCF -h

Checking required modules

usage: VCF-Simplify SimplifyVCF [-h] -toType TOTYPE -inVCF INVCF -out OUT
                                [-keepHeader KEEPHEADER] [-PG PG] [-PI PI]
                                [-unphased UNPHASED] [-samples SAMPLES]
                                [-preHeader PREHEADER] [-infos INFOS]
                                [-formats FORMATS] [-mode MODE]
                                [-gtbase GTBASE]

optional arguments:
  -h, --help            show this help message and exit
  -toType TOTYPE        Type of the output file. Option: haplotype, table
  -inVCF INVCF          sorted vcf file
  -out OUT              name of the output file
  -keepHeader KEEPHEADER
                        Write the HEADER data to a separate output
                        file.Options: &apos;yes&apos; or &apos;no&apos;

VCF To Haplotype:
  -PG PG                FORMAT tag containing the phased genotype of the
                        SAMPLE. Only applicable on &apos;haplotype file output&apos;.
  -PI PI                FORMAT tag representing the unique index of RBphased
                        haplotype block in the SAMPLE. Only applicable on
                        &apos;haplotype file output&apos;. Note: &apos;CHROM&apos; can also be
                        used as PI if VCF is phased chromosome-wide.
  -unphased UNPHASED    include unphased variants in the output. Aavailable
                        options: yes, no

VCF To Table:
  -samples SAMPLES      SAMPLE of interest; write as comma separated names,
                        for e.g: &apos;sampleA,sampleB&apos; or &apos;all&apos;.
  -preHeader PREHEADER
                        Comma separated pre-header fields before the &apos;INFO&apos;
                        field in the input VCF file. Write as comma separated
                        fields, for e.g: &apos;CHR,POS,ID&apos; or &apos;all&apos;. Default:
                        &apos;all&apos;.
  -infos INFOS          INFO tags that are of interest; write as comma
                        separated tags; for e.g: &apos;AC,AF,AN&apos; or &apos;all&apos;.
  -formats FORMATS      FORMAT tags that are of interest; for e.g: &apos;GT,PG,PI&apos;
                        or &apos;all&apos;.
  -mode MODE            Structure of the output table.Options: wide(0),
                        long(1). Default: 0 .
  -gtbase GTBASE        write the GT field as IUPAC base code.Options: no(0),
                        yes(1). Default: 0 .
</pre>

<br>

### Example 01 (VCF to TABLE):
<pre>
python3 VCF-Simplify.py SimplifyVCF -toType table -inVCF input_test.vcf -out simple_table.txt -infos AF,AN,BaseQRankSum,ClippingRankSum -formats PI,GT,PG -preHeader CHROM,POS,REF,ALT,FILTER -mode wide -samples MA605,ms01e  -gtbase yes
</pre>

- Converts the "input_test.vcf" to a TABLE file.
- Uses samples MA605 and ms01e
- Uses FORMAT tags: PI,GT and PG
- Uses INFO tags: AF,AN,BaseQRankSum,ClippingRankSum
- Outputs TABLE in "wide" layout
- Outpts GT field as IUPAC base

**Output in Wide format**
<pre>
CHROM	POS	REF	ALT	FILTER	AF	AN	BaseQRankSum	ClippingRankSum	MA605:PI	MA605:GT	MA605:PG	ms01e:PI	ms01e:GT	ms01e:PG
CHROM	POS	REF	ALT	FILTER	AF	AN	BaseQRankSum	ClippingRankSum	MA605:PI	MA605:GT	MA605:PG	MA611:PI	MA611:GT	MA611:PG	MA622:PI	MA622:GT	MA622:PG
2	15881018	G	A,C	PASS	1.0	8	-0.771	0.0	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881080	A	G	PASS	0.458	6	-0.732	0.0	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881106	C	CA	PASS	0.042	6	0.253	0.0	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881156	A	G	PASS	0.5	6	None	None	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881224	T	G	PASS	0.036	12	1.75	0.0	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881229	C	G	PASS	0.308	10	None	None	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881230	GTT	G,GTTT	PASS	0.346,0.038	10	0.0	0.0	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881246	A	AT	PASS	0.333	8	None	None	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881256	A	T	PASS	0.333	8	None	None	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
2	15881266	T	G	PASS	0.286	12	None	None	.	0/0	0/0	.	0/0	0/0	.	0/0	0/0
</pre>

<br>

**Output in long format:**

    python3 VCF-Simplify.py SimplifyVCF -toType table -inVCF input_test.vcf -out simple_table.txt -infos AF,AN,BaseQRankSum,ClippingRankSum -formats PI,GT,PG -mode long -samples MA605,MA611,MA622

<pre>
CHROM	POS	ID	REF	ALT	QUAL	FILTER	AF	AN	BaseQRankSum	ClippingRankSum	SAMPLE	PI	GT	PG
2	15881018	.	G	A,C	5082.45	PASS	1.0	8	-0.771	0.0	MA605	.	0/0	0/0
2	15881018	.	G	A,C	5082.45	PASS	1.0	8	-0.771	0.0	MA611	.	0/0	0/0
2	15881018	.	G	A,C	5082.45	PASS	1.0	8	-0.771	0.0	MA622	.	0/0	0/0
2	15881080	.	A	G	4336.44	PASS	0.458	6	-0.732	0.0	MA605	.	0/0	0/0
2	15881080	.	A	G	4336.44	PASS	0.458	6	-0.732	0.0	MA611	.	0/0	0/0
2	15881080	.	A	G	4336.44	PASS	0.458	6	-0.732	0.0	MA622	.	0/0	0/0
2	15881106	.	C	CA	33.32	PASS	0.042	6	0.253	0.0	MA605	.	0/0	0/0
2	15881106	.	C	CA	33.32	PASS	0.042	6	0.253	0.0	MA611	.	0/0	0/0
2	15881106	.	C	CA	33.32	PASS	0.042	6	0.253	0.0	MA622	.	0/0	0/0
2	15881156	.	A	G	3595.22	PASS	0.5	6	None	None	MA605	.	0/0	0/0
2	15881156	.	A	G	3595.22	PASS	0.5	6	None	None	MA611	.	0/0	0/0
2	15881156	.	A	G	3595.22	PASS	0.5	6	None	None	MA622	.	0/0	0/0
</pre>

<br>

**Use minimal script for full term simplification:**\
(*Output not shown)

    python3 VCF-Simplify.py SimplifyVCF -to table -inVCF input_test.vcf -out simple_table.txt -keepHeader yes

- Simplified data for all the infos, formats for all the sample
- will output in wide format

<br>

**Include "-keepHeader " to store meta header of the VCF as separate file.**

    python3 VCF-Simplify.py SimplifyVCF -to table -inVCF input_test.vcf -out simple_table.txt -keepHeader yes

<br>

### Example 02 (VCF to Haplotype):

- Converts a VCF to a HAPLOTYPE file.
- The HAPLOTYPE file can be used downstream with tools:
  - phase-Extender
  - phase-Stitcher
- "PG" flag is used to indicate "phased-genotype" field in SAMPLE
- "PI" flag is used to indicate "haplotype-block" index
- Other FORMAT fields can be used with ***PG and PI***
- "CHROM" field can be included with ***PI*** flag, assuming the VCF is
phased chromosome wide.

<br>

    $ python3 VCF-Simplify.py SimplifyVCF -toType haplotype -inVCF input_test.vcf -out simple_haplotype.txt

<pre>
CHROM	POS	all-alleles	ms01e:PI	ms01e:PG_al	ms02g:PI	ms02g:PG_al	ms03g:PI	ms03g:PG_al	ms04h:PI	ms04h:PG_al	MA611:PI	MA611:PG_al	MA605:PI	MA605:PG_al	MA622:PI	MA622:PG_al
2	15881551	A,T	.	.	.	.	.	.	9	T|A	.	.	.	.	.	.
2	15881553	C,A	4	C|A	.	.	.	.	9	C|C	.	.	.	.	.	.
2	15881764	T,C	4	C|T	6	C|T	.	.	9	T|T	.	.	.	.	.	.
2	15881767	C,T	4	C|C	6	T|C	.	.	9	C|C	.	.	.	.	.	.
2	15881810	A,C	4	C|C	6	C|C	.	.	9	C|C	.	.	.	.	.	.
2	15881944	C,T	4	T|C	6	C|C	7	C|T	7	C|T	.	.	.	.	.	.
2	15881974	C,A	4	A|C	6	C|C	7	C|A	7	C|A	.	.	.	.	.	.
2	15881989	C,A	4	C|C	6	A|C	7	C|C	7	C|C	.	.	.	.	.	.
2	15882091	A,T	4	A|T	6	A|T	7	T|A	7	A|A	.	.	.	.	.	.
2	15882148	T,G	4	T|T	6	T|T	7	T|T	7	T|G	.	.	.	.	.	.
2	15882328	T,A	4	A|T	6	T|T	7	T|T	7	T|T	.	.	.	.	.	.
2	15882364	T,G	4	G|T	6	T|T	7	T|T	7	T|T	.	.	.	.	.	.
2	15882451	T,C	4	C|T	4	C|T	7	T|T	7	T|T	.	.	.	.	.	.
2	15882454	T,C	4	C|T	4	C|T	7	T|T	7	T|T	.	.	.	.	.	.
2	15882493	T,C	4	C|T	4	C|T	7	T|T	7	T|T	.	.	.	.	.	.
2	15882505	T,A	4	A|T	4	A|T	7	T|T	7	T|T	.	.	.	.	.	.
2	15882583	G,T	4	G|G	4	G|G	7	G|G	5	G|T	.	.	.	.	.	.
2	15882592	G,A	4	A|G	4	A|G	6	G|A	5	G|A	.	.	.	.	.	.
</pre>

<br>

\
**If your "GT" tag contains the phased genotype you can**

    python3 VCF-Simplify.py SimplifyVCF -to haplotype -inVCF input_test.vcf -out simple_table.txt -PG GT -PI PI

\
**If you want "CHROM" field as "-PI" index**

    python3 VCF-Simplify.py SimplifyVCF -to haplotype -inVCF input_test.vcf -out simple_table.txt -PG GT -PI CHROM

\
**Additionally to output unphased genotypes and store the header of the VCF**

    python3 VCF-Simplify.py SimplifyVCF -to haplotype -inVCF input_test.vcf -out simple_table.txt -PG GT -PI PI -unphased yes -keepHeader yes

<br>


## BuildVCF

<pre>$ python3 VCF-Simplify.py BuildVCF -h

Checking required modules

usage: VCF-Simplify BuildVCF [-h] -fromType FROMTYPE -inFile INFILE -outVCF
                             OUTVCF -vcfHeader VCFHEADER [-GTbase GTBASE]
                             [-samples SAMPLES] [-formats FORMATS]
                             [-infos INFOS]

optional arguments:
  -h, --help            show this help message and exit
  -fromType FROMTYPE    Type of the input file the VCF is being prepared from.
                        Options: haplotype, table
  -inFile INFILE        Sorted table or haplotype file.This haplotype file
                        should be obtained from phase-Stitcher, phase-
                        Extender. The table file should be in the format
                        output by &apos;VCF-Simplify&apos;; only long format table is
                        supported for now.
  -outVCF OUTVCF        Name of the output VCF file.
  -vcfHeader VCFHEADER  A custom VCF header to add to the VCF file. The VCF
                        header should not contain the line with #CHROM ....

Additional flags for &quot;Table To VCF&quot;:
  -GTbase GTBASE        Representation of the GT base is : numeric, IUPAC
  -samples SAMPLES      Name of the samples -&gt; comma separated name of the
                        samples that needs to be converted to VCF format
  -formats FORMATS      Name of the FORMAT tags to write -&gt; comma separated
                        FORMAT tags name.
  -infos INFOS          Name of the INFO tags to write -&gt; comma separated INFO
                        tags name.
</pre>


<br>

### Example 03 (TABLE to VCF):
*Note:
- requires "vcf header" from other VCF or custom VCF, with no `#CHROM` line.
- the type of the data in "GT" should be indicated.

<pre>
python3 VCF-Simplify.py BuildVCF -fromType table -inFile simple_table.txt -vcfHeader vcf_header.txt -outVCF table_toVCF.vcf -GTbase numeric
</pre>

<br>


### Example 04 (Haplotype to VCF):
*Note: This run with a minimal script.

<pre>
python3 VCF-Simplify.py BuildVCF -fromType haplotype -inFile simple_haplotype.txt -vcfHeader vcf_header.txt -outVCF haplotype_toVCF.vcf
</pre>

<br>

## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) Upcoming features:
  - Ability to add `genotype bases` for fields other than "GT".
  - Ability to handle symbolic alleles.
  - Ability to :
    - prepare custom diploid genome.
    - prepare custom GTF, GFF files.
  - Extract gene sequence using ref genome and VCF files for phylogenetic analyses.

<br>
<br>
  
  **Citation:**\
  &ensp; ***Giri, B.K, (2018). VCF-simplify: Tool to build and simplify VCF (variant call format) files.***


