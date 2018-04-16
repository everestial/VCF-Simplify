# VCF-simplify.v2
A python parser to simplify the vcf file into table like format.

There are several tools available to mainpulate and alter VCF file. But, a simple and comprehensive tool that can produce a most simple output required by emperical biologist is still amiss.

This tool takes in sorted vcf file and reports a simplified table output for `INFO` and `FORMAT` field for each `SAMPLE` of interest. With default state (minimal code) all the `INFO`, `FORMAT` for all the `SAMPLE` are simplified. Fields can be further narrowed down using very convenient and comprehensive scripts. **See the examples given below.**

The output table can be created in both "long" and "wide" format, which makes it suitable for **mining data by samples vs position** quite simple. The output can be further filtered downstream with awk and can be loaded onto R and used with tidyr, dplyr where different columns can be accessed by matching `names` or `pre,suf - fixes`.


# Prerequisites:
Python packages and modules:
- argparse (https://docs.python.org/3/library/argparse.html)
- cyvcf2 (https://github.com/brentp/cyvcf2/)
- Python3 (https://www.python.org/)


## Usage (**using the given input test data**): 

   ### Call for available options
    python3 vcf_simplify-v2.py --help    
<br>
    
   ### If no options are provided then all the INFO, FORMAT fields are reported from all the SAMPLE
    python3 vcf_simplify-v2.py --vcf input_test.vcf --out simplified_vcf.txt    
<br>
    
   ### Report wide output and "GT" as nucleotide bases
    python3 vcf_simplify-v2.py --vcf input_test.vcf --out simplified_vcf.txt --infos AF,AN,BaseQRankSum,ClippingRankSum --formats PI,GT,PG --pre_header CHROM,POS,REF,ALT,FILTER --mode wide --samples MA605,ms01e  --gtbase yes    
    
   ### Expected output
```
CHROM	POS	REF	ALT	FILTER	AF	AN	BaseQRankSum	ClippingRankSum	MA605_PI	MA605_GT	MA605_PG	ms01e_PI	ms01e_GT	ms01e_PG
2	15881018	G	A,C	PASS	1.0	8	-0.771	0.0	.	G/G	0/0	.	./.	./.
2	15881080	A	G	PASS	0.458	6	-0.732	0.0	.	A/A	0/0	.	./.	.
2	15881106	C	CA	PASS	0.042	6	0.253	0.0	.	C/C	0/0	.	./.	.
2	15881156	A	G	PASS	0.5	6	None	None	.	A/A	0/0	.	./.	.
2	15881224	T	G	PASS	0.036	12	1.75	0.0	.	T/T	0/0	.	./.	./.
2	15881229	C	G	PASS	0.308	10	None	None	.	C/C	0/0	.	./.	./.
```
<br>



   ### Report simiplified output (all available fields) for sample MA605,ms01e 
python3 vcf_simplify-v2.py --vcf F1.phased_variants.Final02.vcf --out simplified_vcf.txt --samples MA605,ms01e

   ### Expected output
```
CHROM	POS	ID	REF	ALT	QUAL	FILTER	AF	BaseQRankSum	ClippingRankSum	DP	DS	END	ExcessHet	FS	HaplotypeScore	InbreedingCoeff	MLEAC	MLEAF	MQ	MQRankSum	QD	RAW_MQ	ReadPosRankSum	SOR	set	SF	AC	AN	MA605_AD	MA605_DP	MA605_GQ	MA605_GT	MA605_MIN_DP	MA605_PGT	MA605_PID	MA605_PL	MA605_RGQ	MA605_SB	MA605_PG	MA605_PB	MA605_PI	MA605_PM	MA605_PW	MA605_PC	ms01e_AD	ms01e_DP	ms01e_GQ	ms01e_GT	ms01e_MIN_DP	ms01e_PGT	ms01e_PID	ms01e_PL	ms01e_RGQ	ms01e_SB	ms01e_PG	ms01e_PB	ms01e_PI	ms01e_PM	ms01e_PW	ms01e_PC
2	15881018	.	G	A,C	5082.45	PASS	1.0	-0.771	0.0	902	None	None	0.005	0.0	None	0.8	12,1	0.462,0.038	60.29	0.0	33.99	None	0.26	0.657	HignConfSNPs	0,1,2,3,4,5,6	2,0	8	3,0,0	3	9	0/0	None	None	None	0,9,112,9,112,112	None	None	0/0	.	.	.	0/0	.	0,0	0	.	./.	None	None	None	0,0,0,.,.,.	None	None	./.	.	.	.	./.	.
2	15881080	.	A	G	4336.44	PASS	0.458	-0.732	0.0	729	None	None	0.01	0.0	None	0.826	11	0.458	60.0	0.0	34.24	None	-0.414	0.496	HignConfSNPs	4,5,6	0	6	5,0	5	15	0/0	None	None	None	0,15,181	None	None	0/0	.	.	.	0/0	.	.	.	.	.	None	None	None	.	None	None	.	.	.	.	.	.
2	15881106	.	C	CA	33.32	PASS	0.042	0.253	0.0	654	None	None	3.01	0.0	None	-0.047	1	0.042	60.0	0.0	6.66	None	0.253	0.223	HignConfSNPs	4,5,6	0	6	6,0	6	18	0/0	None	None	None	0,18,206	None	None	0/0	.	.	.	0/0	.	.	.	.	.	None	None	None	.	None	None	.	.	.	.	.	.
```
<br>

   ### Report simplified output in "long" format
python3 vcf_simplify-v2.py --vcf F1.phased_variants.Final02.vcf --out simplified_vcf.txt --infos AF,AN,BaseQRankSum,ClippingRankSum --formats PI,GT,PG --pre_header CHROM,POS,REF,ALT,FILTER --mode long --samples MA605,ms01e
     
   ### Expected output
```
CHROM	POS	REF	ALT	FILTER	AF	AN	BaseQRankSum	ClippingRankSum	SAMPLE	PI	GT	PG
2	15881018	G	A,C	PASS	1.0	8	-0.771	0.0	MA605	.	0/0	0/0
2	15881018	G	A,C	PASS	1.0	8	-0.771	0.0	ms01e	.	./.	./.
2	15881080	A	G	PASS	0.458	6	-0.732	0.0	MA605	.	0/0	0/0
2	15881080	A	G	PASS	0.458	6	-0.732	0.0	ms01e	.	.	.
2	15881106	C	CA	PASS	0.042	6	0.253	0.0	MA605	.	0/0	0/0
2	15881106	C	CA	PASS	0.042	6	0.253	0.0	ms01e	.	.	.
2	15881156	A	G	PASS	0.5	6	None	None	MA605	.	0/0	0/0
2	15881156	A	G	PASS	0.5	6	None	None	ms01e	.	.	.
2	15881224	T	G	PASS	0.036	12	1.75	0.0	MA605	.	0/0	0/0
2	15881224	T	G	PASS	0.036	12	1.75	0.0	ms01e	.	./.	./.
2	15881229	C	G	PASS	0.308	10	None	None	MA605	.	0/0	0/0
2	15881229	C	G	PASS	0.308	10	None	None	ms01e	.	./.	./.
```
   
    
## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) Upcoming features:
  - Ability to add `genotype bases` for fields other than "GT".
  - Write the table back to a VCF file.
  
  
  **Citation:** ***Giri, B.K, (2018). VCF-simplify: A vcf simpification tool.***


