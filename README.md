# vcf_simplify-
A python parser to simplify the vcf file into table like format.

Basically this tool takes in sorted vcf file simplifies the vcf to a table like output by representing fields of interest from `INFO` and `FORMAT` field for each `SAMPLE` of interest.
By default all the `INFO`, `FORMAT` for all the `SAMPLE` are simplified


# Prerequisites:
Python packages and modules
- Python3 (https://www.python.org/)
- argparse
- cyvcf2 (https://github.com/brentp/cyvcf2/)

## Usage (**using the given test data in the example folder**): 

    python3 vcf_simplify-v1.py --vcf input_test.vcf --out simplified_vcf.txt --samples MA605,MA622,ms03g --infos AF,AN,BaseQRankSum,ClippingRankSum --formats GT,PI
    
    # Expected output
    CHR	POS	ID	REF	ALT	QUAL	FILTER	AF	AN	BaseQRankSum	MA605_GT	MA605_PI	MA622_GT	MA622_PI	ms03g_GT	ms03g_PI
    2	15881018	.	G	A,C	5082.45	PASS	1.0	8	-0.771	0/0	.	0/0	.	./.	.
    2	15881080	.	A	G	4336.44	PASS	0.458	6	-0.732	0/0	.	0/0	.	.	.
    2	15881106  ...........
    
    
    # if no options are provided then all the INFO, FORMAT field are represented from all the SAMPLE
    python3 vcf_simplify-v1.py --vcf input_test.vcf --out simplified_vcf.txt
    
