# ARG-imputation

This repository contains the code necessary to perform ARG-imputation. For this purpose, it requires the following inputs:

1. An ancestral recombination graph (ARG) inferred by ARG needle in tskit format (check tskit documentation for details here: https://tskit.dev/tskit/docs/stable/introduction.html
2. A vcf file containing the variants that are going to be imputed

# ARG in tskit format
For detailsa about how to obtain an ARG from a VCF file, please check this repository: https://github.com/almejiaga/ARG_needle. Make sure you use the arg2tskit --arg_path ${arg_path} --ts_path {path_for_output} to get an ARG in the tskit format.

# Pre-processing VCF file
To performe the ARG-imputation, we need a set of carriers to use as input. I prepared a bash script that uses bcftools to get a list of homozygous for the reference allele, heterozygous and homozygous for the alternative allele from a VCF file. Run the script in the following way:

bash preprocessingvcf.sh $vcffile
The expected output:

```
CHR	POS	ID	REF	ALT	nHet	nHomAlt	nHomRef	HetSamples	HomSamplesAlt	HomSamplesRef
12	57763698	chr12_57763698_T_TGGGTGGG	T	TGGGTGGG	2	0	2170	Sample_1,Sample_2  Sample_3,sample_4  Sample_5		
```
# 



