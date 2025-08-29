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
# Performing ARG-imputation

After pre-processing the VCF file, it can be used as input along with the ARG file to perform genotype imputation. In addition, you will need a dictionary file between node IDs and your sample IDs (the default is the dictionary I have for the ARGs from CARTaGENE)

python ARG_imputation.py --chr 7 --variant $variant_id --hets $pre-processed_file  --ts $path_to_arg_.trees --dictionary $dictionary_correspondance --map_file $map_used_for_arg_inference

--variant: the variant ID, as stored in the VCF file

--hets: output file of the pre-processed vcf file

--ts: path to the ARG file for that chromosome in .trees format (ts kit format)

--dictionary: dictionary from the ARG leaves to the actual sample IDs (ARG needle by default renames from 0 to 2N samples)

--map_file: map used to generate the ARG, to perform start correction

This script will output a file containing two columns: sample IDs and posterior probability of being a carrier (from 0 to 1, check methods for details)
```
row ID_2	posterior
1	1	1.0
2	2	1.0
3	3	1.0
4	4	1.0
```
#Dependencies for running this repository

BCFTOOLS





