#!/bin/bash

# Check if a VCF file is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <input_vcf_file>"
  exit 1
fi

# Load required modules
module load StdEnv/2020
module load bcftools/1.11

# Get the input VCF file and extract the base name (without .vcf)
input_vcf="$1"
base_name=$(basename "$input_vcf" .vcf)

# Generate the output file name based on the VCF file name
output_file="${base_name}_homhetimputeddata.txt"

# Run the original command, replacing malaria_variant.vcf with the input file and outputting to the new file
paste <(bcftools view "$input_vcf" |\
    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' "$input_vcf" |\
    awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' "$input_vcf" |\
    awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' "$input_vcf" |\
    awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
    \
  <(bcftools view "$input_vcf" | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  <(bcftools view "$input_vcf" | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  <(bcftools view "$input_vcf" | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  | sed 's/,\t/\t/g' | sed 's/,$//g' > "$output_file"

echo "Output saved to $output_file"
