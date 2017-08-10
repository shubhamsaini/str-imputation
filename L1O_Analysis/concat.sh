#!/bin/bash

## A script to concatenate STR and SNP VCF files
## None of the input files should have missing/unphased GT
## All markers should have atleast 2 alleles

### usage: 
### ./concat.sh -s str.vcf.gz -p snp.vcf.gz -o out.vcf.gz
###

while getopts s:p:o: option
do
        case "${option}"
        in
                s) str=${OPTARG};;
                p) snp=${OPTARG};;
                o) out=${OPTARG};;
        esac
done

bcftools index $str
bcftools index $snp
bcftools view $str -m2 -O z -o str.bial.vcf.gz
bcftools index str.bial.vcf.gz
bcftools concat -a str.bial.vcf.gz $snp -O v -o temp.vcf
vcf-sort -c temp.vcf > temp2.vcf
bcftools norm -d any temp2.vcf -O z -o $out
bcftools index $out
rm temp.vcf temp2.vcf str.bial.vcf.gz str.bial.vcf.gz.csi