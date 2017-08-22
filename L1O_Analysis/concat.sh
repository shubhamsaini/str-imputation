#!/bin/bash

### usage: 
### ./combine.str.snp.sh -s str.vcf.gz -p snp.vcf.gz -o out.vcf.gz
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
bcftools view $str -m2 -O z -o str.bial.vcf.gz
bcftools index str.bial.vcf.gz
bcftools query -l str.bial.vcf.gz > sampleID.txt
bcftools view $snp --samples-file sampleID.txt -O z -o snp.reorder.vcf.gz
bcftools index snp.reorder.vcf.gz
bcftools concat -a str.bial.vcf.gz snp.reorder.vcf.gz -O v -o temp.vcf
vcf-sort -c temp.vcf > temp2.vcf
bcftools norm -d any temp2.vcf -O z -o $out
bcftools index $out
rm temp.vcf temp2.vcf str.bial.vcf.gz str.bial.vcf.gz.csi snp.reorder.vcf.gz snp.reorder.vcf.gz.csi