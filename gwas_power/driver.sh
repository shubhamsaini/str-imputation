#!/bin/bash

######
### ./driver.sh chromosome str.vcf.gz snp.vcf.gz numThreads window
#####

CHROM=$1
STR=$2
SNP=$3
NUMTHREADS=$4
WINDOW=$5

bcftools query -f '%POS\n' ${STR} > POS.txt
bcftools query -l ${STR} > samples.txt

#bcftools view ${SNP} --samples-file samples.txt -O z -o snp.chr${CHROM}.reorder.vcf.gz
#bcftools index snp.chr${CHROM}.reorder.vcf.gz
#SNP="snp.chr${CHROM}.reorder.vcf.gz"

nlines=`wc -l < POS.txt`

seq 1 ${nlines} | xargs -I% -P ${NUMTHREADS} ./find.corr.snp.sh % ${CHROM} ${STR} ${SNP} ${WINDOW}
