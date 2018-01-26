#!/bin/bash

######
### ./driver.sh chromosome str.vcf.gz snp.vcf.gz numThreads
#####

CHROM=$1
STR=$2
SNP=$3
NUMTHREADS=$4

bcftools query -f '%POS\n' ${STR} > POS.txt

#bcftools query -l ${STR} > samples.txt
#bcftools view ${SNP} --samples-file samples.txt -O z -o shapeit.chr${CHROM}.reorder.vcf.gz
#bcftools index shapeit.chr${CHROM}.reorder.vcf.gz
#SNP=shapeit.chr${CHROM}.reorder.vcf.gz

nlines=`wc -l < POS.txt`

seq 1 5 | xargs -I% -P ${NUMTHREADS} ./find.corr.snp.sh % ${CHROM} ${STR} ${SNP}
