#!/bin/bash

######
### ./driver.sh 21 hipstr.chr21.phased.vcf.gz shapeit.chr21.with.ref.vcf.gz 4
#####

CHROM=$1
STR=$2
SNP=$3
NUMTHREADS=$4

#bcftools query -f '%POS\n' ${STR} > POS.txt

#bcftools query -l ${STR} > samples.txt
#bcftools view ${SNP} --samples-file samples.txt -O z -o shapeit.chr${CHROM}.reorder.vcf.gz
#bcftools index shapeit.chr${CHROM}.reorder.vcf.gz

nlines=`wc -l < POS.txt`

seq 1 5 | xargs -I% -P ${NUMTHREADS} ./find.corr.snp.sh % ${CHROM} ${STR} shapeit.chr${CHROM}.reorder.vcf.gz
