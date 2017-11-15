#!/bin/bash

i=$1
CHROM=$2
STR=$3
SNP=$4

POS=`sed "${i}q;d" POS.txt`
let "start=$POS-100000"
let "end=$POS+100000"

bcftools query -f '[%TGT\t]\n' ${STR} -r $CHROM:$POS | sed 's/|/\t/g' | datamash transpose > str.${i}.txt
bcftools query -f '[%GT\t]\n' ${SNP} -r $CHROM:$start-$end | sed 's/|/\t/g' | datamash transpose > snp.${i}.txt
bcftools query -f '%ID\n' ${SNP} -r $CHROM:$start-$end > snp.id.${i}.txt
strid=`bcftools query -f '%ID\n' ${STR} -r $CHROM:$POS`

python corr.py ${strid} str.${i}.txt snp.${i}.txt snp.id.${i}.txt output.txt

rm str.${i}.txt snp.${i}.txt snp.id.${i}.txt
